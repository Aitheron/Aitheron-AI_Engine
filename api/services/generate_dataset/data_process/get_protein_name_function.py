import requests
import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from tenacity import retry, wait_exponential, stop_after_attempt

UNIPROT_BASE = "https://rest.uniprot.org"
HEADERS_JSON = {"Accept": "application/json"}

DB_ALIAS = {
    "ENSEMBL_PROTEIN_ID": "Ensembl_Protein",
    "ENSEMBL_TRANSCRIPT_ID": "Ensembl_Transcript",
    "ENSEMBL_GENE_ID": "Ensembl",
    "UniProtKB": "UniProtKB",
    "UNIPROTKB": "UniProtKB",
    "Ensembl_Protein": "Ensembl_Protein",
    "Ensembl_Transcript": "Ensembl_Transcript",
    "Ensembl": "Ensembl",
}

@retry(stop=stop_after_attempt(5), wait=wait_exponential(min=1, max=20))
def run_uniprot_id_mapping(from_db: str, to_db: str, ids: list[str]) -> str:

    from_db = DB_ALIAS.get(from_db, from_db)
    to_db = DB_ALIAS.get(to_db, to_db)

    payload_ids = " ".join([i for i in ids if isinstance(i, str) and i.strip()])
    if not payload_ids:
        raise ValueError("Lista de IDs vazia após saneamento.")

    r = requests.post(
        f"{UNIPROT_BASE}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": payload_ids},
        headers=HEADERS_JSON,
        timeout=60,
    )
    r.raise_for_status()
    j = r.json()
    job_id = j.get("jobId")
    if not job_id:
        raise RuntimeError(f"Resposta sem jobId: {j}")
    return job_id

@retry(stop=stop_after_attempt(10), wait=wait_exponential(min=1, max=30))
def poll_uniprot_job(job_id: str) -> dict:
    r = requests.get(f"{UNIPROT_BASE}/idmapping/status/{job_id}", headers=HEADERS_JSON, timeout=60)
    r.raise_for_status()
    return r.json()

@retry(stop=stop_after_attempt(5), wait=wait_exponential(min=1, max=20))
def fetch_uniprot_id_mapping_results(job_id: str) -> list[dict]:
    results: list[dict] = []
    next_url = f"{UNIPROT_BASE}/idmapping/results/{job_id}?format=json&size=500"
    while next_url:
        r = requests.get(next_url, timeout=120, headers=HEADERS_JSON)
        r.raise_for_status()
        j = r.json()
        results.extend(j.get("results", []))
        next_link = None
        for link in j.get("links", []):
            if link.get("rel") == "next":
                next_link = link.get("href")
                break
        next_url = next_link
    return results

def _safe_taxon_id(d: dict):
    try:
        return d.get("organism", {}).get("taxonId")
    except Exception:
        return None

def _normalize_target(tgt) -> dict | None:

    if isinstance(tgt, str):
        return {"primaryAccession": tgt}

    if isinstance(tgt, dict):
        if "primaryAccession" not in tgt:
            acc = tgt.get("uniProtkbId") or tgt.get("accession") or tgt.get("primaryAccession")
            if isinstance(acc, str):
                tgt["primaryAccession"] = acc
        return tgt if tgt.get("primaryAccession") else None

    return None

def _choose_best_target(targets: list[dict], restrict_taxon: int | None) -> dict | None:

    if not targets:
        return None

    if restrict_taxon is not None:
        filtered = [t for t in targets if _safe_taxon_id(t) == restrict_taxon]
    else:
        filtered = list(targets)

    candidates = filtered or targets
    reviewed = [t for t in candidates if t.get("reviewed") is True]
    return (reviewed or candidates)[0] if (reviewed or candidates) else None

def _map_chunk_ids_to_uniprot(from_db: str, ids: list[str], restrict_taxon: int | None) -> dict[str, str]:

    job_id = run_uniprot_id_mapping(from_db, "UniProtKB", ids)

    while True:
        status = poll_uniprot_job(job_id)
        js = status.get("jobStatus")
        if js in (None, "FINISHED"):
            break
        if js == "FAILED":
            raise RuntimeError("UniProt ID mapping failed")
        time.sleep(2)

    raw_results = fetch_uniprot_id_mapping_results(job_id)

    grouped: dict[str, list[dict]] = {}
    for item in raw_results:
        src = item.get("from")
        tgt = _normalize_target(item.get("to"))
        if not src or not tgt or not isinstance(tgt, dict) or not tgt.get("primaryAccession"):
            continue
        grouped.setdefault(src, []).append(tgt)

    chosen: dict[str, str] = {}
    for src, targets in grouped.items():
        winner = _choose_best_target(targets, restrict_taxon)
        if winner and winner.get("primaryAccession"):
            chosen[src] = winner["primaryAccession"]

    return chosen

def map_ids_to_uniprot(from_db: str, ids: list[str], restrict_taxon: int | None = 9606, chunk_size: int = 300) -> dict[str, str]:

    clean_ids = [i for i in ids if isinstance(i, str) and i.strip()]
    if not clean_ids:
        return {}

    seen = set()
    unique_ids = []
    for i in clean_ids:
        if i not in seen:
            seen.add(i)
            unique_ids.append(i)

    out: dict[str, str] = {}
    for i in range(0, len(unique_ids), chunk_size):
        chunk = unique_ids[i : i + chunk_size]
        mapped = _map_chunk_ids_to_uniprot(from_db, chunk, restrict_taxon)
        out.update(mapped)

    return out

@retry(stop=stop_after_attempt(5), wait=wait_exponential(min=1, max=20))
def fetch_uniprot_entry_json(accession: str) -> dict | None:
    r = requests.get(f"{UNIPROT_BASE}/uniprotkb/{accession}.json", headers=HEADERS_JSON, timeout=60)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()

def extract_name_function_features(entry: dict):

    protein_description = pd.NA
    try:
        protein_description = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except Exception:
        protein_description = entry.get("uniProtkbId") or pd.NA

    function_texts = []
    for c in entry.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            for t in c.get("texts", []):
                v = t.get("value")
                if v:
                    function_texts.append(v)
    protein_function = " | ".join(function_texts) if function_texts else pd.NA

    feature_summaries = []
    for f in entry.get("features", []):
        ftype = f.get("type")
        desc = f.get("description") or ""
        loc = f.get("location", {})
        begin = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if begin and end:
            feature_summaries.append(f"{ftype}: {desc} [{begin}-{end}]".strip())
        else:
            feature_summaries.append(f"{ftype}: {desc}".strip())
    protein_features = " ; ".join(feature_summaries) if feature_summaries else pd.NA

    return protein_description, protein_function, protein_features

def fetch_details_for_accessions(accessions: list[str], max_workers: int = 20) -> dict[str, dict]:
    unique_accessions = sorted({a for a in accessions if isinstance(a, str) and a.strip()})
    results: dict[str, dict] = {}

    def task(acc: str):
        entry = fetch_uniprot_entry_json(acc)
        if not entry:
            return acc, None
        protein_desc, func, feats = extract_name_function_features(entry)
        return acc, {"ProteinDesc": protein_desc, "ProteinFunction": func, "ProteinFeatures": feats}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(task, acc) for acc in unique_accessions]
        for fut in as_completed(futures):
            acc, info = fut.result()
            if info:
                results[acc] = info
    return results

def enrich_with_uniprot(df: pd.DataFrame, output_path: str | None = None, ensp_col: str = "ENSP", enst_col: str = "ENST", max_workers: int = 20) -> pd.DataFrame:
    df_out = df.copy()

    if "ProteinName" not in df_out.columns:
        df_out["ProteinName"] = pd.NA
    if "ProteinFunction" not in df_out.columns:
        df_out["ProteinFunction"] = pd.NA
    if "ProteinFeatures" not in df_out.columns:
        df_out["ProteinFeatures"] = pd.NA
    if "ProteinDesc" not in df_out.columns:
        df_out["ProteinDesc"] = pd.NA

    ensp_values = df_out[ensp_col].dropna().astype(str).str.strip() if ensp_col in df_out.columns else pd.Series(dtype=str)
    enst_values = df_out[enst_col].dropna().astype(str).str.strip() if enst_col in df_out.columns else pd.Series(dtype=str)

    unique_ensps = sorted({v for v in ensp_values if v})
    unique_ensts = sorted({v for v in enst_values if v})

    ensp_to_uniprot = map_ids_to_uniprot("Ensembl_Protein", unique_ensps) if unique_ensps else {}
    enst_to_uniprot = map_ids_to_uniprot("Ensembl_Transcript", unique_ensts) if unique_ensts else {}

    row_to_accession: dict[int, str] = {}
    for idx, row in df_out.iterrows():
        accession = None
        ensp_val = str(row.get(ensp_col)) if ensp_col in df_out.columns else None
        enst_val = str(row.get(enst_col)) if enst_col in df_out.columns else None

        if ensp_val and ensp_val in ensp_to_uniprot:
            accession = ensp_to_uniprot[ensp_val]
        elif enst_val and enst_val in enst_to_uniprot:
            accession = enst_to_uniprot[enst_val]

        if accession:
            row_to_accession[idx] = accession

    accession_list = list(row_to_accession.values())
    accession_details = fetch_details_for_accessions(accession_list, max_workers=max_workers)

    name_series = {}
    func_series = {}
    feat_series = {}
    protein_desc_series = {}
    for idx, acc in row_to_accession.items():
        info = accession_details.get(acc)
        if not info:
            continue
        name_series[idx] = acc
        protein_desc_series[idx] = info["ProteinDesc"]
        func_series[idx] = info["ProteinFunction"]
        feat_series[idx] = info["ProteinFeatures"]

    if name_series:
        df_out.loc[list(name_series.keys()), "ProteinName"] = pd.Series(name_series)
    if func_series:
        df_out.loc[list(func_series.keys()), "ProteinFunction"] = pd.Series(func_series)
    if feat_series:
        df_out.loc[list(feat_series.keys()), "ProteinFeatures"] = pd.Series(feat_series)
    if protein_desc_series:
        df_out.loc[list(protein_desc_series.keys()), "ProteinDesc"] = pd.Series(protein_desc_series)

    if output_path:
        df_out.to_csv(output_path, index=False)

    return df_out

if __name__ == "__main__":
    input_path = "./files/clinvar_BRCA1_GRCh38_with_alleles_and_proteins1.csv"
    df_in = pd.read_csv(input_path, dtype=str)

    df_out = enrich_with_uniprot(
        df_in,
        output_path="./files/clinvar_BRCA1_GRCh38_with_uniprot.csv",
        ensp_col="ENSP",
        enst_col="ENST",
        max_workers=20,
    )
    
    len(df_out)

    df2 = pd.read_csv("./files/clinvar_BRCA2_GRCh38_with_alleles_and_proteins1.csv", dtype=str)
    df2_out = enrich_with_uniprot(
        df2,
        output_path="./files/clinvar_BRCA2_GRCh38_with_uniprot.csv",
        ensp_col="ENSP",
        enst_col="ENST",
        max_workers=20,
    )
