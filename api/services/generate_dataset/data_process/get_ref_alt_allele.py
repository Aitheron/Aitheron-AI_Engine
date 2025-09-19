import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tenacity import retry, wait_exponential, stop_after_attempt

ENSEMBL_VEP_ID_URL = "https://rest.ensembl.org/vep/human/id/{rsid}"
ENSEMBL_LOOKUP_ENST_URL = "https://rest.ensembl.org/lookup/id/{enst}?expand=1"
ENSEMBL_SEQUENCE_ID_URL = "https://rest.ensembl.org/sequence/id/{item_id}?type={seqtype}"

HEADERS_JSON = {"Accept": "application/json"}
HEADERS_FASTA = {"Accept": "text/x-fasta"}

def normalize_rs(value) -> str | None:
    if value is None:
        return None
    base_str = str(value).strip()
    if base_str.lower() in {"-1", "nan", "none", ""}:
        return None
    if base_str.lower().startswith("rs"):
        base = base_str[2:]
    else:
        base = base_str
    return f"rs{base}"

def _has_protein_fields(t: dict) -> bool:
    return any(t.get(k) is not None for k in ("protein_start", "protein_end", "amino_acids"))

def _impact_rank(val: str | None) -> int:
    order = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    return order.get((val or "").upper(), 0)

def _choose_transcript(transcripts: list[dict]) -> tuple[dict | None, str]:

    if not transcripts:
        return None, "fallback"

    # 1) protein_fields
    with_protein = [t for t in transcripts if _has_protein_fields(t)]
    if with_protein:
        candidates = with_protein
        basis = "protein_fields"
    else:
        candidates = transcripts
        basis = None

    # 2) impact
    max_imp = max((_impact_rank(t.get("impact")) for t in candidates), default=0)
    by_impact = [t for t in candidates if _impact_rank(t.get("impact")) == max_imp]
    if len(by_impact) == 1:
        return by_impact[0], (basis or "impact")

    # 3) menor distance
    with_dist = [t for t in by_impact if isinstance(t.get("distance"), int)]
    if with_dist:
        chosen = min(with_dist, key=lambda x: x.get("distance", 10**9))
        return chosen, (basis or "distance")

    # 4) fallback
    return by_impact[0], (basis or "fallback")

NONCODING_UPDOWN_TERMS = {"upstream_gene_variant", "downstream_gene_variant"}
CODING_TERMS = {
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "start_lost", "protein_altering_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
}

def _is_updown_with_distance(t: dict) -> bool:

    terms = {s.lower() for s in (t.get("consequence_terms") or [])}
    if not terms:
        return False
    only_updown = terms.issubset(NONCODING_UPDOWN_TERMS)
    has_distance = isinstance(t.get("distance"), int)
    return only_updown and has_distance

def _is_coding_tc(t: dict) -> bool:
    """True se o item tem evidência de efeito codificante."""
    terms = {s.lower() for s in (t.get("consequence_terms") or [])}
    if terms & CODING_TERMS:
        return True
    if _has_protein_fields(t):
        return True
    if any(t.get(k) is not None for k in ("cds_start", "cds_end", "protein_start", "protein_end")):
        return True
    return False

@retry(stop=stop_after_attempt(20), wait=wait_exponential(min=1, max=60))
def fetch_vep_record(rsid: str, session: requests.Session, timeout=30) -> dict | None:
    url = ENSEMBL_VEP_ID_URL.format(rsid=rsid)
    r = session.get(url, headers=HEADERS_JSON, timeout=timeout)
    if r.status_code != 200:
        raise Exception(f"HTTP {r.status_code} on VEP")
    try:
        data = r.json()
    except Exception:
        raise Exception("JSON parse error on VEP")
    if not isinstance(data, list) or not data:
        raise Exception("Empty VEP response")
    return data[0]

@retry(stop=stop_after_attempt(10), wait=wait_exponential(min=1, max=30))
def fetch_translation_id_for_enst(enst: str, session: requests.Session, timeout=30) -> str | None:

    url = ENSEMBL_LOOKUP_ENST_URL.format(enst=enst)
    r = session.get(url, headers=HEADERS_JSON, timeout=timeout)
    if r.status_code != 200:
        raise Exception(f"HTTP {r.status_code} on lookup ENST")
    j = r.json()
    tr = j.get("Translation") or j.get("translation")
    if isinstance(tr, dict) and tr.get("id"):
        return tr["id"]
    return None

@retry(stop=stop_after_attempt(10), wait=wait_exponential(min=1, max=30))
def fetch_fasta(item_id: str, seqtype: str, session: requests.Session, timeout=30) -> str | None:

    if not item_id:
        return None
    url = ENSEMBL_SEQUENCE_ID_URL.format(item_id=item_id, seqtype=seqtype)
    r = session.get(url, headers=HEADERS_FASTA, timeout=timeout)
    if r.status_code != 200:
        raise Exception(f"HTTP {r.status_code} on sequence {seqtype}")
    text = r.text
    if not text or not text.startswith(">"):
        raise Exception("Empty FASTA")
    return text

def extract_ref_alt_from_rec(rec: dict, rsid: str) -> tuple[str | None, str | None]:

    allele_string = rec.get("allele_string")
    if not allele_string:
        coloc = rec.get("colocated_variants") or []
        for c in coloc:
            if c.get("id", "").lower() == rsid.lower() and c.get("allele_string"):
                allele_string = c["allele_string"]
                break
        if not allele_string and coloc:
            allele_string = coloc[0].get("allele_string")

    if not allele_string or "/" not in allele_string:
        return None, None

    parts = allele_string.split("/")
    if len(parts) < 2:
        return None, None

    ref = parts[0]
    alts_str = ",".join(parts[1:])
    return ref, alts_str

def process_item(args):
    idx, raw_val, session = args

    rsid = normalize_rs(raw_val)
    if not rsid:
        return idx, None, None, {}

    rec = fetch_vep_record(rsid, session)

    ref, alts = extract_ref_alt_from_rec(rec, rsid)

    tcs = rec.get("transcript_consequences") or []
    pc = [t for t in tcs if (t.get("biotype") == "protein_coding")]

    # Split: remover itens que são SÓ upstream/downstream com distance
    pc_non_updown = [t for t in pc if not _is_updown_with_distance(t)]
    pc_updown     = [t for t in pc if _is_updown_with_distance(t)]

    enrich = {
        "IsCoding": False,
        "TranscriptSelectionBasis": pd.NA,
        "ENST": pd.NA,
        "ENSP": pd.NA,
        "ConsequenceTerms": pd.NA,
        "Impact": pd.NA,
        "VariantAllele": pd.NA,
        "ProteinStart": pd.NA,
        "ProteinEnd": pd.NA,
        "CDNAStart": pd.NA,
        "CDNAEnd": pd.NA,
        "CDSStart": pd.NA,
        "CDSEnd": pd.NA,
        "Codons": pd.NA,
        "AminoAcids": pd.NA,
        "Strand": pd.NA,
        "FASTA_CDS": pd.NA,
        "FASTA_Protein": pd.NA,
    }

    chosen = None
    basis = "fallback"
    coding = False

    if pc_non_updown:
        chosen, basis = _choose_transcript(pc_non_updown)
        coding = bool(chosen and _is_coding_tc(chosen))
    elif pc_updown:
        chosen = pc_updown[0]
        basis = "non_coding_up/down_with_distance"
        coding = False
    else:
        # sem protein_coding — mantém não-codificante
        return idx, ref, alts, enrich

    # Se não-codificante: só marque o basis e retorne (demais colunas vazias)
    if not coding:
        enrich["IsCoding"] = False
        enrich["TranscriptSelectionBasis"] = basis
        return idx, ref, alts, enrich

    # Codificante: preencher campos e buscar ENSP/FASTA
    enst = chosen.get("transcript_id")
    ensp = chosen.get("protein_id")

    enrich.update({
        "IsCoding": True,
        "TranscriptSelectionBasis": basis,
        "ENST": enst or pd.NA,
        "ENSP": ensp or pd.NA,
        "ConsequenceTerms": "|".join(chosen.get("consequence_terms", [])) or pd.NA,
        "Impact": chosen.get("impact") or pd.NA,
        "VariantAllele": chosen.get("variant_allele") or pd.NA,
        "ProteinStart": chosen.get("protein_start") if chosen.get("protein_start") is not None else pd.NA,
        "ProteinEnd": chosen.get("protein_end") if chosen.get("protein_end") is not None else pd.NA,
        "CDNAStart": chosen.get("cdna_start") if chosen.get("cdna_start") is not None else pd.NA,
        "CDNAEnd": chosen.get("cdna_end") if chosen.get("cdna_end") is not None else pd.NA,
        "CDSStart": chosen.get("cds_start") if chosen.get("cds_start") is not None else pd.NA,
        "CDSEnd": chosen.get("cds_end") if chosen.get("cds_end") is not None else pd.NA,
        "Codons": chosen.get("codons") or pd.NA,
        "AminoAcids": chosen.get("amino_acids") or pd.NA,
        "Strand": chosen.get("strand") if chosen.get("strand") is not None else pd.NA
    })

    if enst:
        try:
            if not ensp:
                ensp = fetch_translation_id_for_enst(enst, session)
                if ensp:
                    enrich["ENSP"] = ensp
        except Exception:
            pass
        try:
            fasta_cds = fetch_fasta(enst, "cds", session)
            if fasta_cds:
                enrich["FASTA_CDS"] = fasta_cds
        except Exception:
            pass

    if ensp:
        try:
            fasta_prot = fetch_fasta(ensp, "protein", session)
            if fasta_prot:
                enrich["FASTA_Protein"] = fasta_prot
        except Exception:
            pass

    return idx, ref, alts, enrich

def update_df_inplace_with_vep(
    df: pd.DataFrame,
    output_path: str,
    rs_col: str = "RS# (dbSNP)",
    ref_col: str = "ReferenceAllele",
    alt_col: str = "AlternateAllele",
    max_workers: int = 100,
):
    df_out = df.copy()
    session = requests.Session()

    if ref_col not in df_out.columns:
        df_out[ref_col] = pd.NA
    if alt_col not in df_out.columns:
        df_out[alt_col] = pd.NA

    new_cols = [
        "IsCoding", "TranscriptSelectionBasis", "ENST", "ENSP",
        "ConsequenceTerms", "Impact", "VariantAllele", "ProteinStart",
        "ProteinEnd", "CDNAStart", "CDNAEnd", "CDSStart", "CDSEnd",
        "Codons", "AminoAcids", "Strand",
        "FASTA_CDS", "FASTA_Protein",
    ]
    for column in new_cols:
        if column not in df_out.columns:
            df_out[column] = pd.NA

    tasks = [(idx, raw_val, session) for idx, raw_val in df_out[rs_col].items()]

    allele_results = {}
    enrich_by_idx: dict[int, dict] = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_item, task) for task in tasks]
        for fut in as_completed(futures):
            idx, ref, alts, enrich = fut.result()
            if ref is not None:
                allele_results[idx] = (ref, alts)
            if enrich:
                enrich_by_idx[idx] = enrich

    if allele_results:
        ref_series = pd.Series({i: r for i, (r, _) in allele_results.items()})
        alt_series = pd.Series({i: a for i, (_, a) in allele_results.items()})
        df_out.loc[ref_series.index, ref_col] = ref_series
        df_out.loc[alt_series.index, alt_col] = alt_series

    for c in new_cols:
        series = pd.Series({i: d.get(c) for i, d in enrich_by_idx.items() if d.get(c) is not None})
        if not series.empty:
            df_out.loc[series.index, c] = series

    df_out.to_csv(output_path, index=False)
    return df_out

if __name__ == "__main__":
    input_files_path = [
        './files/clinvar_BRCA1_GRCh38_filtered.tsv',
        './files/clinvar_BRCA2_GRCh38_filtered.tsv'
    ]

    df_brca1 = pd.read_csv(input_files_path[0], sep="\t", dtype=str)
    df_brca1_upd = update_df_inplace_with_vep(
        df=df_brca1,
        output_path='./files/clinvar_BRCA1_GRCh38_with_alleles_and_proteins.csv',
        max_workers=20,
    )

    df_brca2 = pd.read_csv(input_files_path[1], sep="\t", dtype=str)
    df_brca2_upd = update_df_inplace_with_vep(
        df=df_brca2,
        output_path='./files/clinvar_BRCA2_GRCh38_with_alleles_and_proteins1.csv',
        max_workers=20,
    )

    def _count_filled(df):
        ref_ok = df["ReferenceAllele"].notna().sum()
        alt_ok = df["AlternateAllele"].notna().sum()
        coding_ok = (df["IsCoding"] == True).sum()
        return ref_ok, alt_ok, coding_ok

    r1, a1, c1 = _count_filled(df_brca1_upd)
    r2, a2, c2 = _count_filled(df_brca2_upd)
    print(f"BRCA1 -> ref preenchidos: {r1}, alt preenchidos: {a1}, IsCoding=True: {c1}")
    print(f"BRCA2 -> ref preenchidos: {r2}, alt preenchidos: {a2}, IsCoding=True: {c2}")
