import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tenacity import retry, wait_exponential, stop_after_attempt

# Endpoint do VEP por rsID
ENSEMBL_VEP_ID_URL = "https://rest.ensembl.org/vep/human/id/{rsid}"
HEADERS = {"Accept": "application/json"}

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

@retry(stop=stop_after_attempt(20), wait=wait_exponential(min=1, max=60))
def fetch_alleles_for_rsid(rsid: str, session: requests.Session, timeout=30):

    url = ENSEMBL_VEP_ID_URL.format(rsid=rsid)
    try:
        request = session.get(url, headers=HEADERS, timeout=timeout)
        data = request.json()
    except Exception:
        raise Exception
        return None

    if not isinstance(data, list) or not data:
        return None

    rec = data[0]
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
        return None

    parts = allele_string.split("/")
    if len(parts) < 2:
        return None

    ref = parts[0]
    alts_str = ",".join(parts[1:])
    return ref, alts_str

def process_item(args):
    idx, raw_val, session = args

    rsid = normalize_rs(raw_val)
    if not rsid:
        return idx, None, None

    result = fetch_alleles_for_rsid(rsid, session)
    print("Aqui", result)
    if result:
        ref, alts = result
        return idx, ref, alts
    return idx, None, None

def update_df_inplace_with_vep(
    df: pd.DataFrame,
    output_path: str,
    rs_col: str = "RS# (dbSNP)",
    ref_col: str = "ReferenceAllele",
    alt_col: str = "AlternateAllele",
):

    df_out = df.copy()
    session = requests.Session()
    cache: dict[str, tuple[str, str] | None] = {}
    print("Aquiii")
    if ref_col not in df_out.columns:
        df_out[ref_col] = pd.NA
    if alt_col not in df_out.columns:
        df_out[alt_col] = pd.NA

    tasks = [
        (idx, raw_val, session) for idx, raw_val in df_out[rs_col].items()
    ]

    results = {}
    with ThreadPoolExecutor(max_workers=100) as executor:
        futures = [executor.submit(process_item, task) for task in tasks]
        i = 0
        for future in as_completed(futures):
            i += 1
            print(f"Processed: {i}")
            idx, ref, alts = future.result()
            if ref is not None:
                results[idx] = (ref, alts)
    print(results)
    ref_series = pd.Series({idx: ref for idx, (ref, _) in results.items()})
    alt_series = pd.Series({idx: alts for idx, (_, alts) in results.items()})
    df_out[ref_col].update(ref_series)
    df_out[alt_col].update(alt_series)
    
    df_out.to_csv(output_path)
    return df_out
    
if __name__ == "__main__":

    input_files_path = [
        './files/clinvar_BRCA1_GRCh38_filtered.tsv',
        './files/clinvar_BRCA2_GRCh38_filtered.tsv'
    ]

    df_brca1 = pd.read_csv(
        input_files_path[0], sep="\t", dtype=str,
    )
    
    df_brca1_upd = update_df_inplace_with_vep(
        df=df_brca1,
        output_path='./files/clinvar_BRCA1_GRCh38_with_alleles_3.csv'
    )

    df_brca2 = pd.read_csv(
        input_files_path[1], sep="\t", dtype=str,
    )
    
    df_brca2_upd = update_df_inplace_with_vep(
        df=df_brca2,
        output_path='./files/clinvar_BRCA2_GRCh38_with_alleles_3.csv'
    )

    def _count_filled(df):
        ref_ok = df["ReferenceAllele"].notna().sum()
        alt_ok = df["AlternateAllele"].notna().sum()
        return ref_ok, alt_ok

    #ref1, alt1 = _count_filled(df_brca1_upd)
    #print(f"BRCA1 -> ref preenchidos: {ref1}, alt preenchidos: {alt1}")
    ref2, alt2 = _count_filled(df_brca2_upd)
    print(f"BRCA2 -> ref preenchidos: {ref2}, alt preenchidos: {alt2}")