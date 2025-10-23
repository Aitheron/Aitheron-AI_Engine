import pandas as pd
from typing import Iterable
from core.settings import IMPACT_RANK, TERM_TO_FLAG

DROP_COLS = {
    "VariationID","ClinSigSimple","ClinicalSignificance","Origin","OriginSimple",
    "Assembly","ReviewStatus","RCVaccession","RS# (dbSNP)","TranscriptSelectionBasis",
    "ENST","ENSP","CDSStart","CDSEnd","AminoAcids","Strand","HGVSc","HGVSp",
    "FASTA_CDS","FASTA_Protein","ProteinName","ProteinFunction","ProteinFeatures","ProteinDesc"
}

BASE_COL_ORDER = [
    "GeneSymbol","Chromosome","Start","End",
    "ReferenceAllele","AlternateAllele","Type",
    "ConsequenceTerms","Impact","ImpactRank",
    "IsCoding","Protein_Start","Protein_End",
    "CDNAStart","CDNAEnd","Codons",
    "len_ref","len_alt","length_change","Label"
]

def _as_str(s):
    return "" if pd.isna(s) else str(s)

def _normalize_is_coding(val):
    if pd.isna(val):
        return 0
    if isinstance(val, (int, float)):
        return 1 if int(val) == 1 else 0
    v = str(val).strip().lower()
    truthy = {"verdadeiro"}
    return 1 if v in truthy else 0

def _parse_terms(s: str) -> set:
    s = _as_str(s)
    parts = [p.strip().lower() for p in s.replace(",", "|").split("|") if p.strip()]
    return set(parts)

def _pick_alt_allele(is_coding: int, alt: str, variant_allele: str) -> str:
    alt = _as_str(alt)
    va = _as_str(variant_allele)
    if is_coding == 1 and va and ("," not in va) and ("|" not in va):
        return va
    if not alt:
        return ""
    return alt.split(",")[0].strip()

def _len_allele(x: str) -> int:
    return len(_as_str(x).replace("-", ""))

def _map_label(clinsig: str) -> int:
    s = _as_str(clinsig).lower()
    if any(k in s for k in ["pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"]):
        return 1
    if "uncertain" in s:
        return 1
    if any(k in s for k in ["benign", "likely_benign", "benign/likely_benign"]):
        return 0
    return 1

def _impact_to_rank(impact: str):
    s = _as_str(impact).upper()
    return IMPACT_RANK.get(s, "")

def _apply_term_flags(terms: set) -> dict:
    d = { col: 0 for col in TERM_TO_FLAG.values() }
    for term, col in TERM_TO_FLAG.items():
        if term.lower() in terms:
            d[col] = 1
    return d

def build_training_csv(input_csv: str, output_csv: str) -> pd.DataFrame:
    df = pd.read_csv(input_csv, dtype=str, keep_default_na=False)
    cols = set(df.columns)

    if "Stop" in cols and "End" not in cols:
        df = df.rename(columns={"Stop": "End"})
    if "ProteinStart" in df.columns:
        df = df.rename(columns={"ProteinStart": "Protein_Start"})
    if "ProteinEnd" in df.columns:
        df = df.rename(columns={"ProteinEnd": "Protein_End"})

    for c in ["GeneSymbol","Chromosome","Start","End","ReferenceAllele","AlternateAllele",
              "Type","ConsequenceTerms","Impact","IsCoding","Protein_Start","Protein_End",
              "CDNAStart","CDNAEnd","Codons","VariantAllele","ClinicalSignificance"]:
        if c not in df.columns:
            df[c] = ""

    df["IsCoding"] = df["IsCoding"].apply(_normalize_is_coding)
    df["AlternateAllele"] = [
        _pick_alt_allele(ic, alt, va)
        for ic, alt, va in zip(df["IsCoding"], df["AlternateAllele"], df["VariantAllele"])
    ]

    df["len_ref"] = df["ReferenceAllele"].apply(_len_allele)
    df["len_alt"] = df["AlternateAllele"].apply(_len_allele)
    df["length_change"] = df["len_alt"] - df["len_ref"]

    df["Label"] = df["ClinicalSignificance"].apply(_map_label)

    df["ImpactRank"] = df["Impact"].apply(_impact_to_rank)

    term_sets = df["ConsequenceTerms"].apply(_parse_terms)
    flag_rows = term_sets.apply(_apply_term_flags).tolist()
    flags_df = pd.DataFrame(flag_rows)
    for c in flags_df.columns:
        df[c] = flags_df[c].astype(int)

    to_drop: Iterable[str] = [c for c in DROP_COLS if c in df.columns]
    if to_drop:
        df = df.drop(columns=to_drop)

    numeric_can_be_blank = ["Start","End","Protein_Start","Protein_End","CDNAStart","CDNAEnd"]
    for c in numeric_can_be_blank:
        if c in df.columns:
            df[c] = df[c].where(df[c].astype(str).str.len() > 0, "")

    if "VariantAllele" in df.columns:
        df = df.drop(columns=["VariantAllele"])
    if "ClinicalSignificance" in df.columns:
        df = df.drop(columns=["ClinicalSignificance"])

    final_cols = [c for c in BASE_COL_ORDER if c in df.columns]
    other_cols = [c for c in df.columns if c not in final_cols]
    df = df[final_cols + other_cols]

    df.to_csv(output_csv, index=False)
    return df

if __name__ == "__main__":
    inp = "./files/clinvar_BRCA1_BRCA2_GRCh38_merged.csv"
    out = "training_dataset.csv"
    df_out = build_training_csv(inp, out)
    print(f"OK: {len(df_out):,} linhas -> {out}")
