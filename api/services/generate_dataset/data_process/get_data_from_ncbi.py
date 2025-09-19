import os, gzip, requests
import pandas as pd

from api.services.settings import (
    OUTDIR,
    CLINVAR_TAB_URL
)

os.makedirs(OUTDIR, exist_ok=True)

def download_variant_summary(dst_path=f"{OUTDIR}/variant_summary.txt.gz", chunk=1<<20):
    with requests.get(CLINVAR_TAB_URL, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dst_path, "wb") as f:
            for b in r.iter_content(chunk_size=chunk):
                if b:
                    f.write(b)
    return dst_path

def load_variant_summary(gz_path=f"{OUTDIR}/variant_summary.txt.gz"):
    with gzip.open(gz_path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", dtype=str, low_memory=False)
    return df

def filter_gene(df, gene="BRCA1",
                assembly="GRCh38",
                clin_sigs=("Pathogenic","Likely pathogenic","Benign","Likely benign","Uncertain significance"),
                var_types=("single nucleotide variant","Indel","Deletion","Duplication","Insertion")):
    g = df.copy()
    g = g[g["Assembly"] == assembly]
    g = g[g["GeneSymbol"].fillna("").str.contains(rf"\b{gene}\b", case=False, regex=True)]
    if clin_sigs:
        cs_pat = "|".join([s.lower() for s in clin_sigs])
        g = g[g["ClinicalSignificance"].str.lower().str.contains(cs_pat, na=False)]
    if var_types:
        vt = set(v.lower() for v in var_types)
        g = g[g["Type"].str.lower().isin(vt)]
    cols = ["VariationID","AlleleID","GeneSymbol","ClinicalSignificance","ClinSigSimple",
            "Type","Origin","OriginSimple","Assembly","Chromosome","Start","Stop",
            "ReferenceAllele","AlternateAllele","ReviewStatus","RCVaccession","RS# (dbSNP)"]
    cols = [c for c in cols if c in g.columns]
    return g[cols].sort_values(["GeneSymbol","Chromosome","Start"])

def download_ncbi_gene_summary(gene):
    gz_path = download_variant_summary()
    df = load_variant_summary(gz_path)

    gene_df = filter_gene(df, gene=gene)

    gene_df.to_csv(f"{OUTDIR}/clinvar_{gene}_GRCh38.tsv", sep="\t", index=False)
    
    return gene_df

if __name__ == "__main__":
    download_ncbi_gene_summary()
