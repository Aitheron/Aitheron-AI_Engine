import pandas as pd
from pathlib import Path

from core.settings import OUTDIR, VALID_CONFIDENCE_NCBI
from data_prep.utils.validations import is_ready_file
from data_prep import (
    download_ncbi_gene_summary,
    clean_data_per_submit_confidence,
    update_df_inplace_with_vep,
    enrich_with_uniprot,
    build_training_dataset
)

def _drop_invalid_rs(df: pd.DataFrame, rs_col: str = "RS# (dbSNP)") -> pd.DataFrame:
    if rs_col in df.columns:
        mask_ok = df[rs_col].astype(str).str.strip().ne("-1")
        return df.loc[mask_ok].copy()
    return df

def run_generate_dataset_service(
    genes: list,
    *,
    force: bool = False,
    is_training_dt: bool = False
):
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("No genes provided.")
    genes_sorted = sorted(genes)
    base_name = f"clinvar_{'_'.join(genes_sorted)}_GRCh38"
    merged_csv_path = str(Path(OUTDIR) / f"{base_name}_merged.csv")
    training_csv_path = str(Path(OUTDIR) / f"{base_name}_training.csv")

    def _ensure_training_from_merged():
        build_training_dataset(merged_csv_path, training_csv_path)
        if is_ready_file(training_csv_path):
            df_tr = pd.read_csv(training_csv_path, dtype=str)
            df_tr = _drop_invalid_rs(df_tr)
            df_tr.to_csv(training_csv_path, index=False)

    if not force and is_ready_file(merged_csv_path):
        if is_training_dt:
            if not is_ready_file(training_csv_path):
                _ensure_training_from_merged()
            return {"dataset": training_csv_path, "kind": "training"}
        return {"dataset": merged_csv_path, "kind": "merged"}

    out_paths_by_gene = {}
    for gene in genes_sorted:
        initial_df_path = f"{OUTDIR}/clinvar_{gene}_GRCh38.tsv"
        if force or not is_ready_file(initial_df_path):
            download_ncbi_gene_summary(gene=gene)

        cleaned_df_path = f"{OUTDIR}/clinvar_{gene}_GRCh38_filtered.tsv"
        if force or not is_ready_file(cleaned_df_path):
            cleaned_df_path = clean_data_per_submit_confidence(
                input_files=[initial_df_path],
                valid_confidence=VALID_CONFIDENCE_NCBI
            )
        cleaned_df = pd.read_csv(cleaned_df_path, sep="\t", dtype=str)
        with_alleles_df_path = f"{OUTDIR}/clinvar_{gene}_GRCh38_with_alleles_and_proteins.csv"
        if force or not is_ready_file(with_alleles_df_path):
            with_alleles_df = update_df_inplace_with_vep(
                df=cleaned_df,
                output_path=with_alleles_df_path,
                max_workers=20,
            )
        else:
            with_alleles_df = pd.read_csv(with_alleles_df_path, dtype=str)
        enriched_uniprot_df_path = f"{OUTDIR}/clinvar_{gene}_GRCh38_with_uniprot.csv"
        if force or not is_ready_file(enriched_uniprot_df_path):
            enrich_with_uniprot(
                with_alleles_df,
                output_path=enriched_uniprot_df_path,
                ensp_col="ENSP",
                enst_col="ENST",
                max_workers=20,
            )

        if is_ready_file(enriched_uniprot_df_path):
            df_tmp = pd.read_csv(enriched_uniprot_df_path, dtype=str)
            df_tmp = _drop_invalid_rs(df_tmp)
            df_tmp.to_csv(enriched_uniprot_df_path, index=False)

        out_paths_by_gene[gene] = {
            "with_uniprot": enriched_uniprot_df_path,
        }

    csv_paths = [out_paths_by_gene[g]["with_uniprot"] for g in genes_sorted]
    dfs = []
    for g, p in zip(genes_sorted, csv_paths):
        df = pd.read_csv(p, dtype=str)
        if "GeneSymbol" not in df.columns:
            df["GeneSymbol"] = g
        dfs.append(df)

    merged_df = pd.concat(dfs, ignore_index=True)
    merged_df = _drop_invalid_rs(merged_df)

    merged_df.to_csv(merged_csv_path, index=False)

    if is_training_dt:
        _ensure_training_from_merged()
        return {"dataset": training_csv_path, "kind": "training"}

    return {"dataset": merged_csv_path, "kind": "merged"}

