import pandas as pd
from pathlib import Path

from services.settings import OUTDIR, VALID_CONFIDENCE_NCBI
from services.generate_dataset.data_process.utils.validations import is_ready_file
from services.generate_dataset.data_process import (
    download_ncbi_gene_summary,
    clean_data_per_submit_confidence,
    update_df_inplace_with_vep,
    enrich_with_uniprot
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
    both: bool = False
):
    out_paths = {}

    merged_base = Path(OUTDIR) / f"clinvar_{'_'.join(genes)}_GRCh38_merged"
    merged_csv_path = str(merged_base.with_suffix(".csv"))
    merged_excel_path = str(merged_base.with_suffix(".xlsx"))

    # curto-circuito se já existe cache dos dois formatos
    if both and not force and is_ready_file(merged_csv_path) and is_ready_file(merged_excel_path):
        out_paths = {
            gene: {
                "initial": f"{OUTDIR}/clinvar_{gene}_GRCh38.tsv",
                "cleaned": f"{OUTDIR}/clinvar_{gene}_GRCh38_filtered.tsv",
                "with_alleles": f"{OUTDIR}/clinvar_{gene}_GRCh38_with_alleles_and_proteins.csv",
                "with_uniprot": f"{OUTDIR}/clinvar_{gene}_GRCh38_with_uniprot.csv",
            }
            for gene in genes
        }

        # limpeza de registros com RS# (dbSNP) == "-1"
        for gene in genes:
            p = out_paths[gene]["with_uniprot"]
            if is_ready_file(p):
                df_tmp = pd.read_csv(p, dtype=str)
                df_tmp = _drop_invalid_rs(df_tmp)
                df_tmp.to_csv(p, index=False)
        if is_ready_file(merged_csv_path):
            df_merge = pd.read_csv(merged_csv_path, dtype=str)
            df_merge = _drop_invalid_rs(df_merge)
            df_merge.to_csv(merged_csv_path, index=False)
            df_merge.to_excel(merged_excel_path, index=False)

        return {"by_gene": out_paths, "merged_csv": merged_csv_path, "merged_excel": merged_excel_path}

    for gene in genes:
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

        # limpeza de RS# == "-1" para os arquivos individuais
        if is_ready_file(enriched_uniprot_df_path):
            df_tmp = pd.read_csv(enriched_uniprot_df_path, dtype=str)
            df_tmp = _drop_invalid_rs(df_tmp)
            df_tmp.to_csv(enriched_uniprot_df_path, index=False)

        out_paths[gene] = {
            "initial": initial_df_path,
            "cleaned": cleaned_df_path,
            "with_alleles": with_alleles_df_path,
            "with_uniprot": enriched_uniprot_df_path,
        }

    merged_csv_path = None
    merged_excel_path = None

    # gera merge em CSV (canônico p/ treino) e XLSX (download/visualização)
    if both:
        csv_paths = [out_paths[gene]["with_uniprot"] for gene in genes]
        dfs = []
        for g, p in zip(genes, csv_paths):
            df = pd.read_csv(p, dtype=str)
            if "GeneSymbol" not in df.columns:
                df["GeneSymbol"] = g
            dfs.append(df)

        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df = _drop_invalid_rs(merged_df)

        merged_base = Path(OUTDIR) / f"clinvar_{'_'.join(genes)}_GRCh38_merged"
        merged_csv_path = str(merged_base.with_suffix(".csv"))
        merged_excel_path = str(merged_base.with_suffix(".xlsx"))

        merged_df.to_csv(merged_csv_path, index=False)
        merged_df.to_excel(merged_excel_path, index=False)

    return {"by_gene": out_paths, "merged_csv": merged_csv_path, "merged_excel": merged_excel_path}


if __name__ == "__main__":
    run_generate_dataset_service(genes=["BRCA1", "BRCA2"])
