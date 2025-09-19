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

def run_generate_dataset_service(
    genes : list,
    *,
    force : bool = False,
    both: bool = False
):
    out_paths = {}
    for gene in genes:
        initial_df_path = f'{OUTDIR}/clinvar_{gene}_GRCh38.tsv'
        if force or not is_ready_file(initial_df_path):
            print("Entrou no if")
            download_ncbi_gene_summary(gene=gene)

        cleaned_df_path = f'{OUTDIR}/clinvar_{gene}_GRCh38_filtered.tsv'
        if force or not is_ready_file(cleaned_df_path):
            print("Entrou no if")
            cleaned_df_path = clean_data_per_submit_confidence(
                input_files=[initial_df_path],
                valid_confidence=VALID_CONFIDENCE_NCBI
            )
        cleaned_df = pd.read_csv(cleaned_df_path, sep="\t", dtype=str)

        with_alleles_df_path = f'{OUTDIR}/clinvar_{gene}_GRCh38_with_alleles_and_proteins.csv'
        if force or not is_ready_file(with_alleles_df_path):
            print("Entrou no if")
            with_alleles_df = update_df_inplace_with_vep(
                df=cleaned_df,
                output_path=with_alleles_df_path,
                max_workers=20,
            )
        else:
            with_alleles_df = pd.read_csv(with_alleles_df_path)

        enriched_uniprot_df_path = f'{OUTDIR}/clinvar_{gene}_GRCh38_with_uniprot.csv'
        if force or not is_ready_file(enriched_uniprot_df_path):
            print("Entrou no if")
            enrich_with_uniprot(
                with_alleles_df,
                output_path=enriched_uniprot_df_path,
                ensp_col="ENSP",
                enst_col="ENST",
                max_workers=20,
            )

        out_paths[gene] = {
            "initial": initial_df_path,
            "cleaned": cleaned_df_path,
            "with_alleles": with_alleles_df_path,
            "with_uniprot": enriched_uniprot_df_path,
        }

    merged_excel_path = None
    if both:
        csv_paths = [out_paths[gene]["with_uniprot"] for gene in genes]
        dfs = []
        for g, p in zip(genes, csv_paths):
            df = pd.read_csv(p)
            if "GeneSymbol" not in df.columns:
                df["GeneSymbol"] = g
            dfs.append(df)

        merged_df = pd.concat(dfs, ignore_index=True)
        merged_excel_path = Path(OUTDIR) / f"clinvar_{'_'.join(genes)}_GRCh38_merged.xlsx"

        with pd.ExcelWriter(merged_excel_path, engine="openpyxl") as writer:
            merged_df.to_excel(writer, index=False, sheet_name="merged_genes")

    return {"by_gene": out_paths, "merged_excel": merged_excel_path}

if __name__ == "__main__":
    run_generate_dataset_service(genes=['BRCA1'])