import io
from pathlib import Path
import pandas as pd
import pytest

import services.generate_dataset_service as svc


@pytest.fixture(autouse=True)
def patch_outdir(tmp_path, monkeypatch):
    monkeypatch.setattr("core.settings.OUTDIR", str(tmp_path), raising=False)
    monkeypatch.setattr("services.generate_dataset_service.OUTDIR", str(tmp_path), raising=False)
    monkeypatch.setattr(
        "data_prep.utils.validations.is_ready_file",
        lambda p: Path(p).exists()
    )
    return tmp_path

def test_drop_invalid_rs_filters_rows():
    df = pd.DataFrame({"RS# (dbSNP)": ["1", "-1", "3"], "x": [1, 2, 3]})
    out = svc._drop_invalid_rs(df)
    assert list(out["RS# (dbSNP)"]) == ["1", "3"]


def test_returns_merged_if_exists_without_training(tmp_path, monkeypatch):
    outdir = Path(svc.OUTDIR if hasattr(svc, "OUTDIR") else tmp_path)
    genes = ["BRCA2", "BRCA1"]
    base = f"clinvar_{'_'.join(sorted(genes))}_GRCh38"
    merged = outdir / f"{base}_merged.csv.gz"
    merged.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"a": "1"}]).to_csv(merged, index=False)

    res = svc.run_generate_dataset_service(genes=genes, force=False, is_training_dt=False)
    assert res["kind"] == "merged"
    assert Path(res["dataset"]).name == merged.name


def test_returns_training_from_merged_and_filters_invalid_rs(tmp_path, monkeypatch):
    outdir = Path(svc.OUTDIR if hasattr(svc, "OUTDIR") else tmp_path)
    genes = ["BRCA1"]
    base = f"clinvar_{'_'.join(sorted(genes))}_GRCh38"
    merged = outdir / f"{base}_merged.csv.gz"
    merged.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"RS# (dbSNP)": ["-1", "2"]}]).to_csv(merged, index=False)

    # build_training_dataset escreve um CSV "sujo" que depois será filtrado pelo _ensure_training_from_merged
    def fake_build_training_dataset(in_path, out_path):
        pd.DataFrame({"RS# (dbSNP)": ["-1", "10"], "y": [0, 1]}).to_csv(out_path, index=False)
    monkeypatch.setattr("services.generate_dataset_service.build_training_dataset", fake_build_training_dataset)

    res = svc.run_generate_dataset_service(genes=genes, force=False, is_training_dt=True)
    assert res["kind"] == "training"
    training_path = Path(res["dataset"])
    assert training_path.exists()
    df_tr = pd.read_csv(training_path, dtype=str)
    # linha com -1 deve ter sido removida
    assert list(df_tr["RS# (dbSNP)"]) == ["10"]


def test_force_full_flow_creates_merged_and_sets_gene_column(tmp_path, monkeypatch):
    # mocks leves que apenas escrevem arquivos pequenos
    def fake_download_ncbi_gene_summary(gene):
        p = Path(tmp_path) / f"clinvar_{gene}_GRCh38.tsv"
        pd.DataFrame({"RS# (dbSNP)": ["1"], "GeneSymbol": [gene]}).to_csv(p, sep="\t", index=False)

    def fake_clean_data_per_submit_confidence(input_files, valid_confidence):
        src = Path(input_files[0])
        out = src.with_name(src.name.replace(".tsv", "_filtered.tsv"))
        df = pd.read_csv(src, sep="\t", dtype=str)
        df.to_csv(out, sep="\t", index=False)
        return str(out)

    def fake_update_df_inplace_with_vep(df, output_path, max_workers):
        out = Path(output_path)
        df2 = df.copy()
        df2["dummy"] = "ok"
        df2.to_csv(out, index=False)
        return df2

    def fake_enrich_with_uniprot(df, output_path, ensp_col, enst_col, max_workers):
        pd.DataFrame({"RS# (dbSNP)": ["1"], "GeneSymbol": df.get("GeneSymbol", pd.Series(["?"]))}).to_csv(output_path, index=False)

    monkeypatch.setattr("services.generate_dataset_service.download_ncbi_gene_summary", fake_download_ncbi_gene_summary)
    monkeypatch.setattr("services.generate_dataset_service.clean_data_per_submit_confidence", fake_clean_data_per_submit_confidence)
    monkeypatch.setattr("services.generate_dataset_service.update_df_inplace_with_vep", fake_update_df_inplace_with_vep)
    monkeypatch.setattr("services.generate_dataset_service.enrich_with_uniprot", fake_enrich_with_uniprot)

    res = svc.run_generate_dataset_service(genes=["BRCA1", "BRCA2"], force=True, is_training_dt=False)
    assert res["kind"] == "merged"
    merged_path = Path(res["dataset"])
    assert merged_path.exists()
    df = pd.read_csv(merged_path, dtype=str)
    # deve conter GeneSymbol preenchido para cada gene
    assert set(df["GeneSymbol"]) <= {"BRCA1", "BRCA2"}
