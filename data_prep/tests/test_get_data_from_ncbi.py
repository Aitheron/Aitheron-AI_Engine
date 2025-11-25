import gzip
from io import BytesIO
from pathlib import Path
import pandas as pd

import data_prep.get_data_from_ncbi as mod

def test_filter_gene_basic():
    df = pd.DataFrame({
        "Assembly": ["GRCh38","GRCh37","GRCh38"],
        "GeneSymbol": ["BRCA1","BRCA1","TP53"],
        "ClinicalSignificance": ["Pathogenic","Benign","Benign"],
        "Type": ["single nucleotide variant","Deletion","Insertion"],
        "Chromosome": ["17","17","17"],
        "Start": ["1","2","3"],
        "Stop": ["1","2","3"],
    })
    out = mod.filter_gene(df, gene="BRCA1", assembly="GRCh38")
    assert set(out["GeneSymbol"]) == {"BRCA1"}
    assert (out["Assembly"] == "GRCh38").all()

def test_download_ncbi_gene_summary_writes_tsv(tmp_path, monkeypatch):
    monkeypatch.setattr("core.settings.OUTDIR", str(tmp_path), raising=False)
    monkeypatch.setattr("data_prep.get_data_from_ncbi.OUTDIR", str(tmp_path), raising=False)

    content = "Assembly\tGeneSymbol\tClinicalSignificance\tType\tChromosome\tStart\tStop\n" \
              "GRCh38\tBRCA1\tPathogenic\tsingle nucleotide variant\t17\t100\t100\n"
    buf = BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content.encode("utf-8"))
    buf.seek(0)

    class FakeResp:
        def __init__(self, data): self.data=data
        def __enter__(self): return self
        def __exit__(self, *a): pass
        def raise_for_status(self): pass
        def iter_content(self, chunk_size): yield self.data.getvalue()

    monkeypatch.setattr(
        "data_prep.get_data_from_ncbi.requests.get",
        lambda *a, **k: FakeResp(buf)
    )

    out_df = mod.download_ncbi_gene_summary("BRCA1")
    tsv = Path(tmp_path) / "clinvar_BRCA1_GRCh38.tsv"
    assert Path(tsv).exists()
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    assert len(df) == 1
    assert df.iloc[0]["GeneSymbol"] == "BRCA1"
