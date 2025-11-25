import io
import pandas as pd
import numpy as np

def _post_predict(client, gene, filename, content_bytes, content_type="application/octet-stream"):
    files = {
        "file": (filename, content_bytes, content_type)
    }
    data = {"gene": gene}
    return client.post("/api/ml/predict", files=files, data=data)

def test_predict_success_excel_response(client, monkeypatch):
    # mock ALLOWED_GENES
    monkeypatch.setattr("router.predict.ALLOWED_GENES", {"BRCA1", "BRCA2"})

    # mock do serviço para retornar DataFrame
    def fake_predict_service(gene, file_path):
        return pd.DataFrame({"col": ["a", "b"]})
    monkeypatch.setattr(
        "router.predict.run_predict_variants_service",
        fake_predict_service
    )

    fasta = b">seq\nACGTACGT\n"
    r = _post_predict(client, "BRCA1", "test.fasta", fasta, "text/plain")
    assert r.status_code == 200
    cd = r.headers.get("content-disposition", "")
    assert 'filename="prediction_BRCA1.xlsx"' in cd

    # valida XLSX básico
    buf = io.BytesIO(r.content)
    df = pd.read_excel(buf, sheet_name="data", dtype=str)
    assert list(df.columns) == ["col"]
    assert df.shape == (2, 1)

def test_predict_invalid_gene(client, monkeypatch):
    monkeypatch.setattr("router.predict.ALLOWED_GENES", {"BRCA1", "BRCA2"})
    fasta = b">s\nAC\n"
    r = _post_predict(client, "TP53", "x.fasta", fasta)
    assert r.status_code == 422

    detail = r.json()["detail"]
    assert isinstance(detail, list)
    err = detail[0]
    assert err["loc"][-1] == "gene"
    assert "BRCA1" in err["msg"] and "BRCA2" in err["msg"]


def test_predict_empty_file(client, monkeypatch):
    monkeypatch.setattr("router.predict.ALLOWED_GENES", {"BRCA1", "BRCA2"})
    r = _post_predict(client, "BRCA1", "empty.fasta", b"")
    assert r.status_code == 422
    assert r.json()["detail"] == "File should not be empty."

def test_predict_invalid_extension(client, monkeypatch):
    monkeypatch.setattr("router.predict.ALLOWED_GENES", {"BRCA1", "BRCA2"})
    r = _post_predict(client, "BRCA1", "bad.txt", b">x\nAA\n")
    assert r.status_code == 422
    assert "Invalid file extension" in r.json()["detail"]

def test_predict_service_returns_non_dataframe(client, monkeypatch):
    monkeypatch.setattr("router.predict.ALLOWED_GENES", {"BRCA1", "BRCA2"})
    def bad_service(gene, file_path):
        return {"not": "a dataframe"}
    monkeypatch.setattr(
        "router.predict.run_predict_variants_service",
        bad_service
    )
    r = _post_predict(client, "BRCA2", "ok.fasta", b">x\nAAA\n")
    assert r.status_code == 500
    assert r.json()["detail"] == "Internal error: prediction did not return a DataFrame."
