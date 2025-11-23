import io
import pandas as pd

def _post_download(client, body):
    return client.post("/api/datasets/download", json=body)

def test_download_excel_success_training(tmp_path, client, monkeypatch):
    # cria CSV temporário
    csv_path = tmp_path / "out.csv"
    pd.DataFrame([{"a": "1"}, {"a": "2"}]).to_csv(csv_path, index=False)

    # mock do serviço para retornar caminho do CSV
    def fake_run_generate_dataset_service(genes, force, is_training_dt):
        # genes devem chegar ordenados e sem duplicados no router
        assert genes == ["BRCA1", "BRCA2"]  # garantia do teste
        assert is_training_dt is True
        return {"dataset": str(csv_path)}
    monkeypatch.setattr(
        "router.datasets.run_generate_dataset_service",
        fake_run_generate_dataset_service
    )

    body = {"genes": ["BRCA2", "BRCA1"], "is_training_dt": True}
    r = _post_download(client, body)
    assert r.status_code == 200
    cd = r.headers.get("content-disposition", "")
    assert "attachment" in cd
    # nome deve refletir genes ordenados e sufixo "training"
    assert 'filename="dataset_BRCA1_BRCA2_training.xlsx"' in cd

    # conteúdo deve ser um XLSX válido com a planilha "data"
    content = io.BytesIO(r.content)
    df = pd.read_excel(content, sheet_name="data", dtype=str)
    assert df.shape == (2, 1)
    assert list(df.columns) == ["a"]

def test_download_excel_success_merged(tmp_path, client, monkeypatch):
    csv_path = tmp_path / "merged.csv"
    pd.DataFrame([{"x": "y"}]).to_csv(csv_path, index=False)

    def fake_run_generate_dataset_service(genes, force, is_training_dt):
        assert genes == ["BRCA1"]
        assert is_training_dt is False
        return {"dataset": str(csv_path)}
    monkeypatch.setattr(
        "router.datasets.run_generate_dataset_service",
        fake_run_generate_dataset_service
    )

    body = {"genes": ["BRCA1"], "is_training_dt": False}
    r = _post_download(client, body)
    assert r.status_code == 200
    cd = r.headers.get("content-disposition", "")
    assert 'filename="dataset_BRCA1_merged.xlsx"' in cd

def test_download_excel_missing_csv(client, monkeypatch, tmp_path):
    # serviço retorna caminho inexistente -> 500 CSV not found
    def fake_service(*args, **kwargs):
        return {"dataset": str(tmp_path / "does_not_exist.csv")}
    monkeypatch.setattr(
        "router.datasets.run_generate_dataset_service",
        fake_service
    )

    r = _post_download(client, {"genes": ["BRCA1"]})
    assert r.status_code == 500
    assert r.json()["detail"] == "CSV not found."

def test_download_excel_service_exception(client, monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("explode")
    monkeypatch.setattr(
        "router.datasets.run_generate_dataset_service",
        boom
    )
    r = _post_download(client, {"genes": ["BRCA1"]})
    assert r.status_code == 500
    # detalhe deve carregar a mensagem de erro
    assert "explode" in r.json()["detail"]
