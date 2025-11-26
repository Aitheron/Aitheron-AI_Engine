import types, sys

from client.api_client import APIClient, resolve_api_base_url

def test_resolve_api_base_url_streamlit_secrets(monkeypatch):
    fake_st = types.ModuleType("streamlit")
    fake_st.secrets = {"API_BASE_URL": "https://svc.local///"}
    monkeypatch.setitem(sys.modules, "streamlit", fake_st)

    assert resolve_api_base_url() == "https://svc.local"


def test_resolve_api_base_url_streamlit_secrets(monkeypatch):
    # simula st.secrets
    fake_st = types.SimpleNamespace(secrets={"API_BASE_URL": "https://svc.local///"})
    monkeypatch.setitem(globals(), "streamlit", fake_st)
    monkeypatch.setenv("API_BASE_URL", "")

    # injeta módulo streamlit fake na importação interna
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "streamlit":
            return fake_st
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    try:
        assert resolve_api_base_url() == "https://svc.local"
    finally:
        monkeypatch.setenv("API_BASE_URL", "")


def test_download_dataset_bytes_posts_correct_payload(monkeypatch):
    sent = {}

    class FakeResp:
        def __init__(self, content=b"OK"): self.content = content
        def raise_for_status(self): 
            """Fake method used in tests: does nothing to simulate a successful HTTP response."""
            pass

    def fake_post(url, json=None, timeout=None, **_):
        sent["url"] = url
        sent["json"] = json
        sent["timeout"] = timeout
        return FakeResp(b"EXCEL")

    monkeypatch.setattr("client.api_client.requests.post", fake_post)
    cli = APIClient(base_url="http://localhost:9999")

    out = cli.download_dataset_bytes(
        genes=["BRCA1", "BRCA2"],
        force=False,
        is_training_dt=True,
        both=True
    )

    assert out == b"EXCEL"
    assert sent["url"].endswith("/api/datasets/download")
    # checa chaves importantes do payload
    assert sent["json"]["genes"] == ["BRCA1", "BRCA2"]
    assert sent["json"]["is_training_dt"] is True
    assert sent["timeout"] == 3600


def test_predict_from_fasta_file_posts_multipart(monkeypatch):
    captured = {}

    class FakeResp:
        def __init__(self, content=b"XLSX"): self.content = content
        def raise_for_status(self): 
            """Fake method used in tests: does nothing to simulate a successful HTTP response."""
            pass

    def fake_post(url, files=None, data=None, timeout=None, **_):
        captured["url"] = url
        captured["files"] = files
        captured["data"] = data
        captured["timeout"] = timeout
        return FakeResp()

    monkeypatch.setattr("client.api_client.requests.post", fake_post)
    cli = APIClient(base_url="http://api:8000")

    out = cli.predict_from_fasta_file("BRCA1", "p.fa", b">x\nACGT\n")
    assert out == b"XLSX"
    assert captured["url"].endswith("/api/ml/predict")
    assert captured["data"] == {"gene": "BRCA1"}
    assert "file" in captured["files"]
    name, payload, mime = captured["files"]["file"]
    assert name == "p.fa" and isinstance(payload, (bytes, bytearray)) and "octet-stream" in mime
    assert captured["timeout"] == 300
