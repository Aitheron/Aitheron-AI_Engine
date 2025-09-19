import os
import requests
from enum import Enum

class GeneSelection(str, Enum):
    BRCA1 = "BRCA1"
    BRCA2 = "BRCA2"
    BOTH  = "BOTH"

def resolve_api_base_url(default: str = "http://localhost:8000") -> str:
    env_val = os.getenv("API_BASE_URL")
    if env_val:
        return env_val.rstrip("/")
    try:
        import streamlit as st
        try:
            api_val = st.secrets["API_BASE_URL"]
            return str(api_val).rstrip("/")
        except Exception:
            pass
    except Exception:
        pass
    return default.rstrip("/")

class APIClient:
    def __init__(self, base_url: str | None = None):
        self.base_url = (base_url or resolve_api_base_url()).rstrip("/")

    def download_dataset_bytes(self, gene: GeneSelection, force: bool = False) -> bytes:
        url = f"{self.base_url}/api/datasets/download"

        if gene is GeneSelection.BOTH:
            payload = {"genes": ["BRCA1", "BRCA2"], "force": force, "both": True}
        else:
            payload = {"genes": [gene.value], "force": force, "both": False}

        print("payload ->", payload)

        r = requests.post(url, json=payload, timeout=120)
        r.raise_for_status()
        return r.content
