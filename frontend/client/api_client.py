import os
import requests
from enum import Enum
from typing import Dict, Any

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

    def download_dataset_bytes(
        self,
        genes: GeneSelection,
        force: bool = False,
        is_training_dt: bool = False,
        both : bool = True
    ) -> bytes:
        url = f"{self.base_url}/api/datasets/download"

        if is_training_dt:
            payload = {
                "genes": genes,
                "force": force, 
                "both": True,
                "is_training_dt" : True
            }
        else:
            payload = payload = {
                "genes": genes,
                "force": force, 
                "both": True,
                "is_training_dt" : False
            }

        r = requests.post(url, json=payload, timeout=3600)
        r.raise_for_status()
        return r.content

    def predict_from_fasta_file(
        self,
        gene: str,
        file_name: str,
        file_bytes: bytes,
        timeout: int = 300
    ) -> Dict[str, Any]:

        url = f"{self.base_url}/api/ml/predict"
        files = {
            "file": (file_name, file_bytes, "application/octet-stream")
        }
        data = {"gene": gene}
        r = requests.post(url, files=files, data=data, timeout=timeout)
        r.raise_for_status()
        return r.content