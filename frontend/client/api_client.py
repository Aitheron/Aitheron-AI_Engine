import os
import requests
from enum import Enum

class GeneSelection(str, Enum):
    BRCA1 = "BRCA1"
    BRCA2 = "BRCA2"
    BOTH  = "BOTH"

def resolve_api_base_url(default: str = "http://localhost:8000") -> str:
    """Safe resolver for API base URL without crashing when secrets.toml is absent."""
    # 1) env var wins
    env_val = os.getenv("API_BASE_URL")
    if env_val:
        return env_val.rstrip("/")

    # 2) try Streamlit secrets (if present)
    try:
        # late import to avoid hard dependency at import time
        import streamlit as st  # type: ignore
        try:
            api_val = st.secrets["API_BASE_URL"]  # will raise if no secrets file
            return str(api_val).rstrip("/")
        except Exception:
            pass
    except Exception:
        pass

    # 3) default
    return default.rstrip("/")

class APIClient:
    def __init__(self, base_url: str | None = None):
        self.base_url = (base_url or resolve_api_base_url()).rstrip("/")

    def download_dataset_bytes(self, gene: GeneSelection, force: bool = False) -> bytes:
        url = f"{self.base_url}/api/datasets/download"
        params = {"gene": gene.value, "force": str(force).lower()}
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        return r.content
