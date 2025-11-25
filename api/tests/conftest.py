# api/tests/conftest.py
import sys
from pathlib import Path
import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[2]
API_DIR = ROOT / "api"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from main import create_app

@pytest.fixture(scope="session")
def app():
    return create_app()

@pytest.fixture(scope="session")
def client(app):
    return TestClient(app)
