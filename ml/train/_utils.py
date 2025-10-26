import os
import json
import joblib
import torch

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def save_json(obj, path: str):
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)

def save_joblib(obj, path: str):
    ensure_dir(os.path.dirname(path))
    joblib.dump(obj, path)

def save_torch_model(model, path: str):
    ensure_dir(os.path.dirname(path))
    torch.save(model.state_dict(), path)
