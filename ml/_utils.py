import os
import json
import joblib
import numpy as np

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def save_json(obj, path: str):
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)

def save_joblib(obj, path: str):
    ensure_dir(os.path.dirname(path))
    joblib.dump(obj, path)

def unpack_class_probs_from_cumulative(p_cum: np.ndarray) -> np.ndarray:
    Km1 = p_cum.shape[1]
    K = Km1 + 1
    out = np.zeros((p_cum.shape[0], K), dtype=np.float32)
    out[:, 0] = 1.0 - p_cum[:, 0]
    for k in range(1, Km1):
        out[:, k] = p_cum[:, k - 1] - p_cum[:, k]
    out[:, Km1] = p_cum[:, Km1 - 1]
    out = np.clip(out, 0.0, 1.0)
    row_sum = out.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0.0] = 1.0
    out = out / row_sum
    return out

LABELS = {0:"Benigno", 1:"Possivelmente Benigno", 2:"VUS", 3:"Patogênico"}

def entropy(p):
    p = np.clip(p, 1e-12, 1.0) # Clip evita log de 0
    return -(p * np.log(p)).sum()

def pretty_line(head, idx, pvec, k):
    top1 = float(pvec.max())
    top2 = float(np.partition(pvec, -2)[-2]) if pvec.size >= 2 else 0.0
    margin = top1 - top2
    H = entropy(pvec)
    Hn = H / np.log(len(pvec))  # entropia normalizada
    return (
        f"{head} idx {idx:>5} | "
        f"P0={pvec[0]:.4f} P1={pvec[1]:.4f} P2={pvec[2]:.4f} P3={pvec[3]:.4f} | "
        f"argmax={k} ({LABELS[int(k)]}) | "
        f"top1={top1:.4f} top2={top2:.4f} margin={margin:.4f} conf={1.0-Hn:.4f}"
    )
