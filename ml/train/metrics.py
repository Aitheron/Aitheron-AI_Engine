import numpy as np
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score, precision_score,
    recall_score, matthews_corrcoef, balanced_accuracy_score, accuracy_score
)

def compute_binary_metrics(y_true: np.ndarray, y_prob: np.ndarray, threshold: float) -> dict:
    y_pred = (y_prob >= threshold).astype(int)
    metrics = {}
    try:
        metrics["auroc"] = roc_auc_score(y_true, y_prob)
    except Exception:
        metrics["auroc"] = float("nan")
    try:
        metrics["auprc"] = average_precision_score(y_true, y_prob)
    except Exception:
        metrics["auprc"] = float("nan")

    metrics["recall"] = recall_score(y_true, y_pred, zero_division=0)  # sensibilidade
    metrics["precision"] = precision_score(y_true, y_pred, zero_division=0)
    metrics["f1"] = f1_score(y_true, y_pred, zero_division=0)
    metrics["balanced_acc"] = balanced_accuracy_score(y_true, y_pred)
    metrics["mcc"] = matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 else float("nan")
    metrics["accuracy"] = accuracy_score(y_true, y_pred)
    metrics["threshold"] = threshold
    return metrics

def pick_best_threshold_by_f1(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    candidates = np.linspace(0.1, 0.9, 17)
    best_thr, best_f1 = 0.5, -1.0
    for t in candidates:
        f1 = f1_score(y_true, (y_prob >= t).astype(int), zero_division=0)
        if f1 > best_f1:
            best_f1, best_thr = f1, t
    return float(best_thr)
