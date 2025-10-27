import numpy as np
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score, precision_score,
    recall_score, matthews_corrcoef, balanced_accuracy_score, accuracy_score, confusion_matrix, precision_recall_curve
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

def pick_threshold_for_recall(y_true, y_prob, target_recall=0.999):
    from sklearn.metrics import precision_recall_curve
    _, recall, thresholds = precision_recall_curve(y_true, y_prob)
    recall_thr = recall[1:]
    if (recall_thr >= target_recall).any():
        idx = np.where(recall_thr >= target_recall)[0][-1]
        return float(thresholds[idx])
    return 1.0

def confusion_matrix_dict(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp), "matrix": [[int(tn), int(fp)], [int(fn), int(tp)]]}

def head_reports(y_val, prob0, prob1, mask_brca1, mask_brca2, thr0, thr1):
    rep = {}
    if mask_brca1.any():
        y_b1 = y_val[mask_brca1]
        p_b1 = prob0[mask_brca1]
        pred_b1 = (p_b1 >= thr0).astype(int)
        rep["BRCA1"] = compute_binary_metrics(y_b1, p_b1, thr0)
        rep["BRCA1_confusion"] = confusion_matrix_dict(y_b1, pred_b1)
    if mask_brca2.any():
        y_b2 = y_val[mask_brca2]
        p_b2 = prob1[mask_brca2]
        pred_b2 = (p_b2 >= thr1).astype(int)
        rep["BRCA2"] = compute_binary_metrics(y_b2, p_b2, thr1)
        rep["BRCA2_confusion"] = confusion_matrix_dict(y_b2, pred_b2)
    rep["thresholds"] = {"BRCA1": float(thr0), "BRCA2": float(thr1)}
    return rep

def pick_threshold_under_recall(y_true, y_prob, recall_target=0.999, objective="precision", tie_break="highest_threshold"):
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    p, r, thr = precision_recall_curve(y_true, y_prob)  # r tem len = len(thr)+1
    r_thr = r[1:]; p_thr = p[1:]  # alinhar com thr

    idxs = np.where(r_thr >= recall_target)[0]
    if len(idxs) == 0:
        return 1.0, {"precision": float("nan"), "recall": float(np.max(r_thr) if len(r_thr) else 0.0),
                     "f1": float("nan"), "accuracy": float("nan")}

    # candidatos
    thr_cand = thr[idxs]
    p_cand   = p_thr[idxs]
    r_cand   = r_thr[idxs]

    if objective == "precision":
        scores = p_cand
    elif objective == "f1":
        scores = (2 * p_cand * r_cand) / np.clip(p_cand + r_cand, 1e-12, None)
    elif objective == "accuracy":
        # calcula acc por limiar
        scores = np.array([accuracy_score(y_true, (y_prob >= t).astype(int)) for t in thr_cand])
    else:
        scores = p_cand  # padrão: precision

    best = np.argmax(scores)
    if tie_break == "highest_threshold":
        best_score = scores[best]
        ties = np.where(scores == best_score)[0]
        if len(ties) > 1:
            best = ties[np.argmax(thr_cand[ties])]

    t_star = float(thr_cand[best])
    y_pred = (y_prob >= t_star).astype(int)
    return t_star, {
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "accuracy": float(accuracy_score(y_true, y_pred))
    }