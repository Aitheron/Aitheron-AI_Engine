import os
import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
    precision_score,
    recall_score,
    confusion_matrix,
)

def _plot_roc_pr_for_gene(gene, y_true, y_score, threshold, out_dir):
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    prec, rec, _ = precision_recall_curve(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)

    os.makedirs(out_dir, exist_ok=True)

    y_pred = (y_score >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    fpr_op = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    tpr_op = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    rec_op = recall_score(y_true, y_pred)
    prec_op = precision_score(y_true, y_pred) if (y_pred.sum() > 0) else 0.0

    plt.figure()
    plt.plot(fpr, tpr, label=f"AUROC = {roc_auc:.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.scatter(
        [fpr_op],
        [tpr_op],
        color="red",
        label=f"Threshold = {threshold:.3f}\nTPR = {tpr_op:.3f}, FPR = {fpr_op:.3f}",
    )
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC – {gene}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"roc_{gene}.png"), dpi=200)
    plt.close()

    plt.figure()
    plt.plot(rec, prec, label=f"AUPRC = {auprc:.4f}")
    plt.scatter(
        [rec_op],
        [prec_op],
        color="red",
        label=f"Threshold = {threshold:.3f}\nRecall = {rec_op:.3f}, Precision = {prec_op:.3f}",
    )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"Precision–Recall – {gene}")
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"pr_{gene}.png"), dpi=200)
    plt.close()

    return roc_auc, auprc, {
        "fpr_op": float(fpr_op),
        "tpr_op": float(tpr_op),
        "recall_op": float(rec_op),
        "precision_op": float(prec_op),
    }

def generate_evaluation_graphs(
    final_dir: str,
    y_full,
    g_full,
    p0_full,
    p1_full,
    thresholds: dict[str, float],
):

    os.makedirs(os.path.join(final_dir, "graphs"), exist_ok=True)

    y_patho = (y_full == 3).astype(np.int64)
    mask_b1 = (g_full == "BRCA1")
    mask_b2 = (g_full == "BRCA2")

    p0_patho = p0_full[:, 3]
    p1_patho = p1_full[:, 3]

    data = {
        "BRCA1": {
            "y_true": y_patho[mask_b1],
            "y_score": p0_patho[mask_b1],
        },
        "BRCA2": {
            "y_true": y_patho[mask_b2],
            "y_score": p1_patho[mask_b2],
        },
    }

    curves_dir = os.path.join(final_dir, "graphs")
    results = {}

    for gene, d in data.items():
        y_true = d["y_true"]
        y_score = d["y_score"]
        thr = thresholds[gene]
        roc_auc, auprc, op = _plot_roc_pr_for_gene(
            gene, y_true, y_score, thr, curves_dir
        )
        results[gene] = {
            "auroc_from_curve": float(roc_auc),
            "auprc_from_curve": float(auprc),
            "operating_point": op,
            "threshold": float(thr),
        }

    print("Curvas salvas em:", curves_dir)
    print("Resultados:", results)
    return results
