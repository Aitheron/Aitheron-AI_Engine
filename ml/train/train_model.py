import os
import math
import json
import numpy as np
import pandas as pd
import joblib
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from sklearn.model_selection import StratifiedKFold

from ml.train._settings import Config
from ml.train.data_loader import load_training_csv, prepare_feature_sets, Preprocessor
from ml.train.model import MultitaskMLP
from ml.train.evaluation_graph import generate_evaluation_graphs
from ml.train.metrics import (
    compute_binary_metrics,
    pick_best_threshold_by_f1,
    pick_threshold_under_recall,
    head_reports,
)
from ml._utils import (
    save_json,
    unpack_class_probs_from_cumulative,
)

from data_prep.prepare_training_dataset import build_training_dataset

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def make_mask_for_heads(gene_series: pd.Series) -> dict[str, np.ndarray]:
    genes = gene_series.to_numpy()
    mask0 = (genes == "BRCA1").astype(np.float32)
    mask1 = (genes == "BRCA2").astype(np.float32)
    return {"BRCA1": mask0, "BRCA2": mask1}

def to_cumulative_targets(y_ord: np.ndarray, n_classes: int = 4) -> np.ndarray:
    K = n_classes
    t = np.zeros((len(y_ord), K - 1), dtype=np.float32)
    for k in range(K - 1):
        t[:, k] = (y_ord > k).astype(np.float32)
    return t

def pos_weights_for_head_ordinal(y_ord: np.ndarray, g: np.ndarray, gene_name: str, n_classes: int = 4):
    mask = (g == gene_name)
    y_head = y_ord[mask]
    if y_head.size == 0:
        return torch.ones(n_classes - 1, dtype=torch.float32, device=device)
    K = n_classes
    t = to_cumulative_targets(y_head, n_classes=K)
    pos = t.sum(axis=0)
    neg = t.shape[0] - pos
    w = np.ones(K - 1, dtype=np.float32)
    for i in range(K - 1):
        if pos[i] <= 0:
            w[i] = 1.0
        else:
            w[i] = float(neg[i] / max(pos[i], 1.0))
    return torch.tensor(w, dtype=torch.float32, device=device)

def train_one_fold(cfg: Config, fold_idx: int, X_train, y_train, g_train, X_val, y_val, g_val, input_dim: int, pre: Preprocessor, outdir_fold: str) -> dict:
    model = MultitaskMLP(
        in_dim=input_dim,
        hidden_dims=cfg.train.hidden_dims,
        dropout=cfg.train.dropout,
        n_heads=2,
        n_classes=4
    ).to(device)

    pw0 = pos_weights_for_head_ordinal(y_train, g_train, "BRCA1", n_classes=4)
    pw1 = pos_weights_for_head_ordinal(y_train, g_train, "BRCA2", n_classes=4)

    criterion0 = nn.BCEWithLogitsLoss(pos_weight=pw0, reduction="none")
    criterion1 = nn.BCEWithLogitsLoss(pos_weight=pw1, reduction="none")
    optimizer = torch.optim.AdamW(model.parameters(), lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)

    y_tr_cum = to_cumulative_targets(y_train, n_classes=4)

    ds_train = TensorDataset(
        torch.tensor(X_train, dtype=torch.float32),
        torch.tensor(y_tr_cum, dtype=torch.float32),
        torch.tensor((g_train == "BRCA2").astype(np.int64))
    )
    dl_train = DataLoader(ds_train, batch_size=cfg.train.batch_size, shuffle=True, drop_last=False)

    ds_val = TensorDataset(
        torch.tensor(X_val, dtype=torch.float32),
        torch.tensor(y_val, dtype=torch.long),
    )
    dl_val = DataLoader(ds_val, batch_size=1024, shuffle=False, drop_last=False)

    best_score = -math.inf
    best_state = None
    patience = cfg.train.early_stopping_patience
    no_improve = 0

    mask_va = make_mask_for_heads(pd.Series(g_val))
    m0_va_np = mask_va["BRCA1"]; m1_va_np = mask_va["BRCA2"]

    eps = 1e-9

    for epoch in range(cfg.train.epochs):
        model.train()
        epoch_loss = 0.0
        seen_brca1 = 0
        seen_brca2 = 0
        for xb, yb_cum, gb in dl_train:
            xb = xb.to(device); yb_cum = yb_cum.to(device); gb = gb.to(device)

            logits0, logits1 = model(xb)
            l0 = criterion0(logits0, yb_cum)
            l1 = criterion1(logits1, yb_cum)

            m0 = (gb == 0).unsqueeze(1).float()
            m1 = (gb == 1).unsqueeze(1).float()

            b0 = int(m0.sum().item()); b1 = int(m1.sum().item())
            seen_brca1 += b0; seen_brca2 += b1

            loss0 = (l0 * m0).sum() / (m0.sum() * l0.shape[1] + eps)
            loss1 = (l1 * m1).sum() / (m1.sum() * l1.shape[1] + eps)

            loss = 0.0
            if torch.isfinite(loss0): loss = loss + loss0
            if torch.isfinite(loss1): loss = loss + loss1

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if isinstance(loss, torch.Tensor):
                epoch_loss += loss.item()

        print(f"[Fold {fold_idx}] ep{epoch+1} seen -> BRCA1={seen_brca1} BRCA2={seen_brca2}")

        model.eval()
        with torch.no_grad():
            xb_val = torch.tensor(X_val, dtype=torch.float32, device=device)
            lo0, lo1 = model(xb_val)
            p0_cum = torch.sigmoid(lo0).cpu().numpy()
            p1_cum = torch.sigmoid(lo1).cpu().numpy()
            p0_full = unpack_class_probs_from_cumulative(p0_cum)
            p1_full = unpack_class_probs_from_cumulative(p1_cum)
            prob0_patho = p0_full[:, 3]
            prob1_patho = p1_full[:, 3]

        y_val_patho = (y_val == 3).astype(np.int64)

        mask_b1 = np.asarray(m0_va_np, dtype=bool)
        if mask_b1.any():
            thr0, _ = pick_threshold_under_recall(
                y_val_patho[mask_b1], prob0_patho[mask_b1],
                recall_target=0.99,
                objective="precision"
            )
        thr1 = pick_best_threshold_by_f1(y_val_patho[m1_va_np > 0.5], prob1_patho[m1_va_np > 0.5]) if (m1_va_np > 0.5).any() else 0.5

        scores = []
        if (m0_va_np > 0.5).any():
            m0_metrics = compute_binary_metrics(y_val_patho[m0_va_np > 0.5], prob0_patho[m0_va_np > 0.5], thr0)
            scores.append(m0_metrics["auroc"])
        if (m1_va_np > 0.5).any():
            m1_metrics = compute_binary_metrics(y_val_patho[m1_va_np > 0.5], prob1_patho[m1_va_np > 0.5], thr1)
            scores.append(m1_metrics["auroc"])
        mean_auroc = np.nanmean(scores) if scores else -math.inf

        print(f"[Fold {fold_idx}] Época {epoch+1}/{cfg.train.epochs} - Loss treino: {epoch_loss:.4f} - Val AUROC médio: {mean_auroc:.4f} (best: {best_score:.4f})")

        if mean_auroc > best_score:
            print(f"[Fold {fold_idx}]   >> Novo melhor modelo na época {epoch+1}")
            best_score = mean_auroc
            best_state = {
                "model": {k: v.cpu() for k, v in model.state_dict().items()},
                "thr0": float(thr0),
                "thr1": float(thr1),
                "epoch": epoch
            }
            no_improve = 0
        else:
            no_improve += 1
            if no_improve >= patience:
                print(f"[Fold {fold_idx}]   >> Early stopping após {patience} épocas sem melhora")
                break

    os.makedirs(outdir_fold, exist_ok=True)
    final_model_path = os.path.join(outdir_fold, f"mlp_fold{fold_idx}.pth")
    torch.save(best_state["model"], final_model_path)

    model.load_state_dict(best_state["model"])
    model.eval()
    with torch.no_grad():
        xb_val = torch.tensor(X_val, dtype=torch.float32, device=device)
        lo0, lo1 = model(xb_val)
        p0_cum = torch.sigmoid(lo0).cpu().numpy()
        p1_cum = torch.sigmoid(lo1).cpu().numpy()
        p0_full = unpack_class_probs_from_cumulative(p0_cum)
        p1_full = unpack_class_probs_from_cumulative(p1_cum)
        prob0_patho = p0_full[:, 3]
        prob1_patho = p1_full[:, 3]

    thr0 = best_state["thr0"]; thr1 = best_state["thr1"]
    report = {"fold": fold_idx, "best_epoch": best_state["epoch"], "mean_val_auroc": best_score}

    if (m0_va_np > 0.5).any():
        report["BRCA1"] = compute_binary_metrics((y_val == 3).astype(np.int64)[m0_va_np > 0.5], prob0_patho[m0_va_np > 0.5], thr0)
    if (m1_va_np > 0.5).any():
        report["BRCA2"] = compute_binary_metrics((y_val == 3).astype(np.int64)[m1_va_np > 0.5], prob1_patho[m1_va_np > 0.5], thr1)

    report["thresholds"] = {"BRCA1": thr0, "BRCA2": thr1}

    mask_b1_final = np.asarray(m0_va_np, dtype=bool)
    mask_b2_final = np.asarray(m1_va_np, dtype=bool)
    rep_heads = head_reports((y_val == 3).astype(np.int64), prob0_patho, prob1_patho, mask_b1_final, mask_b2_final, thr0, thr1)
    report.update(rep_heads)

    save_json(report, os.path.join(outdir_fold, f"metrics_fold{fold_idx}.json"))
    return report

def pos_weight_for_head_full(y, g, gene_name, device):
    return torch.ones(3, dtype=torch.float32, device=device)

def train_full_data(
    cfg: Config,
    X_full,
    y_full,
    g_full,
    init_state_dict: dict,
    pre_full: Preprocessor,
    final_epochs: int,
    thresholds_final: dict,
    final_outdir: str
) -> dict:
    model = MultitaskMLP(
        in_dim=X_full.shape[1],
        hidden_dims=cfg.train.hidden_dims,
        dropout=cfg.train.dropout,
        n_heads=2,
        n_classes=4
    ).to(device)
    model.load_state_dict(init_state_dict)

    pw0 = pos_weights_for_head_ordinal(y_full, g_full, "BRCA1", n_classes=4)
    pw1 = pos_weights_for_head_ordinal(y_full, g_full, "BRCA2", n_classes=4)
    criterion0 = nn.BCEWithLogitsLoss(pos_weight=pw0, reduction="none")
    criterion1 = nn.BCEWithLogitsLoss(pos_weight=pw1, reduction="none")
    optimizer = torch.optim.AdamW(model.parameters(), lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)

    y_full_cum = to_cumulative_targets(y_full, n_classes=4)

    ds = TensorDataset(
        torch.tensor(X_full, dtype=torch.float32),
        torch.tensor(y_full_cum, dtype=torch.float32),
        torch.tensor((g_full == "BRCA2").astype(np.int64))
    )
    dl = DataLoader(ds, batch_size=cfg.train.batch_size, shuffle=True, drop_last=False)

    eps = 1e-9
    for epoch in range(final_epochs):
        model.train()
        epoch_loss = 0.0
        for xb, yb_cum, gb in dl:
            xb = xb.to(device); yb_cum = yb_cum.to(device); gb = gb.to(device)
            log0, log1 = model(xb)
            l0 = criterion0(log0, yb_cum)
            l1 = criterion1(log1, yb_cum)
            m0 = (gb == 0).unsqueeze(1).float()
            m1 = (gb == 1).unsqueeze(1).float()
            loss0 = (l0 * m0).sum() / (m0.sum() * l0.shape[1] + eps)
            loss1 = (l1 * m1).sum() / (m1.sum() * l1.shape[1] + eps)
            loss = 0.0
            if torch.isfinite(loss0):
                loss = loss + loss0
            if torch.isfinite(loss1):
                loss = loss + loss1
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if isinstance(loss, torch.Tensor):
                epoch_loss += loss.item()
        print(f"[FULL] ep{epoch+1}/{final_epochs} loss={epoch_loss:.4f}")

    os.makedirs(final_outdir, exist_ok=True)
    torch.save(model.state_dict(), os.path.join(final_outdir, "model.pth"))
    joblib.dump(pre_full, os.path.join(final_outdir, "preprocessor.pkl"))

    model.eval()
    with torch.no_grad():
        xb = torch.tensor(X_full, dtype=torch.float32, device=device)
        lo0, lo1 = model(xb)
        p0_cum = torch.sigmoid(lo0).cpu().numpy()
        p1_cum = torch.sigmoid(lo1).cpu().numpy()
        p0_full = unpack_class_probs_from_cumulative(p0_cum)
        p1_full = unpack_class_probs_from_cumulative(p1_cum)

    mask_b1 = (g_full == "BRCA1")
    mask_b2 = (g_full == "BRCA2")
    thr0 = float(thresholds_final["BRCA1"])
    thr1 = float(thresholds_final["BRCA2"])

    rep = {
        "final_train_epochs": int(final_epochs),
        "thresholds": {"BRCA1": thr0, "BRCA2": thr1},
    }

    y_full_patho = (y_full == 3).astype(np.int64)
    p0_patho = p0_full[:, 3]
    p1_patho = p1_full[:, 3]

    if mask_b1.any():
        rep["BRCA1"] = compute_binary_metrics(
            y_full_patho[mask_b1], p0_patho[mask_b1], thr0
        )
    if mask_b2.any():
        rep["BRCA2"] = compute_binary_metrics(
            y_full_patho[mask_b2], p1_patho[mask_b2], thr1
        )

    rep_heads_final = head_reports(
        y_full_patho, p0_patho, p1_patho, mask_b1, mask_b2, thr0, thr1
    )
    rep.update(
        {
            "BRCA1_confusion": rep_heads_final.get("BRCA1_confusion", {}),
            "BRCA2_confusion": rep_heads_final.get("BRCA2_confusion", {}),
        }
    )
    print(f"[FINAL] BRCA1 Confusion: {rep['BRCA1_confusion'].get('matrix')}")
    print(f"[FINAL] BRCA2 Confusion: {rep['BRCA2_confusion'].get('matrix')}")

    print("Generating evaluation graphs...")
    generate_evaluation_graphs(
        final_dir=final_outdir,
        y_full=y_full,
        g_full=g_full,
        p0_full=p0_full,
        p1_full=p1_full,
        thresholds=thresholds_final,
    )

    save_json(rep, os.path.join(final_outdir, "metrics.json"))
    save_json(rep, os.path.join(final_outdir, "cv_report.json"))
    return rep

def main():
    cfg = Config()
    train_file_path = Path(cfg.paths.train_csv_path)

    if train_file_path.exists():
        df = load_training_csv(cfg.paths.train_csv_path)
        print(f"Loaded training dataset from {train_file_path}")
    else:
        print("Creating training dataset before training the model !")
        inp = "./files/clinvar_BRCA1_BRCA2_GRCh38_merged.csv.gz"
        df = build_training_dataset(
            input_csv=inp,
            output_csv=train_file_path
        )

    strat = (df[cfg.cols.gene_col].astype(str) + "_" + df[cfg.cols.label_col].astype(int).astype(str))
    y_all, g_all, idx_all, used_cols, (num_cols, bin_cols, cat_cols) = prepare_feature_sets(df, cfg.cols)

    skf = StratifiedKFold(
        n_splits=cfg.cv.n_splits,
        shuffle=cfg.cv.shuffle,
        random_state=cfg.cv.random_state
    )

    outdir = Path(cfg.paths.artifacts_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    fold_reports = []
    folds_info = []

    for fold_idx, (tr_idx, va_idx) in enumerate(skf.split(df, strat)):
        print(f"==== Iniciando Fold {fold_idx} ====")
        print(f"Treino: {len(tr_idx)} exemplos, Validação: {len(va_idx)} exemplos")

        df_tr = df.iloc[tr_idx].copy()
        pre_fold = Preprocessor(num_cols, bin_cols, cat_cols)
        pre_fold.fit(df_tr)

        X_tr = pre_fold.transform(df.iloc[tr_idx])
        X_va = pre_fold.transform(df.iloc[va_idx])

        y_tr, y_va = y_all[tr_idx], y_all[va_idx]
        g_tr = g_all.iloc[tr_idx].to_numpy()
        g_va = g_all.iloc[va_idx].to_numpy()

        fold_dir = os.path.join(cfg.paths.artifacts_dir, f"fold_{fold_idx}")
        os.makedirs(fold_dir, exist_ok=True)
        joblib.dump(pre_fold, os.path.join(fold_dir, "preprocessor.pkl"))

        if cfg.flags.use_mlp:
            report = train_one_fold(cfg, fold_idx, X_tr, y_tr, g_tr, X_va, y_va, g_va, X_tr.shape[1], pre_fold, fold_dir)
            fold_reports.append(report)
            folds_info.append({
                "fold": fold_idx,
                "score": report.get("mean_val_auroc", float("-inf")),
                "dir": fold_dir
            })

    agg = {}
    for gene in ["BRCA1","BRCA2"]:
        keys = ["auroc","auprc","recall","precision","f1","balanced_acc","mcc","accuracy"]
        for k in keys:
            vals = []
            for rep in fold_reports:
                if gene in rep and k in rep[gene]:
                    vals.append(rep[gene][k])
            if vals:
                agg[f"{gene}_{k}_mean"] = float(np.nanmean(vals))
                agg[f"{gene}_{k}_std"]  = float(np.nanstd(vals))

    if folds_info:
        winner = max(folds_info, key=lambda x: x["score"])
        winner_fold = int(winner["fold"])
        winner_dir = winner["dir"]
        agg["winner_fold"] = {"fold": winner_fold, "score_name": "mean_val_auroc", "score_value": float(winner["score"])}

        metrics_winner_path = os.path.join(winner_dir, f"metrics_fold{winner_fold}.json")
        print(f"Fold vencedor: {winner_fold}")
        with open(metrics_winner_path, "r") as f:
            metrics_winner = json.load(f)

        final_thresholds = {
            "BRCA1": float(metrics_winner["thresholds"]["BRCA1"]),
            "BRCA2": float(metrics_winner["thresholds"]["BRCA2"]),
        }
        agg["final_thresholds"] = final_thresholds

        if fold_reports:
            best_epochs = [int(r.get("best_epoch", 0)) for r in fold_reports]
            final_epochs = int(np.median(best_epochs)) + 1
        else:
            final_epochs = cfg.train.epochs

        pre_full = Preprocessor(num_cols, bin_cols, cat_cols)
        pre_full.fit(df)
        X_full = pre_full.transform(df)
        y_full = y_all
        g_full = g_all.to_numpy()
        init_state = torch.load(os.path.join(winner_dir, f"mlp_fold{winner_fold}.pth"), map_location=device)
        final_outdir = os.path.join(cfg.paths.artifacts_dir, "final_model_")
        final_rep = train_full_data(
            cfg=cfg,
            X_full=X_full,
            y_full=y_full,
            g_full=g_full,
            init_state_dict=init_state,
            pre_full=pre_full,
            final_epochs=final_epochs,
            thresholds_final=final_thresholds,
            final_outdir=final_outdir
        )
        agg["final_model"] = {
            "path_model": os.path.join(final_outdir, "model.pth"),
            "path_preprocessor": os.path.join(final_outdir, "preprocessor.pkl"),
            "path_metrics": os.path.join(final_outdir, "metrics.json"),
            "train_epochs": int(final_rep.get("final_train_epochs", final_epochs)),
            "thresholds": final_rep.get("thresholds", {}),
            "per_head": {
                "BRCA1": final_rep.get("BRCA1", {}),
                "BRCA2": final_rep.get("BRCA2", {}),
            },
            "confusion": {
                "BRCA1_confusion": final_rep.get("BRCA1_confusion", {}),
                "BRCA2_confusion": final_rep.get("BRCA2_confusion", {}),
            }
        }

    save_json(agg, os.path.join(cfg.paths.artifacts_dir, "cv_report.json"))
    print("Done. Summary:", agg)

if __name__ == "__main__":
    main()
