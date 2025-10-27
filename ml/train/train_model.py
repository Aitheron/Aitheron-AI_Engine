import os
import math
import numpy as np
import pandas as pd
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

from sklearn.model_selection import StratifiedKFold

from ml.train._settings import Config
from ml.train.data_loader import load_training_csv, prepare_feature_sets, Preprocessor
from ml.train.model import MultitaskMLP
from ml.train.metrics import compute_binary_metrics, pick_best_threshold_by_f1, pick_threshold_for_recall, pick_threshold_under_recall, head_reports
from ml.train._utils import save_json, save_torch_model

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def make_mask_for_heads(gene_series: pd.Series) -> dict[str, np.ndarray]:
    # mapeia BRCA1 -> head 0, BRCA2 -> head 1
    genes = gene_series.to_numpy()
    mask0 = (genes == "BRCA1").astype(np.float32)
    mask1 = (genes == "BRCA2").astype(np.float32)
    return {"BRCA1": mask0, "BRCA2": mask1}

def train_one_fold(cfg: Config, fold_idx: int, X_train, y_train, g_train, X_val, y_val, g_val, input_dim: int, pre: Preprocessor, outdir_fold: str) -> dict:
    model = MultitaskMLP(
        in_dim=input_dim,
        hidden_dims=cfg.train.hidden_dims,
        dropout=cfg.train.dropout,
        n_heads=2
    ).to(device)

    # class weighting por cabeça (pos_weight para patogênicos)
    def pos_weight_for_head(y, g, gene_name):
        mask = (g == gene_name)
        y_head = y[mask]
        if y_head.sum() == 0 or y_head.sum() == len(y_head):
            return torch.tensor(1.0, device=device)
        neg = (y_head == 0).sum()
        pos = (y_head == 1).sum()
        return torch.tensor(neg / max(pos, 1), dtype=torch.float32, device=device)

    pw0 = pos_weight_for_head(y_train, g_train, "BRCA1")
    pw1 = pos_weight_for_head(y_train, g_train, "BRCA2")

    # reduction='none' para mascarar por lote
    criterion0 = nn.BCEWithLogitsLoss(pos_weight=pw0, reduction="none")
    criterion1 = nn.BCEWithLogitsLoss(pos_weight=pw1, reduction="none")

    optimizer = torch.optim.AdamW(model.parameters(), lr=cfg.train.lr, weight_decay=cfg.train.weight_decay)

    # inclui o gene no dataset para mascarar por lote (0=BRCA1, 1=BRCA2)
    ds_train = TensorDataset(
        torch.tensor(X_train, dtype=torch.float32),
        torch.tensor(y_train, dtype=torch.float32).unsqueeze(1),
        torch.tensor((g_train == "BRCA2").astype(np.int64))
    )
    dl_train = DataLoader(ds_train, batch_size=cfg.train.batch_size, shuffle=True, drop_last=False)

    ds_val = TensorDataset(
        torch.tensor(X_val, dtype=torch.float32),
        torch.tensor(y_val, dtype=torch.float32).unsqueeze(1),
    )
    dl_val = DataLoader(ds_val, batch_size=1024, shuffle=False, drop_last=False)

    best_score = -math.inf
    best_state = None
    patience = cfg.train.early_stopping_patience
    no_improve = 0

    # masks para heads (validação no conjunto inteiro, sem shuffle)
    mask_va = make_mask_for_heads(pd.Series(g_val))
    m0_va_np = mask_va["BRCA1"]; m1_va_np = mask_va["BRCA2"]

    eps = 1e-9

    for epoch in range(cfg.train.epochs):
        model.train()
        epoch_loss = 0.0
        seen_brca1 = 0
        seen_brca2 = 0
        for xb, yb, gb in dl_train:
            xb = xb.to(device); yb = yb.to(device); gb = gb.to(device)

            logits0, logits1 = model(xb)

            # perdas elemento a elemento
            l0 = criterion0(logits0, yb)  # (B,1)
            l1 = criterion1(logits1, yb)  # (B,1)

            # aplica loss apenas nos exemplos da cabeça correspondente (por lote)
            m0 = (gb == 0).unsqueeze(1).float()  # BRCA1
            m1 = (gb == 1).unsqueeze(1).float()  # BRCA2

            b0 = int(m0.sum().item())
            b1 = int(m1.sum().item())
            seen_brca1 += b0
            seen_brca2 += b1

            loss0 = (l0 * m0).sum() / (m0.sum() + eps)
            loss1 = (l1 * m1).sum() / (m1.sum() + eps)

            loss = 0.0
            if torch.isfinite(loss0): loss = loss + loss0
            if torch.isfinite(loss1): loss = loss + loss1

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if isinstance(loss, torch.Tensor):
                epoch_loss += loss.item()
        print(f"[Fold {fold_idx}] ep{epoch+1} seen -> BRCA1={seen_brca1} BRCA2={seen_brca2}")

        # validação
        model.eval()
        with torch.no_grad():
            xb_val = torch.tensor(X_val, dtype=torch.float32, device=device)
            logits0, logits1 = model(xb_val)
            prob0 = torch.sigmoid(logits0).squeeze(1).cpu().numpy()
            prob1 = torch.sigmoid(logits1).squeeze(1).cpu().numpy()

        # thresholds por cabeça (apenas onde a cabeça se aplica)
        mask_b1 = np.asarray(mask_va["BRCA1"], dtype=bool)
        if mask_b1.any():
            thr0, _ = pick_threshold_under_recall(
                y_val[mask_b1], prob0[mask_b1],
                recall_target=0.99,
                objective="precision"
            )
        thr1 = pick_best_threshold_by_f1(y_val[m1_va_np > 0.5], prob1[m1_va_np > 0.5]) if (m1_va_np > 0.5).any() else 0.5

        # métrica principal para early stopping: média AUROC das duas cabeças (onde aplicável)
        scores = []
        if (m0_va_np > 0.5).any():
            m0 = compute_binary_metrics(y_val[m0_va_np > 0.5], prob0[m0_va_np > 0.5], thr0)
            scores.append(m0["auroc"])
        if (m1_va_np > 0.5).any():
            m1 = compute_binary_metrics(y_val[m1_va_np > 0.5], prob1[m1_va_np > 0.5], thr1)
            scores.append(m1["auroc"])
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

    # salva
    os.makedirs(outdir_fold, exist_ok=True)
    final_model_path = os.path.join(outdir_fold, f"mlp_fold{fold_idx}.pth")
    torch.save(best_state["model"], final_model_path)

    # métricas finais no val com thresholds escolhidos
    model.load_state_dict(best_state["model"])
    model.eval()
    with torch.no_grad():
        xb_val = torch.tensor(X_val, dtype=torch.float32, device=device)
        logits0, logits1 = model(xb_val)
        prob0 = torch.sigmoid(logits0).squeeze(1).cpu().numpy()
        prob1 = torch.sigmoid(logits1).squeeze(1).cpu().numpy()

    thr0 = best_state["thr0"]; thr1 = best_state["thr1"]
    report = {"fold": fold_idx, "best_epoch": best_state["epoch"], "mean_val_auroc": best_score}

    if (m0_va_np > 0.5).any():
        report["BRCA1"] = compute_binary_metrics(y_val[m0_va_np > 0.5], prob0[m0_va_np > 0.5], thr0)
    if (m1_va_np > 0.5).any():
        report["BRCA2"] = compute_binary_metrics(y_val[m1_va_np > 0.5], prob1[m1_va_np > 0.5], thr1)

    report["thresholds"] = {"BRCA1": thr0, "BRCA2": thr1}

    # Calculo da matriz de confusão por Fold
    mask_b1_final = np.asarray(m0_va_np, dtype=bool)
    mask_b2_final = np.asarray(m1_va_np, dtype=bool)
    rep_heads = head_reports(y_val, prob0, prob1, mask_b1_final, mask_b2_final, thr0, thr1)
    report.update(rep_heads)

    save_json(report, os.path.join(outdir_fold, f"metrics_fold{fold_idx}.json"))
    return report

def main():
    cfg = Config()

    df = load_training_csv(cfg.paths.train_csv_path)

    # estratificação por (GeneSymbol, label)
    strat = (df[cfg.cols.gene_col].astype(str) + "_" + df[cfg.cols.label_col].astype(int).astype(str))

    y_all, g_all, idx_all, used_cols, (num_cols, bin_cols, cat_cols) = prepare_feature_sets(df, cfg.cols)

    skf = StratifiedKFold(
        n_splits=cfg.cv.n_splits,
        shuffle=cfg.cv.shuffle,
        random_state=cfg.cv.random_state
    )

    outdir = Path(cfg.paths.artifacts_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # salvar preprocessor global por fold (treinado no train-fold)
    fold_reports = []
    y_all, g_all, idx_all, used_cols, (num_cols, bin_cols, cat_cols) = prepare_feature_sets(df, cfg.cols)

    skf = StratifiedKFold(
        n_splits=cfg.cv.n_splits,
        shuffle=cfg.cv.shuffle,
        random_state=cfg.cv.random_state
    )

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

        if cfg.flags.use_mlp:
            report = train_one_fold(cfg, fold_idx, X_tr, y_tr, g_tr, X_va, y_va, g_va, X_tr.shape[1], pre_fold, fold_dir)
            fold_reports.append(report)

    # agregação simples das métricas
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

    # Visualização da matriz de confusão
    conf_tot = {
        "BRCA1": {"tn": 0, "fp": 0, "fn": 0, "tp": 0},
        "BRCA2": {"tn": 0, "fp": 0, "fn": 0, "tp": 0},
    }
    for rep in fold_reports:
        for gene in ["BRCA1", "BRCA2"]:
            key = f"{gene}_confusion"
            if key in rep:
                for k in ["tn", "fp", "fn", "tp"]:
                    conf_tot[gene][k] += int(rep[key][k])

    brca1_mat = [[conf_tot["BRCA1"]["tn"], conf_tot["BRCA1"]["fp"]],
                 [conf_tot["BRCA1"]["fn"], conf_tot["BRCA1"]["tp"]]]
    brca2_mat = [[conf_tot["BRCA2"]["tn"], conf_tot["BRCA2"]["fp"]],
                 [conf_tot["BRCA2"]["fn"], conf_tot["BRCA2"]["tp"]]]

    print(f"[FINAL] BRCA1 Confusion: {brca1_mat}")
    print(f"[FINAL] BRCA2 Confusion: {brca2_mat}")

    # coloca no relatório agregado
    agg["BRCA1_confusion_final"] = {
        "tn": conf_tot["BRCA1"]["tn"], "fp": conf_tot["BRCA1"]["fp"],
        "fn": conf_tot["BRCA1"]["fn"], "tp": conf_tot["BRCA1"]["tp"],
        "matrix": brca1_mat,
        "total_examples": int(sum(sum(row) for row in brca1_mat))
    }
    agg["BRCA2_confusion_final"] = {
        "tn": conf_tot["BRCA2"]["tn"], "fp": conf_tot["BRCA2"]["fp"],
        "fn": conf_tot["BRCA2"]["fn"], "tp": conf_tot["BRCA2"]["tp"],
        "matrix": brca2_mat,
        "total_examples": int(sum(sum(row) for row in brca2_mat))
    }

    save_json(agg, os.path.join(cfg.paths.artifacts_dir, "cv_report.json"))
    print("Done. Summary:", agg)

if __name__ == "__main__":
    main()
