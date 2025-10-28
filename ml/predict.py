import os
import json
import joblib
import numpy as np
import pandas as pd
import torch

from ml.train._settings import Config
from ml.train.model import MultitaskMLP

from ml._utils import (
    entropy,
    unpack_class_probs_from_cumulative
)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def _load_thresholds(dir_path: str) -> tuple[float, float]:
    name = next(f for f in os.listdir(dir_path) if f.startswith("metrics") and f.endswith(".json"))
    with open(os.path.join(dir_path, name), "r") as f:
        rep = json.load(f)
    return float(rep["thresholds"]["BRCA1"]), float(rep["thresholds"]["BRCA2"])

def _top_stats(p):
    top1 = p.max(axis=1)
    idx1 = p.argmax(axis=1)
    p2 = p.copy()
    rows = np.arange(len(p2))
    p2[rows, idx1] = -1.0
    top2 = p2.max(axis=1)
    margin = top1 - top2
    return top1, idx1, top2, margin

class Predictor:
    # carrega preprocessor, thresholds e pesos; instancia o modelo no primeiro predict (para pegar o input_dim correto)
    def __init__(self, final_model_path: str, cfg: Config):
        self.cfg = cfg
        self.final_model_path = final_model_path
        self.fold_dir = os.path.dirname(final_model_path)

        pp_path = os.path.join(self.fold_dir, "preprocessor.pkl")
        if not os.path.exists(pp_path):
            raise FileNotFoundError(f"Preprocessor não encontrado em {pp_path}")
        self.pre = joblib.load(pp_path)

        self.thr_b1, self.thr_b2 = _load_thresholds(self.fold_dir)

        self._state_dict = torch.load(self.final_model_path, map_location=device)
        self.model = None  # definido no primeiro predict

    def _ensure_model(self, input_dim: int):
        if self.model is None:
            self.model = MultitaskMLP(
                in_dim=input_dim,
                hidden_dims=self.cfg.train.hidden_dims,
                dropout=self.cfg.train.dropout,
                n_heads=2,
                n_classes=4
            ).to(device)
            self.model.load_state_dict(self._state_dict)
            self.model.eval()

    @torch.no_grad()
    def predict_df(self, df_in: pd.DataFrame, batch_size: int = 4096) -> pd.DataFrame:
        gene_col = self.cfg.cols.gene_col
        if gene_col not in df_in.columns:
            raise ValueError(f"Coluna obrigatória ausente: {gene_col}")

        X = self.pre.transform(df_in)
        self._ensure_model(X.shape[1])

        n = len(df_in)
        preds = np.zeros(n, dtype=int)
        probs = np.zeros(n, dtype=float)
        class_idx = np.zeros(n, dtype=int)
        conf = np.zeros(n, dtype=float)
        margin = np.zeros(n, dtype=float)
        top1_arr = np.zeros(n, dtype=float)
        top2_arr = np.zeros(n, dtype=float)
        p0_store = np.zeros((n, 4), dtype=float)
        p1_store = np.zeros((n, 4), dtype=float)
        genes = df_in[gene_col].to_numpy()

        LABELS = {0: "Benigno", 1: "Possivelmente Benigno", 2: "VUS", 3: "Patogênico"}

        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            xb = torch.tensor(X[start:end], dtype=torch.float32, device=device)
            log0, log1 = self.model(xb)
            p0_cum = torch.sigmoid(log0).cpu().numpy()
            p1_cum = torch.sigmoid(log1).cpu().numpy()
            p0_full = unpack_class_probs_from_cumulative(p0_cum)
            p1_full = unpack_class_probs_from_cumulative(p1_cum)
            k0 = p0_full.argmax(axis=1)
            k1 = p1_full.argmax(axis=1)

            top1_0, _, top2_0, margin_0 = _top_stats(p0_full)
            top1_1, _, top2_1, margin_1 = _top_stats(p1_full)

            H0 = entropy(p0_full); H1 = entropy(p1_full)
            H0n = H0 / np.log(4.0); H1n = H1 / np.log(4.0)
            conf0 = 1.0 - H0n
            conf1 = 1.0 - H1n

            p0_patho = p0_full[:, 3]
            p1_patho = p1_full[:, 3]

            g = genes[start:end]
            m1 = (g == "BRCA1")
            m2 = (g == "BRCA2")

            probs[start:end][m1] = p0_patho[m1]
            probs[start:end][m2] = p1_patho[m2]

            preds[start:end][m1] = (probs[start:end][m1] >= self.thr_b1).astype(int)
            preds[start:end][m2] = (probs[start:end][m2] >= self.thr_b2).astype(int)

            class_idx[start:end][m1] = k0[m1]
            class_idx[start:end][m2] = k1[m2]

            conf[start:end][m1] = conf0[m1]
            conf[start:end][m2] = conf1[m2]
            margin[start:end][m1] = margin_0[m1]
            margin[start:end][m2] = margin_1[m2]
            top1_arr[start:end][m1] = top1_0[m1]
            top1_arr[start:end][m2] = top1_1[m2]
            top2_arr[start:end][m1] = top2_0[m1]
            top2_arr[start:end][m2] = top2_1[m2]
            p0_store[start:end][m1] = p0_full[m1]
            p1_store[start:end][m2] = p1_full[m2]

            for j in range(end - start):
                if m1[j]:
                    p = p0_full[j]
                    k = int(k0[j])
                    print(f"BRCA1 idx {start+j:>5} | P0={p[0]:.4f} P1={p[1]:.4f} P2={p[2]:.4f} P3={p[3]:.4f} | argmax={k} ({LABELS[k]}) | top1={top1_0[j]:.4f} top2={top2_0[j]:.4f} margin={margin_0[j]:.4f} conf={conf0[j]:.4f}")
                elif m2[j]:
                    p = p1_full[j]
                    k = int(k1[j])
                    print(f"BRCA2 idx {start+j:>5} | P0={p[0]:.4f} P1={p[1]:.4f} P2={p[2]:.4f} P3={p[3]:.4f} | argmax={k} ({LABELS[k]}) | top1={top1_1[j]:.4f} top2={top2_1[j]:.4f} margin={margin_1[j]:.4f} conf={conf1[j]:.4f}")

        out = df_in.copy()
        label_map = {0: "Benigno", 1: "Possivelmente Benigno", 2: "VUS", 3: "Patogênico"}
        final_class = [label_map[int(k)] for k in class_idx]
        out["Classificação Clinica"] = final_class
        col = out.pop("Classificação Clinica")
        i = out.columns.get_loc("impact")
        out.insert(i + 1, "Classificação Clinica", col)

        out["Confiança do Modelo"] = np.round(conf, 4)
        out["Risco da mutação ser Patogênica"] = np.round(probs, 4)

        out.to_csv("predicted.csv", index=False)
        return out

if __name__ == "__main__":
    cfg = Config()
    predictor = Predictor("./artifacts/final_model/model.pth", cfg)
    
    df = pd.read_csv("./patient_variants_annotated.csv", sep=",")
    to_predict_df = df.copy()
    new_values = ["BRCA1"] * len(to_predict_df)
    to_predict_df["GeneSymbol"] = new_values
    to_predict_df.drop(
        columns=[
            "Origin",
            "Type",
            "ConsequenceTerms",
            "impact"
        ]
    )

    df_out = predictor.predict_df(to_predict_df, batch_size=4096)
