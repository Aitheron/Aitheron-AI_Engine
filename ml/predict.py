import os
import json
import joblib
import numpy as np
import pandas as pd
import torch

from ml.train._settings import Config
from ml.train.model import MultitaskMLP

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def _load_thresholds(dir_path: str) -> tuple[float, float]:
    name = next(f for f in os.listdir(dir_path) if f.startswith("metrics") and f.endswith(".json"))
    with open(os.path.join(dir_path, name), "r") as f:
        rep = json.load(f)
    return float(rep["thresholds"]["BRCA1"]), float(rep["thresholds"]["BRCA2"])

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
                n_heads=2
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
        genes = df_in[gene_col].to_numpy()

        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            xb = torch.tensor(X[start:end], dtype=torch.float32, device=device)
            log0, log1 = self.model(xb)
            p0 = torch.sigmoid(log0).squeeze(1).cpu().numpy()
            p1 = torch.sigmoid(log1).squeeze(1).cpu().numpy()

            g = genes[start:end]
            m1 = (g == "BRCA1")
            m2 = (g == "BRCA2")

            probs[start:end][m1] = p0[m1]
            probs[start:end][m2] = p1[m2]

            preds[start:end][m1] = (probs[start:end][m1] >= self.thr_b1).astype(int)
            preds[start:end][m2] = (probs[start:end][m2] >= self.thr_b2).astype(int)

        out = df_in.copy()
        final_preds = []
        for pred in preds:
            if pred == 1:
                pred = "Patogênico"
            else:
                pred = "Benigno"
            final_preds.append(pred)
        out["Classificação Clinica"] = final_preds            # 1=patogênico, 0=benigno
        col = out.pop("Classificação Clinica")
        i = out.columns.get_loc("impact")  # posição da coluna "Impacto"
        out.insert(i + 1, "Classificação Clinica", col)
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
