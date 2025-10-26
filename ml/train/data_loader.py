import os
import joblib
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer

class Preprocessor:
    def __init__(self, numeric_cols, binary_cols, categorical_cols):
        self.numeric_cols = numeric_cols
        self.binary_cols = binary_cols
        self.categorical_cols = categorical_cols
        self.num_imputer = SimpleImputer(strategy="constant", fill_value=0.0)
        self.num_scaler = StandardScaler(with_mean=True, with_std=True)
        self.ohe = OneHotEncoder(handle_unknown="ignore", sparse_output=False)

    def fit(self, df: pd.DataFrame):
        # conversão numérica segura ("" -> NaN)
        X_num = df[self.numeric_cols].apply(pd.to_numeric, errors="coerce")
        X_num = self.num_imputer.fit_transform(X_num)
        self.num_scaler.fit(X_num)

        if self.categorical_cols:
            X_cat = df[self.categorical_cols].astype(str)
            self.ohe.fit(X_cat)
        return self

    def transform(self, df: pd.DataFrame):
        X_num = df[self.numeric_cols].apply(pd.to_numeric, errors="coerce")
        X_num = self.num_imputer.transform(X_num)
        X_num = self.num_scaler.transform(X_num)

        # binárias já são 0/1 no CSV final
        X_bin = df[self.binary_cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()

        if self.categorical_cols:
            X_cat = df[self.categorical_cols].astype(str)
            X_cat = self.ohe.transform(X_cat)
            X = np.hstack([X_num, X_bin, X_cat])
        else:
            X = np.hstack([X_num, X_bin])
        return X

    def save(self, path: str):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        joblib.dump({
            "numeric_cols": self.numeric_cols,
            "binary_cols": self.binary_cols,
            "categorical_cols": self.categorical_cols,
            "num_imputer": self.num_imputer,
            "num_scaler": self.num_scaler,
            "ohe": self.ohe
        }, path)

    @classmethod
    def load(cls, path: str):
        payload = joblib.load(path)
        obj = cls(
            payload["numeric_cols"],
            payload["binary_cols"],
            payload["categorical_cols"]
        )
        obj.num_imputer = payload["num_imputer"]
        obj.num_scaler = payload["num_scaler"]
        obj.ohe = payload["ohe"]
        return obj

def load_training_csv(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    return df

def prepare_feature_sets(df: pd.DataFrame, cols_cfg):
    # todas as flags binárias Is* (mais IsCoding explicitamente)
    binary_cols = [col for col in df.columns if col.startswith(cols_cfg.is_prefix)]
    if "IsCoding" in df.columns and "IsCoding" not in binary_cols:
        binary_cols.append("IsCoding")

    numeric_cols = [col for col in cols_cfg.numeric_cols if col in df.columns]
    categorical_cols = [col for col in cols_cfg.cat_cols if col in df.columns]

    # alvo e “cabeça” (gene)
    y = df[cols_cfg.label_col].astype(int)
    gene = df[cols_cfg.gene_col].astype(str)

    used_feature_names = numeric_cols + binary_cols + categorical_cols

    # NÃO faz fit/transform aqui para evitar leakage
    return y.to_numpy(), gene, df.index, used_feature_names, (numeric_cols, binary_cols, categorical_cols)

