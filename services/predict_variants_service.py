import pandas as pd
import numpy as np

from data_prep import (
    detect_variants_from_fasta,
    enrich_patient_variants_with_vep
)

from services._utils import drop_cols, reorder_columns
from ml import get_predictor


def run_predict_variants_service(gene, file_path):

    initial_df = detect_variants_from_fasta(
        target_gene=gene,
        fasta_file_path=file_path
    )

    with_vep_df = enrich_patient_variants_with_vep(
        df=initial_df
    )

    to_predict_df = with_vep_df.copy()
    new_values = [gene] * len(to_predict_df)
    to_predict_df["GeneSymbol"] = new_values

    predictor = get_predictor()
    predicted_values = predictor.predict_df(to_predict_df=to_predict_df)
    drop_cols(
        df=predicted_values,
        keep=["IsCoding"],
        inplace=True,
        also_drop=[
            "ImpactRank",
            "len_ref",
            "len_alt",
            "length_change"
        ]
    )

    cols_pct = [
    "Confiança do Modelo (%)",
    "Risco da mutação ser Patogênica (%)",
]

    predicted_values[cols_pct] = (
        predicted_values[cols_pct]
        .apply(pd.to_numeric, errors="coerce")
        .replace([np.inf, -np.inf], np.nan)
        .mul(100)
        .round(2)
    )
    
    moves = [
        ("Classificação Clinica", "impact"),
        ("Confiança do Modelo (%)", "Classificação Clinica"),
        ("Risco da mutação ser Patogênica (%)", "Confiança do Modelo (%)"),
    ]
    predicted_values = reorder_columns(
        df=predicted_values,
        moves=moves,
        after=True
    )
    predicted_values["IsCoding"] = predicted_values["IsCoding"].map({
        1: "Sim", 0: "Não"
    }).fillna("Não")

    predicted_values.rename(columns={
        "impact": "Impacto",
        "ConsequenceTerms": "Efeito no transcrito",
        "IsCoding" : "Ocorreu em região codificante ?",
        "GeneSymbol" : "Gene",
        "Type" : "Tipo de Mutação"
    }, inplace=True)

    predicted_values.to_csv("./teste.csv", index=False)
    return predicted_values

if __name__ == "__main__":
    run_predict_variants_service(
        gene="BRCA1",
        file_path='./patient.fasta'
    )