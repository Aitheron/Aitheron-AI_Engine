import pandas as pd
import numpy as np

import services.predict_variants_service as psvc

class DummyPredictor:
    def __init__(self, out_df):
        self.out_df = out_df
        self.called = False

    def predict_df(self, to_predict_df):
        # valida que GeneSymbol foi injetado
        assert (to_predict_df["GeneSymbol"] == to_predict_df["GeneSymbol"].iloc[0]).all()
        self.called = True
        return self.out_df.copy()


def test_predict_happy_path_formats_output(monkeypatch, tmp_path):
    # detect_variants_from_fasta -> df mínimo
    def fake_detect_variants_from_fasta(target_gene, fasta_file_path):
        return pd.DataFrame({"pos": [1, 2], "Type": ["SNV", "SNV"]})

    # enrich_patient_variants_with_vep -> adiciona colunas necessárias
    def fake_enrich_patient_variants_with_vep(df):
        out = df.copy()
        out["ConsequenceTerms"] = ["missense_variant", "synonymous_variant"]
        out["IsCoding"] = [1, 0]
        out["impact"] = ["Pathogenic", "Benign"]
        return out

    # predictor -> retorna probabilidades em 0..1
    pred_out = pd.DataFrame({
        "pos": [1, 2],
        "Type": ["SNV", "SNV"],
        "ConsequenceTerms": ["missense_variant", "synonymous_variant"],
        "IsCoding": [1, 0],
        "impact": ["Pathogenic", "Benign"],
        "Confiança do Modelo (%)": [0.9, 0.42],
        "Risco da mutação ser Patogênica (%)": [0.95, 0.05],
        "GeneSymbol": ["WILL_BE_OVERWRITTEN", "WILL_BE_OVERWRITTEN"],
    })
    predictor = DummyPredictor(pred_out)

    # mocks
    monkeypatch.setattr("services.predict_variants_service.detect_variants_from_fasta", fake_detect_variants_from_fasta)
    monkeypatch.setattr("services.predict_variants_service.enrich_patient_variants_with_vep", fake_enrich_patient_variants_with_vep)
    monkeypatch.setattr("services.predict_variants_service.get_predictor", lambda: predictor)

    # drop_cols e reorder_columns como no-op
    monkeypatch.setattr("services.predict_variants_service.drop_cols", lambda df, keep, inplace, also_drop: df)
    monkeypatch.setattr("services.predict_variants_service.reorder_columns", lambda df, moves, after: df)

    df = psvc.run_predict_variants_service(gene="BRCA1", file_path=str(tmp_path / "patient.fasta"))

    # predictor foi acionado
    assert predictor.called is True

    # porcentagens devem estar em 0..100 com 2 casas
    assert set(df["Confiança do Modelo (%)"].round(2)) == {90.0, 42.0}
    assert set(df["Risco da mutação ser Patogênica (%)"].round(2)) == {95.0, 5.0}

    # mapeamento de IsCoding para Sim/Não
    assert set(df["Ocorreu em região codificante ?"]) == {"Sim", "Não"}

    # renomeações aplicadas
    for col in ["Impacto", "Efeito no transcrito", "Gene", "Tipo de Mutação"]:
        assert col in df.columns


def test_predict_handles_inf_values(monkeypatch, tmp_path):
    def fake_detect_variants_from_fasta(target_gene, fasta_file_path):
        return pd.DataFrame({"pos": [1]})

    def fake_enrich_patient_variants_with_vep(df):
        return pd.DataFrame({
            "pos": [1],
            "ConsequenceTerms": ["missense_variant"],
            "IsCoding": [1],
            "impact": ["Pathogenic"],
            "Type": ["SNV"],
        })

    out = pd.DataFrame({
        "pos": [1],
        "ConsequenceTerms": ["missense_variant"],
        "IsCoding": [1],
        "impact": ["Pathogenic"],
        "Type": ["SNV"],
        "Confiança do Modelo (%)": [np.inf],
        "Risco da mutação ser Patogênica (%)": [-np.inf],
        "GeneSymbol": ["X"],
    })
    predictor = DummyPredictor(out)

    monkeypatch.setattr("services.predict_variants_service.detect_variants_from_fasta", fake_detect_variants_from_fasta)
    monkeypatch.setattr("services.predict_variants_service.enrich_patient_variants_with_vep", fake_enrich_patient_variants_with_vep)
    monkeypatch.setattr("services.predict_variants_service.get_predictor", lambda: predictor)
    monkeypatch.setattr("services.predict_variants_service.drop_cols", lambda df, keep, inplace, also_drop: df)
    monkeypatch.setattr("services.predict_variants_service.reorder_columns", lambda df, moves, after: df)

    df = psvc.run_predict_variants_service(gene="BRCA2", file_path=str(tmp_path / "patient2.fasta"))
    # após replace de inf/-inf e multiplicação, vira NaN
    assert df["Confiança do Modelo (%)"].isna().iloc[0]
    assert df["Risco da mutação ser Patogênica (%)"].isna().iloc[0]
