import pandas as pd
import data_prep.prepare_training_dataset as mod

def test_build_training_dataset_basic(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "IMPACT_RANK", {"HIGH": 4, "MODERATE": 3})
    monkeypatch.setattr(mod, "TERM_TO_FLAG", {"missense_variant":"IsMissense","frameshift_variant":"IsFrameshift"})
    monkeypatch.setattr(mod, "TYPE_FLAG_COLS", ["IsSNV","IsDeletion","IsDuplication","IsInsertion","IsIndel"])

    inp = tmp_path / "in.csv"
    df = pd.DataFrame({
        "GeneSymbol":["BRCA1","BRCA1"],
        "Stop":["5","6"],
        "ReferenceAllele":["A","AA"],
        "AlternateAllele":["T","-"],
        "Type":["single nucleotide variant","deletion"],
        "ConsequenceTerms":["missense_variant","frameshift_variant"],
        "Impact":["MODERATE","HIGH"],
        "IsCoding":[1,0],
        "VariantAllele":["T",""],
        "ClinicalSignificance":["Pathogenic","Benign"]
    })
    df.to_csv(inp, index=False)

    out = tmp_path / "out.csv"
    res = mod.build_training_dataset(str(inp), str(out))
    assert out.exists()
    assert "Label" in res.columns
    assert "ImpactRank" in res.columns
    assert {"IsMissense","IsFrameshift"}.issubset(set(res.columns))
    assert {"IsSNV","IsDeletion","IsDuplication","IsInsertion","IsIndel"}.issubset(set(res.columns))
    assert "VariantAllele" not in res.columns
    assert "ClinicalSignificance" not in res.columns
    # checa renome de Stop->End
    assert "End" in res.columns and "Stop" not in res.columns
