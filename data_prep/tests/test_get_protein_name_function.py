import pandas as pd
import data_prep.get_protein_name_function as mod

def test_normalize_target_and_choose_best():
    assert mod._normalize_target("P38398") == {"primaryAccession": "P38398"}
    t = {"uniProtkbId": "X1"}
    assert mod._normalize_target(t)["primaryAccession"] == "X1"

    a = {"primaryAccession": "A", "reviewed": True}
    b = {"primaryAccession": "B", "reviewed": False}
    pick = mod._choose_best_target([b, a], restrict_taxon=None)
    assert pick["primaryAccession"] == "A"

def test_extract_name_function_features_basic():
    entry = {
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Prot X"}}},
        "comments": [{"commentType": "FUNCTION", "texts": [{"value": "Func 1"}, {"value": "Func 2"}]}],
        "features": [
            {"type": "DOMAIN", "description": "desc", "location": {"start": {"value": 1}, "end": {"value": 10}}}
        ]
    }
    name, func, feats = mod.extract_name_function_features(entry)
    assert "Prot X" in str(name)
    assert "Func 1" in func and "Func 2" in func
    assert "DOMAIN" in feats and "[1-10]" in feats

def test_map_and_enrich_with_uniprot(tmp_path, monkeypatch):
    df = pd.DataFrame({"ENSP": ["ENSP1",""], "ENST": ["","ENST2"]})

    # mock mapeamentos e detalhes
    monkeypatch.setattr("data_prep.get_protein_name_function.map_ids_to_uniprot",
                        lambda from_db, ids: {"ENSP1":"P11111"} if "Protein" in from_db else {"ENST2":"Q22222"})
    details = {"P11111": {"ProteinDesc":"D1","ProteinFunction":"F1","ProteinFeatures":"Fe1"},
               "Q22222": {"ProteinDesc":"D2","ProteinFunction":"F2","ProteinFeatures":"Fe2"}}
    monkeypatch.setattr("data_prep.get_protein_name_function.fetch_details_for_accessions",
                        lambda accessions, max_workers=20: details)

    out = mod.enrich_with_uniprot(df, output_path=str(tmp_path/"out.csv"), ensp_col="ENSP", enst_col="ENST")
    assert (out["ProteinDesc"].notna()).any()
    assert (out["ProteinFunction"].notna()).any()
    assert (out["ProteinFeatures"].notna()).any()
