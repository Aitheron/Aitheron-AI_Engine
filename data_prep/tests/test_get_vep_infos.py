import pandas as pd
import data_prep.get_vep_infos as mod

def test_normalize_rs_variants():
    assert mod.normalize_rs("rs123") == "rs123"
    assert mod.normalize_rs("123") == "rs123"
    assert mod.normalize_rs("-1") is None
    assert mod.normalize_rs(None) is None

def test_extract_ref_alt_from_rec_paths():
    rec = {"allele_string": "A/T"}
    assert mod.extract_ref_alt_from_rec(rec, "rs1") == ("A","T")

    rec2 = {"colocated_variants":[{"id":"rs2","allele_string":"G/C"}]}
    assert mod.extract_ref_alt_from_rec(rec2, "rs2") == ("G","C")

    rec3 = {"colocated_variants":[{"allele_string":"T/A"}]}
    assert mod.extract_ref_alt_from_rec(rec3, "rsX") == ("T","A")

    assert mod.extract_ref_alt_from_rec({}, "rsZ") == (None, None)

def test_choose_transcript_prefers_protein_and_impact():
    t1 = {"impact":"LOW","protein_start":10}
    t2 = {"impact":"HIGH","protein_start":20}
    chosen, basis = mod._choose_transcript([t1,t2])
    assert chosen is t2
    assert basis in ("protein_fields","impact")

def test_update_df_inplace_with_vep_uses_process_item_and_writes(tmp_path, monkeypatch):
    # cria df mínimo
    df = pd.DataFrame({"RS# (dbSNP)": ["rs1","-1","rs2"]})
    out = tmp_path / "out.csv"

    # mock process_item devolvendo 2 entradas válidas
    def fake_process_item(args):
        idx, raw, session = args
        if raw == "-1":
            return idx, None, None, {}
        return idx, "A", "T", {
            "IsCoding": True,
            "ENST": "ENSTX",
            "ENSP": "ENSPY",
            "ConsequenceTerms": "missense_variant",
            "Impact": "MODERATE",
            "ProteinStart": 5,
            "ProteinEnd": 6
        }
    monkeypatch.setattr("data_prep.get_vep_infos.process_item", fake_process_item)

    mod.update_df_inplace_with_vep(df, str(out), max_workers=2)
    got = pd.read_csv(out, dtype=str)
    assert "ReferenceAllele" in got.columns and "AlternateAllele" in got.columns
    assert got["ReferenceAllele"].notna().sum() >= 2
    assert got["IsCoding"].notna().sum() >= 2
