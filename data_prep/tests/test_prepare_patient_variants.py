import pandas as pd
import data_prep.prepare_patient_variants as mod

def test_to_vcf_record_snv_deletion_insertion(monkeypatch):
    # SNV
    row = {"Chromosome":"17","Type":"SNV","Ref":"C","Alt":"A","Start":430}
    assert mod.to_vcf_record(row, {}) == ("17","C","A",430)

    # Deletion usa base âncora
    monkeypatch.setattr(mod, "fetch_base", lambda chrom, pos, cache: "G")
    row = {"Chromosome":"17","Type":"Deleção","Ref":"T","Alt":"-","Start":100}
    chrom, ref, alt, pos = mod.to_vcf_record(row, {})
    assert pos == 99 and ref.startswith("G") and alt == "G"

    # Insertion usa âncora
    row = {"Chromosome":"17","Type":"Inserção","Ref":"-","Alt":"TT","Start":200}
    chrom, ref, alt, pos = mod.to_vcf_record(row, {})
    assert ref == "G" and alt.startswith("G")

def test_apply_type_flags_and_decide_is_coding(monkeypatch):
    # patcha tabelas de settings para cenário mínimo
    monkeypatch.setattr(mod, "TERM_TO_FLAG", {"missense_variant":"IsMissense"})
    monkeypatch.setattr(mod, "TYPE_FLAG_COLS", ["IsSNV","IsDeletion","IsDuplication","IsInsertion","IsIndel"])

    flags = mod._apply_type_flags("SNV")
    assert flags["IsSNV"] == 1 and sum(flags.values()) == 1

    is_coding = mod.decide_is_coding(
        rec={"transcript_consequences":[]},
        most="missense_variant",
        chosen_tc={"biotype":"protein_coding"},
        observed_terms={"missense_variant"}
    )
    assert is_coding == 1

def test_enrich_patient_variants_with_vep_minimal(monkeypatch):
    monkeypatch.setattr(mod, "IMPACT_RANK", {"MODERATE": 3})
    monkeypatch.setattr(mod, "TERM_TO_FLAG", {"missense_variant":"IsMissense"})
    monkeypatch.setattr(mod, "TYPE_FLAG_COLS", ["IsSNV","IsDeletion","IsDuplication","IsInsertion","IsIndel"])

    df = pd.DataFrame([{
        "Chromosome":"17","Start":430,"End":430,"Ref":"C","Alt":"A","Type":"SNV","Origin":"CS"
    }])

    monkeypatch.setattr(mod, "fetch_base", lambda chrom, pos, cache: "G")
    monkeypatch.setattr(mod, "vep_call", lambda v: [{
        "most_severe_consequence":"missense_variant",
        "transcript_consequences":[{"biotype":"protein_coding","impact":"MODERATE","protein_start":1,"protein_end":1,"consequence_terms":["missense_variant"]}]
    }])

    out = mod.enrich_patient_variants_with_vep(df)
    assert "ConsequenceTerms" in out.columns
    assert "impact" in out.columns
    assert "IsCoding" in out.columns
    assert out.loc[0,"IsCoding"] in (0,1)
    assert out["IsSNV"].iloc[0] == 1
