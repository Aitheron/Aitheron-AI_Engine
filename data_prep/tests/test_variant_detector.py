from textwrap import dedent
import data_prep.variant_detector as mod

def test_read_fasta_string_reads_records(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(dedent(""">id1
        ACGT
        >id2
        TTAA
    """).strip())
    recs = mod.read_fasta_string(str(f))
    ids = [r[0] for r in recs]
    seqs = [r[1] for r in recs]
    assert ids == ["id1","id2"]
    assert seqs == ["ACGT","TTAA"]

from textwrap import dedent
import data_prep.variant_detector as mod

def test_read_fasta_string_reads_records(tmp_path):
    f = tmp_path / "x.fasta"
    f.write_text(">id1\nACGT\n>id2\nTTAA\n")
    recs = mod.read_fasta_string(str(f))
    ids = [r[0] for r in recs]
    seqs = [r[1] for r in recs]
    assert ids == ["id1", "id2"]
    assert seqs == ["ACGT", "TTAA"]

def test_coalesce_microindels_and_parse_cs():
    events = [
        {"Chromosome":"17","Start":10,"End":10,"Ref":"-","Alt":"A","Type":"Inserção","Origin":"CS"},
        {"Chromosome":"17","Start":11,"End":12,"Ref":"TT","Alt":"-","Type":"Deleção","Origin":"CS"},
    ]
    merged = mod.coalesce_microindels(events)
    assert any(e["Type"] == "Indel" for e in merged)

    # casos simples de cs: SNV, deleção e inserção
    ref_seq = "TTTACAGG"
    out_snv = mod.parse_cs_to_variants("*AC", ref_start_1based=1, ref_name="17", ref_seq=ref_seq, locus_start_1b=1)
    out_del = mod.parse_cs_to_variants("-TA", ref_start_1based=1, ref_name="17", ref_seq=ref_seq, locus_start_1b=1)
    out_ins = mod.parse_cs_to_variants("+GG", ref_start_1based=1, ref_name="17", ref_seq=ref_seq, locus_start_1b=1)

    types = {v["Type"] for v in (out_snv + out_del + out_ins)}
    # garante que pelo menos um de cada tipo foi produzido
    assert "SNV" in types
    assert "Deleção" in types
    assert "Inserção" in types

