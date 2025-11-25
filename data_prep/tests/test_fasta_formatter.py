import runpy
from pathlib import Path

def test_fasta_formatter_writes_fasta(tmp_path, monkeypatch):
    ptxt = tmp_path / "patient.txt"
    ptxt.write_text("acgtNN\nacgt-xxx\n")
    monkeypatch.chdir(tmp_path)

    runpy.run_module("data_prep.fasta_formatter", run_name="__main__")

    out = tmp_path / "patient.fasta"
    assert out.exists()
    lines = out.read_text().strip().splitlines()
    assert lines[0].startswith(">sample_")
    seq = "".join(lines[1:])
    assert set(seq).issubset(set("ACGTURYKMSWBDHVN"))
    assert len(seq) == 10
