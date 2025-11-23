import pandas as pd
from pathlib import Path
import data_prep.clean_data as mod

def test_clean_data_filters_and_writes_file(tmp_path):
    src = tmp_path / "in.tsv"
    df = pd.DataFrame({
        "ReviewStatus": ["reviewed by expert panel", "no assertion"],
        "x": [1, 2]
    })
    df.to_csv(src, sep="\t", index=False)

    out_files = mod.clean_data_per_submit_confidence(
        input_files=[str(src)],
        valid_confidence=["reviewed by expert panel"]
    )
    assert len(out_files) == 1
    out = Path(out_files[0])
    assert out.exists()
    got = pd.read_csv(out, sep="\t", dtype=str)
    assert len(got) == 1
    assert list(got["x"]) == ["1"]

def test_clean_data_raises_when_no_confidence_col(tmp_path):
    src = tmp_path / "in.tsv"
    pd.DataFrame({"foo": [1]}).to_csv(src, sep="\t", index=False)
    try:
        mod.clean_data_per_submit_confidence([str(src)], ["reviewed by expert panel"])
        assert False, "should have raised"
    except ValueError as e:
        assert "No confidence/review column" in str(e)
