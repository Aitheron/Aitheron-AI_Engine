from pathlib import Path

CLINVAR_TAB_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

OUTDIR = Path(__file__).resolve().parents[2] / "files"
OUTDIR.mkdir(exist_ok=True)
VALID_CONFIDENCE_NCBI = ["reviewed by expert panel", "practice guideline"]