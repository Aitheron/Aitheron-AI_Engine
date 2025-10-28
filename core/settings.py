from pathlib import Path

# Mapeia cada consequence_term (VEP/SO) -> coluna binária Is*
TERM_TO_FLAG = {
    # coding / protein-level
    "frameshift_variant":                "IsFrameshiftVariant",
    "missense_variant":                  "IsMissenseVariant",
    "synonymous_variant":                "IsSynonymousVariant",
    "inframe_insertion":                 "IsInframeInsertion",
    "inframe_deletion":                  "IsInframeDeletion",
    "protein_altering_variant":          "IsProteinAlteringVariant",
    "coding_sequence_variant":           "IsCodingSequenceVariant",
    "stop_gained":                       "IsStopGained",
    "stop_lost":                         "IsStopLost",
    "stop_retained_variant":             "IsStopRetainedVariant",
    "start_lost":                        "IsStartLost",

    # splices
    "splice_acceptor_variant":           "IsSpliceAcceptorVariant",
    "splice_donor_variant":              "IsSpliceDonorVariant",
    "splice_donor_5th_base_variant":     "IsSpliceDonor5thBaseVariant",
    "splice_donor_region_variant":       "IsSpliceDonorRegionVariant",
    "splice_polypyrimidine_tract_variant":"IsSplicePolypyrimidineTractVariant",
    "splice_region_variant":             "IsSpliceRegionVariant",

    # não-codificantes
    "intron_variant":                    "IsIntronVariant",
    "downstream_gene_variant":           "IsDownstreamGeneVariant",
    "upstream_gene_variant":             "IsUpstreamGeneVariant",
    "3_prime_UTR_variant":               "Is3primeUTRVariant",
    "5_prime_UTR_variant":               "Is5primeUTRVariant",
}

IMPACT_RANK = {"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}

CLINVAR_TAB_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

OUTDIR = Path(__file__).resolve().parents[1] / "filess"
OUTDIR.mkdir(exist_ok=True)
VALID_CONFIDENCE_NCBI = ["reviewed by expert panel", "practice guideline"]

TYPE_FLAG_COLS = [
    "IsSNV", "IsDeletion", "IsDuplication", "IsInsertion", "IsIndel"
]