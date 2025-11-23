import pandas as pd

from pages import predict

def test_classify_from_label_basic():
    # verifica mapeamentos diretos
    assert predict.classify_from_label({predict.CLASS_COL: "Patogênico"}) == "Patogênico"
    assert predict.classify_from_label({predict.CLASS_COL: "Benigno"}) == "Benigno"
    assert predict.classify_from_label({predict.CLASS_COL: "VUS"}) == "VUS"
    assert predict.classify_from_label({predict.CLASS_COL: "Possivelmente Benigno"}) == "Provavelmente benigno"


def test_prepare_dataframe_splits_coding_and_non():
    df = pd.DataFrame({
        predict.CLASS_COL: ["Patogênico", "Benigno", "VUS"],
        predict.CODING_FLAG_COL: ["Sim", "Não", "Sim"],
        predict.PROT_START: [10, None, 50],
        predict.PROT_END:   [12, None, 55],
        predict.CONSEQUENCE_COL: ["missense_variant", "intergenic_variant", "frameshift_variant"],
    })
    df_all, df_coding, df_non = predict.prepare_dataframe(df)
    # 2 codificantes, 1 não
    assert len(df_coding) == 2
    assert len(df_non) == 1
    # colunas auxiliares foram criadas
    assert "_ClasseRender" in df_all.columns and "_CorRender" in df_all.columns

def test_build_records_and_expand_range_and_stagger():
    dfc = pd.DataFrame({
        predict.PROT_START: [5, 100],
        predict.PROT_END:   [7, 102],
        predict.CONSEQUENCE_COL: ["missense_variant", "frameshift_variant"],
        "_ClasseRender": ["Patogênico", "Benigno"],
        "_CorRender": ["red", "green"],
    })
    recs = predict.build_records(dfc)
    assert len(recs) == 2
    assert recs[0]["start"] == 5 and recs[0]["end"] == 7
    assert isinstance(predict.expand_range(10, 12, pad=2), list)
    assert predict.stagger_offset(0) != predict.stagger_offset(1)

def test_prepare_dataframe_without_coding_col_falls_back_to_protein_cols():
    df = pd.DataFrame({
        predict.CLASS_COL: ["Patogênico", "VUS"],
        predict.PROT_START: [1, None],
        predict.PROT_END:   [2, None],
        predict.CONSEQUENCE_COL: ["missense_variant", "upstream_gene_variant"],
    })
    _, df_coding, df_non = predict.prepare_dataframe(df)
    assert len(df_coding) == 1 and len(df_non) == 1
