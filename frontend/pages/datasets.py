import os
import pandas as pd
import streamlit as st

from utils import sidebar
from client.api_client import APIClient, GeneSelection, resolve_api_base_url

from analysis.plot import (
    preprocess,
    plot_significance_counts,
    plot_hotspots_by_protein,
    plot_top_alleles,
    plot_consequence_by_significance,
    plot_impact_by_significance,
    plot_protein_isoform_mix,
    plot_type_by_significance,
    plot_coding_vs_noncoding_by_gene_and_sig,
    plot_gene_hotspots,
)

st.set_page_config(page_title="BRCA1/2 – Datasets", layout="wide")

DATASET_PATH = '../files/clinvar_BRCA1_BRCA2_GRCh38_merged.xlsx'

st.title("📊 Datasets")
st.caption("Baixe os datasets utilizados (com cache no servidor) e visualize análises padronizadas.")
API_BASE_URL = resolve_api_base_url()
sidebar()

st.subheader("Opções de download")

left, right = st.columns([1, 2])
with left:
    gene_option = st.selectbox(
        "Selecione o dataset",
        options=[GeneSelection.BRCA1, GeneSelection.BRCA2, GeneSelection.BOTH],
        index=0,
        format_func=lambda x: {"BRCA1": "BRCA1", "BRCA2": "BRCA2", "BOTH": "BRCA1 + BRCA2"}[x.value],
    )
with right:
    st.info(
        "Como funciona:\n\n"
        "- Se o dataset já estiver em cache no servidor, ele é devolvido imediatamente.\n"
        "- Caso contrário, ele é gerado e armazenado para futuras requisições.\n"
        "- Com **Forçar reconstrução**, você reprocessa antes do download."
    )

if "download_bytes" not in st.session_state:
    st.session_state.download_bytes = None
    st.session_state.download_name = None

cta_col, dl_col = st.columns([1, 2])
with cta_col:
    if st.button("⬇️ Preparar download"):
        with st.spinner("Preparando o arquivo..."):
            client = APIClient(API_BASE_URL)
            try:
                data = client.download_dataset_bytes(gene_option, force=False)
                filename = {
                    GeneSelection.BRCA1: "dataset_BRCA1.xlsx",
                    GeneSelection.BRCA2: "dataset_BRCA2.xlsx",
                    GeneSelection.BOTH:  "dataset_BRCA1_BRCA2.xlsx",
                }[gene_option]
                st.session_state.download_bytes = data
                st.session_state.download_name = filename
                st.success("Dataset pronto para download!")
            except Exception as e:
                st.error(f"Falha ao obter o dataset: {e}")

with dl_col:
    if st.session_state.download_bytes:
        st.download_button(
            label="📥 Clique para baixar",
            data=st.session_state.download_bytes,
            file_name=st.session_state.download_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            use_container_width=True
        )

@st.cache_data(show_spinner=False)
def _read_fixed_xlsx(path: str) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    frames = [xls.parse(s, dtype=str) for s in xls.sheet_names]
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

st.divider()
st.subheader("📈 Interpretação da Base de Dados")

if not os.path.exists(DATASET_PATH):
    st.error(
        "Arquivo de dataset não encontrado.\n\n"
        f"Tente definir a variável de ambiente DATASET_PATH ou copiar o arquivo para:\n\n{os.path.abspath(DATASET_PATH)}"
    )
    st.stop()

with st.spinner("Lendo a planilha fixa..."):
    df_raw = _read_fixed_xlsx(DATASET_PATH)

if df_raw.empty:
    st.warning("A planilha foi lida, mas está vazia.")
    st.stop()

df = preprocess(df_raw)

st.subheader("Filtros")
col1, col2, col3 = st.columns(3)
with col1:
    genes = sorted(df["GeneSymbol"].dropna().unique().tolist())
    gene_sel = st.multiselect("Filtrar por GeneSymbol", options=genes, default=genes)
with col2:
    sigs = sorted(df["ClinicalSignificance"].dropna().unique().tolist())
    sig_sel = st.multiselect("Filtrar por Clinical Significance", options=sigs, default=sigs)
with col3:
    rs_text = st.text_input("Filtrar por RS# (dbSNP) contém...", value="")

mask = pd.Series(True, index=df.index)
mask &= df["GeneSymbol"].isin(gene_sel)
mask &= df["ClinicalSignificance"].isin(sig_sel)
if rs_text:
    mask &= df["RS# (dbSNP)"].fillna("").str.contains(rs_text, case=False, regex=False)

dff = df[mask].copy()
st.caption(f"Linhas após filtro: {len(dff):,}")

st.subheader("1) Distribuição geral")
plot_significance_counts(st, dff)

st.subheader("2) Clinical Significance × Tipo de mutação")
plot_type_by_significance(st, dff)

st.subheader("3) Codificante vs. não-codificante por gene (com Clinical Significance)")
plot_coding_vs_noncoding_by_gene_and_sig(st, dff)

st.subheader("4) Mix de isoformas por gene (somente codificantes)")
plot_protein_isoform_mix(st, dff)

st.subheader("5) Hotspots por posição na proteína (somente codificantes)")
bin_size = st.sidebar.number_input("Tamanho do bin (aa) para hotspots de proteína", min_value=10, max_value=300, value=50, step=10)
plot_hotspots_by_protein(st, dff, bin_size=bin_size)

st.subheader("6) Top pares de alelos (todas as variantes)")
top_n_alleles = st.slider("Top N de pares", 5, 40, 20, 1)
plot_top_alleles(st, dff, top_n=top_n_alleles)

st.subheader("7) Consequence terms × Clinical Significance (somente codificantes)")
top_n_terms = st.slider("Top N de consequence terms", 5, 40, 20, 1)
plot_consequence_by_significance(st, dff, top_n=top_n_terms)

st.subheader("8) Impact × Clinical Significance (somente codificantes)")
plot_impact_by_significance(st, dff)

st.subheader("9) Hotspots por posição no gene (todas as variantes)")
bin_bp = st.sidebar.number_input("Tamanho do bin (bp) para hotspots no gene", min_value=50, max_value=5000, value=100, step=50)
plot_gene_hotspots(st, dff, bin_bp=int(bin_bp))

with st.expander("Ver amostra da tabela filtrada (primeiras 200 linhas)"):
    st.dataframe(dff.head(200))
