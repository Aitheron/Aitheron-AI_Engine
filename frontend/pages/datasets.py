import streamlit as st

from utils import sidebar

from client.api_client import APIClient, GeneSelection, resolve_api_base_url

st.title("📊 Datasets")
st.caption("Baixe os datasets utilizados (com cache no servidor).")

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
    force_rebuild = st.toggle(
        "Forçar reconstrução (reprocessar agora)",
        value=False,
        help="Regera o dataset no servidor antes de baixar."
    )
    

with right:
    st.info(
        "Como funciona:\n\n"
        "- Se o dataset já estiver em cache no servidor, ele é devolvido imediatamente.\n"
        "- Caso contrário, ele é gerado e armazenado para futuras requisições.\n"
        "- Com **Forçar reconstrução**, você reprocessa antes do download."
    )

st.divider()

if "download_bytes" not in st.session_state:
    st.session_state.download_bytes = None
    st.session_state.download_name = None

cta_col, dl_col = st.columns([1, 2])
with cta_col:
    if st.button("⬇️ Preparar download"):
        with st.spinner("Preparando o arquivo..."):
            client = APIClient(API_BASE_URL)
            try:
                data = client.download_dataset_bytes(gene_option, force=force_rebuild)
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
            mime="text/tab-separated-values",
            use_container_width=True
        )
