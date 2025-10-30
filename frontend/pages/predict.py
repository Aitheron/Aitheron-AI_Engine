import streamlit as st
import pandas as pd
from io import BytesIO

from _utils import (
    sidebar,
    parse_single_fasta,
    validate_sequence,
    ALLOWED
)
from client.api_client import APIClient

sidebar()

def reset_state():
    for k in ("header","seq","gene","payload_preview","api_result","result_df","result_bytes","result_filename"):
        st.session_state.pop(k, None)

st.set_page_config(page_title="Predição", page_icon="🧬", layout="wide")
st.title("Upload de FASTA para Predição")

left, right = st.columns([7,3])

with right:
    st.info(
        """
        **Observações:**\n\n
        - Este fluxo **não faz autodetecção** de gene.\n
        - Envie **apenas um** registro FASTA por vez.\n
        - Alfabeto aceito: **A/C/G/T/N**.\n
        - Se enviar FASTA de BRCA1 com gene **BRCA2** selecionado (ou vice-versa), **não há correção automática**.
        """
    )
    st.divider()

with left:
    col_up, col_gene = st.columns([2,1])
    with col_up:
        uploaded = st.file_uploader("Selecione o arquivo FASTA", type=["fasta","fa","fna"], on_change=reset_state)
    with col_gene:
        gene = st.selectbox("Gene", ["BRCA1", "BRCA2"], key="gene_select")

    if uploaded is not None:
        file_bytes = uploaded.read()
        header, seq, err = parse_single_fasta(file_bytes)
        if err:
            st.error(err)
            st.stop()
        ok, msg = validate_sequence(seq)
        if not ok:
            st.error(msg)
            st.stop()

        st.session_state["header"] = header
        st.session_state["seq"] = seq
        st.session_state["gene"] = gene

        st.subheader("Resumo do FASTA")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.success("FASTA válido", icon="✅")
            st.write(f"**Header**: `{header}`")
        with c2:
            st.write(f"**Comprimento**: {len(seq)} bp")
            st.write(f"**Alfabeto**: {'/'.join(sorted(ALLOWED))}")
        with c3:
            st.write(f"**Gene selecionado**: **{gene}**")
            st.caption("Certifique-se de que o gene corresponde ao conteúdo do FASTA.")

        st.divider()

        file_preview = {
            "gene": st.session_state["gene"],
            "fasta_header": st.session_state["header"],
            "sequence": st.session_state["seq"],
        }
        st.session_state["payload_preview"] = file_preview

        st.expander("Prévia do arquivo", expanded=False).json(st.session_state["payload_preview"])
        btn = st.button("📤 Enviar para processamento", type="primary", use_container_width=True)
        if btn:
            client = APIClient()
            try:
                with st.spinner("Processando..."):
                    xlsx_bytes = client.predict_from_fasta_file(
                        gene=st.session_state["gene"],
                        file_name=uploaded.name,
                        file_bytes=file_bytes,
                    )
                filename = "predicted.xlsx"
                st.session_state["result_bytes"] = xlsx_bytes
                st.session_state["result_filename"] = filename or f"prediction_{gene}.xlsx"
                st.success("Processamento concluído.")
            except Exception as e:
                st.error(f"Falha na chamada de API: {e}")

        if "result_bytes" in st.session_state:
            st.subheader("Resultado da Predição")

            df = pd.read_excel(BytesIO(st.session_state["result_bytes"]), dtype=str)
            st.session_state["result_df"] = df
            with st.expander("Ver amostra da tabela filtrada (primeiras 200 linhas)"):
                st.dataframe(df, use_container_width=True, hide_index=True)

    else:
        st.info("Aguardando o upload de um arquivo FASTA.")
