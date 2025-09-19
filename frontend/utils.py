import streamlit as st

def sidebar():
    with st.sidebar:
        st.markdown("## 🧬 Análise BRCA")
        st.caption("Navegue pelas páginas do app:")
        st.page_link("app.py", label="Início", icon="🏠")
        st.page_link("pages/datasets.py", label="Datasets", icon="📊")
        st.page_link("pages/train_model.py", label="Treino (em breve)", icon="🧠")
        st.page_link("pages/predict.py", label="Predição (em breve)", icon="🔮")
        st.divider()