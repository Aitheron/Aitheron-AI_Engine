import streamlit as st
from utils import sidebar

st.set_page_config(page_title="Análise BRCA", page_icon="🧬", layout="wide")

sidebar()

# Hero
st.markdown("# 🧬 Análise de Variantes BRCA")
st.write("Navegue pelo menu ou escolha abaixo o que deseja fazer:")

# “Cards” clicáveis (links oficiais do Streamlit)
c1, c2, c3 = st.columns(3)
with c1:
    st.markdown("### 📊 Datasets")
    st.write("Baixe os datasets do projeto (BRCA1, BRCA2 ou ambos), com cache no servidor.")
    st.page_link("pages/datasets.py", label="Ir para Datasets", icon="➡️")
with c2:
    st.markdown("### 🧠 Treino")
    st.write("Configuração e início do treinamento do modelo. (Em breve)")
    st.page_link("pages/train_model.py", label="Ver tela de Treino", icon="➡️")
with c3:
    st.markdown("### 🔮 Predição")
    st.write("Envio de arquivos para predição de patogenicideade. (Em breve)")
    st.page_link("pages/predict.py", label="Ver tela de Predição", icon="➡️")

st.divider()