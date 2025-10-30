import streamlit as st
from _utils import sidebar

st.set_page_config(page_title="Aitheron", page_icon="🧬", layout="wide")

sidebar()

st.markdown(
    "<h1 style='text-align:center'>🧬 Análise de Variantes Genéticas nos genes BRCA1 e BRCA2</h1>",
    unsafe_allow_html=True
)

st.markdown("""
<style>
.block-container { max-width: 65%; }
.justify { text-align: justify; text-justify: inter-word; hyphens: auto; line-height: 1.65; }
</style>
""", unsafe_allow_html=True)

st.header("Resumo do projeto")

st.markdown("""
<p class="justify">
Este projeto propõe a aplicação de técnicas de Inteligência Artificial (IA) à análise de variantes nos genes BRCA1/BRCA2 no contexto do câncer de mama. Serão desenvolvidos dois modelos de Inteligência Artificial (IA): (i) um classificador de variantes (patogênicas, benignas ou de significado incerto) e (ii) um modelo probabilístico capaz de estimar o potencial patogênico de cada variante.
</p>

<p class="justify">
A base de dados será consolidada a partir do ClinVar e dbSNP, normalizada para o genoma de referência GRCh38 e enriquecida com anotações obtidas via APIs do Ensembl, bem como com informações estruturais do UniProt e recursos do projeto AlphaFold, quando pertinente.
</p>

<p class="justify">
As representações das variantes combinarão atributos tabulares (coordenadas, tipo e consequência, metadados de curadoria) com informações derivadas de sequência (janelas de DNA/proteína). O segundo modelo dará ênfase ao impacto funcional, considerando posição e transcrito afetados, tipo de evento (missense, nonsense, frameshift, alterações de splicing) e efeito esperado na proteína. Alterações truncantes e mutações associadas à perda de função tendem a aumentar a probabilidade prevista.
</p>

<p class="justify">
O desempenho dos modelos será avaliado com métricas apropriadas, incluindo a calibração probabilística, de forma a possibilitar a priorização interpretável de variantes e oferecer um fluxo reprodutível de análise, útil à pesquisa translacional em oncologia de precisão.
</p>

<br>
<br>
""", unsafe_allow_html=True)


st.header("Navegue pelo menu ou escolha abaixo o que deseja fazer:")

st.markdown("<br>", unsafe_allow_html=True)


# “Cards” clicáveis (links oficiais do Streamlit)
c1, c2, c3 = st.columns(3)
with c1:
    st.markdown("### 📊 Datasets")
    st.write("Baixe os datasets do projeto (BRCA1, BRCA2 ou ambos), com cache no servidor.")
    st.page_link("pages/datasets.py", label="Ir para Datasets", icon="➡️")
with c2:
    st.markdown("### 🧠 Treino")
    st.write("Configuração e início do treinamento do modelo. (Em breve)")
    st.page_link("pages/train_model.py", label="Ir para tela de Treino", icon="➡️")
with c3:
    st.markdown("### 🔮 Predição")
    st.write("Envio de arquivos para predição de patogenicideade. (Em breve)")
    st.page_link("pages/predict.py", label="Ir para tela de Predição", icon="➡️")

st.divider()