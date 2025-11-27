import streamlit as st
from _utils import sidebar
import base64  # certifique-se de ter isso no topo do arquivo

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
O câncer de mama permanece um desafio de alta incidência e complexidade biológica. No contexto hereditário, alterações em BRCA1 e BRCA2 estão associadas a risco aumentado e motivam a classificação de patogenicidade de variantes nesses genes, com impacto direto em aconselhamento genético e tomada de decisão clínica.
</p>

<p class="justify">
Este trabalho propõe um pipeline de IA para classificar a patogenicidade de variantes em BRCA1 e BRCA2 a partir de um sequenciamento genético (fasta) do paciente. As variantes são detectadas por alinhamento à referência, anotadas via Ensembl VEP e transformadas em features binárias e numéricas. O classificador é um MLP multitarefa (uma cabeça por gene) com classificação ordinal (CORAL), que modela quatro classes ordenadas por meio de logits cumulativos, respeitando a estrutura ordinal do problema. O sistema é calibrado para alto recall na classe patogênica, priorizando a redução de falsos negativos clinicamente críticos.
</p>

<p class="justify">
A acurácia estrutural oferecida pela AlphaFold ampliou o acesso a modelos tridimensionais, incluindo BRCA1 e BRCA2, permitindo contextualizar variantes em regiões e domínios da proteína. Neste trabalho, essa visualização é utilizada como apoio à interpretação, enquanto a predição de patogenicidade permanece a tarefa principal.
</p>

<p class="justify">
Para treinamento e validação, utilizamos dados públicos; para testes ponta a ponta, o FASTA do “paciente” é sintético, derivado da sequência de referência com mutações in silico que emulam casos reais, permitindo verificar o pipeline sem expor dados clínicos.
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
    st.markdown("### 🔮 Predição")
    st.write("Envio de arquivos para predição de patogenicideade. (Em breve)")
    st.page_link("pages/predict.py", label="Ir para tela de Predição", icon="➡️")
with c3:
    st.markdown("### 📑 Resultados")
    st.write("Resultados e métricas do modelo")
    st.page_link("pages/results.py", label="Ir para tela de Resultados", icon="➡️")

st.divider()

st.markdown("<br>", unsafe_allow_html=True)

col_repo, col_doc = st.columns(2)
with col_repo:
    st.markdown("### 💻 Repositório no GitHub")
    st.write(
        "Acesse o código-fonte completo do Aitheron AI, "
        "incluindo scripts de treino e pipelines."
    )

    st.link_button(
        label="Ir para o repositório",
        url="https://github.com/Aitheron/Aitheron-AI_Engine",
        type="secondary",
        icon=":material/open_in_new:",
    )

with col_doc:
    st.markdown("### 📖 RFC - Aitheron")
    st.write(
        "Consulte o documento técnico (RFC) com a descrição do "
        "pipeline, arquitetura do modelo e decisões de projeto."
    )

    rfc_path = "../files/rfc_aitheron.pdf"

    with open(rfc_path, "rb") as f:
        pdf_bytes = f.read()

    btn_col1, btn_col2 = st.columns(2)

    with btn_col1:
        with st.popover("Pré-visualizar RFC", icon=":material/description:"):
            base64_pdf = base64.b64encode(pdf_bytes).decode("utf-8")
            pdf_html = f"""
            <iframe src="data:application/pdf;base64,{base64_pdf}"
                    width="100%" height="600"
                    type="application/pdf"></iframe>
            """
            st.markdown(pdf_html, unsafe_allow_html=True)

    with btn_col2:
        st.download_button(
            label="Baixar RFC em PDF",
            data=pdf_bytes,
            file_name="aitheron_rfc.pdf",
            mime="application/pdf",
            type="secondary",
            icon=":material/download:",
        )