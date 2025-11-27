import streamlit as st
import pandas as pd

from _utils import (
    sidebar
)

st.set_page_config(page_title="Resultados", page_icon="🧬", layout="wide")
sidebar()

st.title("📑 Resultados dos Modelos BRCA1 e BRCA2")

METRICS = {
    "BRCA1": {
        "auroc": 0.9942128750394601,
        "auprc": 0.9970682720272839,
        "recall": 0.989805375347544,
        "precision": 0.9187096774193548,
        "f1": 0.9529333035913451,
        "balanced_acc": 0.9088371139032803,
        "mcc": 0.8555017073447976,
        "accuracy": 0.9351965601965602,
        "threshold": 0.07836013287305832,
    },
    "BRCA2": {
        "auroc": 0.9957554504044906,
        "auprc": 0.997459930988328,
        "recall": 0.9945504087193461,
        "precision": 0.9793024147182828,
        "f1": 0.9868675164156044,
        "balanced_acc": 0.9784729759474168,
        "mcc": 0.9630901319563844,
        "accuracy": 0.9830212234706617,
        "threshold": 0.15000000000000002,
    },
}

heatmap_path = "../images/confusion_heatmap.png"
pr_brca1_path = "../images/pr_BRCA1.png"
pr_brca2_path = "../images/pr_BRCA2.png"
roc_brca1_path = "../images/roc_BRCA1.png"
roc_brca2_path = "../images/roc_BRCA2.png"

col_left, col_right = st.columns([1, 1])

with col_left:
    st.subheader("Matrizes de confusão")
    st.image(heatmap_path, caption="Matrizes de confusão — BRCA1 e BRCA2", use_container_width=True)

with col_right:
    st.subheader("Métricas globais dos modelos")
    rows = []
    for gene, vals in METRICS.items():
        rows.append(
            {
                "Gene": gene,
                "AUROC": vals["auroc"],
                "AUPRC": vals["auprc"],
                "Acurácia": vals["accuracy"],
                "Recall (Sensibilidade)": vals["recall"],
                "Precisão": vals["precision"],
                "F1-score": vals["f1"],
                "Balanced Accuracy": vals["balanced_acc"],
                "MCC": vals["mcc"],
                "Threshold": vals["threshold"],
            }
        )
    df_metrics = pd.DataFrame(rows).set_index("Gene")
    st.dataframe(df_metrics.style.format("{:.4f}"))

with st.expander("O que cada métrica significa"):
    st.markdown(
    """
    - **AUROC (Area Under the ROC Curve)**  
    Mede a capacidade do modelo de separar classes ao variar o limiar de decisão.  
    Quanto mais perto de 1, melhor o poder discriminativo.

    - **AUPRC (Area Under the Precision–Recall Curve)**  
    Foca na relação entre **Precisão** e **Recall**, especialmente útil quando a classe positiva é rara.  
    Valores próximos de 1 indicam bom desempenho em identificar verdadeiros positivos sem muitos falsos positivos.

    - **Acurácia**  
    Proporção de predições corretas entre todas as amostras (TP + TN / total).  
    Pode ser enganosa se as classes forem muito desbalanceadas.

    - **Recall (Sensibilidade)**  
    Entre todos os casos realmente positivos, quantos o modelo identificou como positivos (TP / TP + FN).  
    Importante para não deixar passar variantes patogênicas.

    - **Precisão**  
    Entre todos os casos que o modelo chamou de positivos, quantos realmente eram positivos (TP / TP + FP).  
    Alta precisão significa poucos falsos positivos.

    - **F1-score**  
    Média harmônica entre Precisão e Recall.  
    Resume o trade-off entre os dois em um único número.

    - **Balanced Accuracy**  
    Média entre sensibilidade da classe positiva e da classe negativa.  
    Corrige o efeito de desbalanceamento entre classes.

    - **MCC (Matthews Correlation Coefficient)**  
    Métrica de correlação entre predições e rótulos verdadeiros, variando de -1 a 1.  
    Próximo de 1 indica forte concordância; próximo de 0, aleatório.

    - **Threshold**  
    Limiar de probabilidade usado para decidir se uma variante é considerada patogênica.  
    Probabilidades acima do threshold são classificadas como positivas.
    """
    )

st.markdown("---")

st.subheader("Curvas Precision–Recall (Patogênico vs Não-Patogênico)")
st.markdown(
    """
Essas curvas mostram o comportamento do modelo ao variar o limiar de decisão, 
relacionando **Recall** (quantos casos patogênicos o modelo encontra) com **Precisão** 
(quantos dos casos que ele chama de patogênicos realmente são patogênicos). 
Curvas mais próximas do canto superior direito indicam melhor desempenho.

O ponto vermelho em cada curva indica o limiar de decisão utilizado na prática para aquele gene.
Nele, o modelo opera com o par (Recall, Precisão) mostrado na legenda, representando o trade-off escolhido
entre não perder casos patogênicos e evitar falsos positivos em excesso.
"""
)

pr_col1, pr_col2 = st.columns(2)

with pr_col1:
    st.image(pr_brca1_path, caption="Precision–Recall — BRCA1", use_container_width=True)

with pr_col2:
    st.image(pr_brca2_path, caption="Precision–Recall — BRCA2", use_container_width=True)

st.markdown("---")

st.subheader("Curvas ROC (Patogênico vs Não-Patogênico)")
st.markdown(
    """
As curvas ROC mostram a relação entre **Taxa de Verdadeiros Positivos (Recall)** 
e **Taxa de Falsos Positivos** para todos os limiares possíveis. 
Curvas próximas do canto superior esquerdo indicam que o modelo consegue 
identificar a maior parte dos casos patogênicos mantendo poucos falsos positivos.

O ponto vermelho marca o ponto de operação atual do modelo, correspondente ao limiar de probabilidade
escolhido para cada gene (BRCA1 e BRCA2). Esse ponto mostra qual combinação de
sensibilidade (TPR) e taxa de falsos positivos (FPR) foi adotada no sistema.
"""
)

roc_col1, roc_col2 = st.columns(2)

with roc_col1:
    st.image(roc_brca1_path, caption="ROC — BRCA1", use_container_width=True)

with roc_col2:
    st.image(roc_brca2_path, caption="ROC — BRCA2", use_container_width=True)
