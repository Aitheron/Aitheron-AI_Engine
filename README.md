# Aitheron — Classificação de Patogenicidade em BRCA1/BRCA2

<div align="justify">

O câncer de mama permanece um desafio de alta incidência e complexidade biológica.
No contexto hereditário, alterações em BRCA1 e BRCA2 estão associadas a risco aumen-
tado e motivam a classificação de patogenicidade de variantes nesses genes, com impacto
direto em aconselhamento genético e tomada de decisão clínica.

Este trabalho propõe um pipeline de IA para classificar a patogenicidade de variantes
em BRCA1 e BRCA2 a partir de um sequenciamento genético (fasta) do paciente. As
variantes são detectadas por alinhamento à referência, anotadas via Ensembl VEP e
transformadas em features binárias e numéricas. O classificador é um MLP multitarefa
(uma cabeça por gene) com classificação ordinal (CORAL), que modela quatro classes
ordenadas por meio de logits cumulativos, respeitando a estrutura ordinal do problema. O
sistema é calibrado para alto recall na classe patogênica, priorizando a redução de falsos
negativos clinicamente críticos.

A acurácia estrutural oferecida pela AlphaFold ampliou o acesso a modelos tridimen-
sionais, incluindo BRCA1 e BRCA2, permitindo contextualizar variantes em regiões e
domínios da proteína. Neste trabalho, essa visualização é utilizada como apoio à interpre-
tação, enquanto a predição de patogenicidade permanece a tarefa principal.

Para treinamento e validação, utilizamos dados públicos; para testes ponta a ponta, o
FASTA do “paciente” é sintético, derivado da sequência de referência com mutações in
silico que emulam casos reais, permitindo verificar o pipeline sem expor dados clínicos.

</div>


---

## 🔧 Requisitos

- Python 3.11+
- (Opcional) venv

## 🚀 Instalação (ambiente local)

Crie e ative um ambiente virtual e instale as dependências do projeto.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 🏋️ Treinar o modelo

O treino completo (ETL → treino → avaliação) é disparado pelo script **`run-training.sh`**.

> Caso veja “Permission denied”, torne o script executável:
```bash
chmod +x ./run-training.sh
```

Agora rode a partir da **raiz do repositório**:
```bash
./run-training.sh
```

Saídas esperadas (exemplo):
- `artifacts/.../model.pth` — checkpoint do modelo
- `artifacts/.../metrics.json` — métricas de validação
- `artifacts/.../preprocessor.pkl` — hiperparâmetros e configurações usadas

---

## 🌐 Subir o sistema (API + Web UI)

Use **`run-saver.sh`** para subir **FastAPI** (porta `8000`) e **Streamlit** (porta `8501`).

> Se necessário, dê permissão:
```bash
chmod +x ./run-saver.sh
```

Execute **sempre da raiz do projeto**:
```bash
./run-saver.sh
```

---

## 🧪 Testes e cobertura

Rodar testes:
```bash
pytest -q
```

Cobertura (com faltantes):
```bash
pytest --cov=. --cov-report=term-missing
```

---

## Amostragem dos dados

| VariationID | GeneSymbol | ClinicalSignificance | ClinSigSimple | Type                       | Origin           | Assembly | Chromosome |   Start   |   Stop    |
|-------------|------------|----------------------|---------------|----------------------------|-----------------|----------|------------|-----------|-----------|
| 209219      | BRCA1      | Benign               | 0             | single nucleotide variant  | germline        | GRCh38   | 17         | 43039471  | 43039471  |
| 55602       | BRCA1      | Pathogenic           | 1             | Deletion                   | germline        | GRCh38   | 17         | 43045706  | 43045767  |
| 209597      | BRCA2      | Benign               | 0             | single nucleotide variant  | germline        | GRCh38   | 13         | 32314943  | 32314943  |
| 51579       | BRCA2      | Pathogenic           | 1             | single nucleotide variant  | germline;unknown| GRCh38   | 13         | 32316463  | 32316463  |

---

## Distribuição das variantes por classificação clínica — BRCA1 e BRCA2

![Gráfico de Distribuição das variantes por classificação clinica](./images/classifications_brca_bars.png)

---

## Tipo e quantidade de mutações por gene - BRCA1 e BRCA2

![Tipo e quantidade de mutações por gene](./images/mutation_types_brca_bars.png)


## 🧩 Saída ordinal com CORAL (Cumulative Link)

Cada cabeça (BRCA1/BRCA2) produz \(K{-}1\) **logits cumulativos**; com \(K=4\), são **3 logits** que, após `sigmoid`, viram probabilidades cumulativas \(q\). As probabilidades por classe são:

```math
\begin{aligned}
P(y{=}0) &= 1 - q_{0} \\
P(y{=}1) &= q_{0} - q_{1} \\
P(y{=}2) &= q_{1} - q_{2} \\
P(y{=}3) &= q_{2}
\end{aligned}
```

A classe final é o `argmax` de \([P(y{=}0),P(y{=}1),P(y{=}2),P(y{=}3)]\). O modelo retorna - para cada predição - a classificação clinica (benigno, possívelmente benigno, vus ou patogênico), a confiança calculada através da entropia e a probabilidade de a mutação ser patogênica.

---

## 🧠 Como calculamos a Confiança (entropia)

A confiança é derivada da **entropia normalizada** dos logits gerados pelo modelo para cada classificação:

```math
\text{Confiança} \;=\; 1 \;-\; \frac{H(p)}{\ln K}
```

onde 
```math 
(p = [P(y{=}0),\,P(y{=}1),\,\dots,\,P(y{=}K{-}1)])
```
e:

```math
H(p) \;=\; - \sum_{k=0}^{K-1} p_k \,\ln p_k
```

Para \(K=4\):
```math
H_{\text{norm}}(p) = \frac{H(p)}{\ln 4}, \qquad \text{Confiança} = 1 - H_{\text{norm}}(p)
```

---

## 🗺️ Arquitetura (visão geral)

- **Frontend (Streamlit)**: upload de FASTA, filtros e visualização (incl. 3D AlphaFold).
- **Backend (FastAPI)**: endpoints para datasets e predição.
- **Serviços de dados**: ETL/anotação com VEP/ClinVar/UniProt/AlphaFold.
- **Modelo**: MLP ordinal multitarefa (CORAL), calibrado para alto recall em Patogênico.

### Decisões de modelagem — justificativas e fontes

| Decisão | Justificativa | Fonte |
|----------|----------------|--------|
| **Hotspots binários (Is\*)** | Evitam ordens artificiais, favorecem interpretabilidade e compatibilidade com ML tabular. | [scikit-learn OneHotEncoder](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OneHotEncoder.html) |
| **MLP multitarefa** | Simples, robusto, boa performance em dados mistos, fácil manutenção e comparação. | [Overview NAR 2023](https://academic.oup.com/nar/article/51/D1/D1095/6848490) |
| **Ordinal (CORAL)** | Mantém coerência entre fronteiras e usa ordenação natural. | [CORAL paper (Cao et al., 2019)](https://arxiv.org/abs/1901.07884) |
| **Cross-validation (10-fold)** | Reduz variância e evita overfitting em datasets pequenos. | [Kohavi, 1995](https://www.cs.cornell.edu/people/tj/publications/joachims_97a.pdf) |
| **VEP como anotador** | Ferramenta padrão de mercado para *consequence_terms* e impactos. | [Ensembl VEP](https://rest.ensembl.org) |
| **AlphaFold para visualização** | Modelos 3D precisos e amplamente aceitos na comunidade. | [AlphaFold DB](https://alphafold.ebi.ac.uk/) |
| **ACMG/AMP diretrizes** | Padrão clínico global para interpretação de variantes genéticas. | [ACMG/AMP 2015](https://pubmed.ncbi.nlm.nih.gov/25741868/) |
| **LGPD conformidade** | Obrigatório para tratamento de dados genéticos sensíveis. | [Lei 13.709/2018](https://www.planalto.gov.br/ccivil_03/_ato2015-2018/2018/lei/L13709.htm) |

---

## 🛠️ Dicas e solução de problemas

- **Permission denied ao rodar `.sh`**  
  Use `chmod +x ./run-training.sh` e/ou `chmod +x ./run-saver.sh`.
- **Portas ocupadas (8000/8501)**  
  Troque as portas nos scripts ou finalize processos que estejam usando essas portas.
- **Erros de import em testes**  
  O `pytest.ini` define `PYTHONPATH` apropriado. Evite rodar os testes dentro de subpastas; execute na raiz do repo.

---

## 📚 Referências principais

- **CORAL** (Cao et al., 2019): classificação ordinal por ligações cumulativas.  
  https://arxiv.org/abs/1901.07884
- **ACMG/AMP (2015)**: diretrizes para interpretação clínica de variantes.  
  https://pubmed.ncbi.nlm.nih.gov/25741868/
- **Ensembl VEP** (REST API e termos de consequência).  
  https://rest.ensembl.org  
  https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
- **AlphaFold DB** (estruturas 3D).  
  https://alphafold.ebi.ac.uk/
- **UniProt** (BRCA1/BRCA2).  
  https://www.uniprot.org/uniprotkb/P38398  
  https://www.uniprot.org/uniprotkb/P51587

---

## 🔐 Nota sobre privacidade (LGPD)

Tratamos dados genéticos como **dados pessoais sensíveis**. O pipeline de demonstração usa dados públicos e/ou FASTA sintético para testes ponta a ponta. Em cenários reais, siga princípios de **finalidade**, **minimização** e **segurança** (Lei 13.709/2018).


## Acesso a aplicação

> https://aitheron.com.br