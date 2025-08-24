# Resumo do Projeto

<div align="justify">

Este projeto prevê a aplicação de técnicas de Inteligência Artificial à análise de variantes nos genes BRCA1/BRCA2 no contexto do câncer de mama. Serão desenvolvidos dois modelos de Inteligência Artificial (IA): (i) um classificador supervisionado para rotular variantes (patogênica, benigna ou de significado incerto) e (ii) um modelo probabilístico que estima o potencial patogênico de cada variante. A base de dados será consolidada a partir de ClinVar e dbSNP, normalizada para GRCh38 e enriquecida por anotações obtidas via APIs do Ensembl e por informações estruturais do UniProt, com uso de recursos como AlphaFold quando pertinente. As representações combinam atributos tabulares (coordenadas, tipo e consequência, metadados de curadoria) e informações derivadas de sequência (janelas de DNA/proteína). O segundo modelo enfatiza o impacto funcional considerando posição e transcrito afetados, tipo de evento (missense, nonsense, frameshift, alterações de splicing) e efeito esperado na proteína; eventos truncantes e alterações com perda de função tendem a elevar a probabilidade prevista. O desempenho será avaliado com métricas apropriadas e calibração de probabilidades, visando priorização interpretável de variantes e um fluxo reprodutível útil à pesquisa translacional em oncologia de precisão.

</div>

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

---

# Perguntas para Especialistas

## Pergunta 1

Quando o ClinVar ou o Ensembl fornece uma anotação em **HGVS**, posso assumir que os alelos listados (ref/alt) já estão alinhados à sequência oficial (**GRCh38**)? E, de forma mais ampla, como garantir o uso da sequência correta (*wild-type*) para os genes **BRCA1/BRCA2**, considerando fontes distintas como **L78833.1** e **GRCh38**, bem como validar ou migrar coordenadas entre diferentes versões do genoma humano (ex.: GRCh37 vs GRCh38)?

---

## Pergunta 2
Quais são as recomendações atuais para padronizar a representação de variantes genéticas em relatórios clínicos e bases computacionais, especialmente no caso de mutações complexas como indels ou rearranjos maiores? Existem convenções internacionais (ex.: HGVS, VCF, guidelines clínicas) consideradas referência obrigatória em contextos médicos e de pesquisa?

---

## Pergunta 3
Em quais cenários variantes do tipo *missense* ou grandes rearranjos genômicos (LGRs) exigem estudos funcionais ou estruturais complementares? Quais são os critérios clínicos mínimos atualmente aceitos para classificá-las como patogênicas?

---

## Pergunta 4
Na prática clínica, quais seriam os erros mais graves e frequentes na interpretação de variantes em **BRCA1/BRCA2** que devo evitar ao conduzir meu projeto?

---

## Pergunta 5
Do ponto de vista clínico, ético e legal, quais são as principais limitações e responsabilidades ao propor um modelo de **Inteligência Artificial** que classifica variantes genéticas (patogênicas, benignas, VUS, etc.) e prediz probabilidades de risco?  
Além da **LGPD**, existe alguma legislação, regulamentação ou normativa específica que trate do uso de dados genômicos e sensíveis em pesquisas médicas no Brasil?
