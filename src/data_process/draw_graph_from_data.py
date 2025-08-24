import numpy as np
import matplotlib.pyplot as plt

from data_process.utils.iterables import read_data

def draw_data_graph_in_bars_style(
    *,
    data: dict,
    ref_column : str,
    graph_title : str,
    x_label : str,
    y_label : str,
    output_path: str
):
    # 1) pega DataFrames
    df1 = data["BRCA1"]
    df2 = data["BRCA2"]

    # 2) checagens mínimas
    if ref_column not in df1.columns or ref_column not in df2.columns:
        raise KeyError(f"Both DataFrames must contain a {ref_column} column.")
    
    # 3) conta por tipo
    c1 = df1[ref_column].value_counts()
    c2 = df2[ref_column].value_counts()

    # 4) categorias (união) e vetores alinhados
    categories = sorted(set(c1.index).union(set(c2.index)))
    idx = np.arange(len(categories))
    y1 = np.array([int(c1.get(cat, 0)) for cat in categories])
    y2 = np.array([int(c2.get(cat, 0)) for cat in categories])

    # 5) barras lado a lado (BRCA1 esquerda; BRCA2 direita)
    width = 0.4
    plt.figure(figsize=(10, 6))
    bars1 = plt.bar(idx - width/2, y1, width=width, label="BRCA1")
    bars2 = plt.bar(idx + width/2, y2, width=width, label="BRCA2")

    # 6) adicionar os valores em cima das barras
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            plt.text(bar.get_x() + bar.get_width()/2, height, str(height),
                     ha="center", va="bottom", fontsize=9)
    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            plt.text(bar.get_x() + bar.get_width()/2, height, str(height),
                     ha="center", va="bottom", fontsize=9)

    plt.xticks(idx, categories, rotation=0)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(graph_title)
    plt.legend()
    plt.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=600)
    plt.close()

    return output_path

if __name__ == "__main__":
    
    brca1_data = read_data('../files/clinvar_BRCA1_GRCh38_filtered.tsv')
    brca2_data = read_data('../files/clinvar_BRCA2_GRCh38_filtered.tsv')
    
    data_obj = {
        "BRCA1" : brca1_data,
        "BRCA2" : brca2_data
    }
    
    draw_data_graph_in_bars_style(
        data=data_obj,
        ref_column="ClinicalSignificance",
        graph_title='Distribuição das variantes por classificação clínica — BRCA1 e BRCA2',
        y_label='Quantidade de variantes',
        x_label='Classificação Clínica',
        output_path='../images/classifications_brca_bars.png'
    )
    
    draw_data_graph_in_bars_style(
        data=data_obj,
        ref_column="Type",
        graph_title="Tipo e quantidade de mutações por gene — BRCA1 e BRCA2",
        y_label="Quantidade de variantes",
        x_label="Tipo de mutação",
        output_path="../images/mutation_types_brca_bars.png"
    )
