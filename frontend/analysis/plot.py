import pandas as pd
import altair as alt

EXPECTED_COLS = [
    "VariationID","GeneSymbol","ClinicalSignificance","ClinSigSimple","Type",
    "Origin","OriginSimple","Assembly","Chromosome","Start","Stop",
    "ReferenceAllele","AlternateAllele","ReviewStatus","RCVaccession","RS# (dbSNP)",
    "IsCoding","TranscriptSelectionBasis","ENST","ENSP","ConsequenceTerms","Impact",
    "VariantAllele","ProteinStart","ProteinEnd","CDNAStart","CDNAEnd","CDSStart","CDSEnd",
    "Codons","AminoAcids","Strand","HGVSc","HGVSp","FASTA_CDS","FASTA_Protein",
    "ProteinName","ProteinFunction","ProteinFeatures","ProteinDesc",
]

def _coerce_numeric(series: pd.Series) -> pd.Series:
    try:
        return pd.to_numeric(series, errors="coerce")
    except Exception:
        return pd.Series(dtype="float64")

def preprocess(df: pd.DataFrame) -> pd.DataFrame:

    if df.empty:
        return df
    out = df.copy()
    out.columns = out.columns.str.strip()

    for c in EXPECTED_COLS:
        if c not in out.columns:
            out[c] = pd.NA

    # strings
    for c in ["ClinicalSignificance","GeneSymbol","ReferenceAllele","AlternateAllele","ProteinName","Type"]:
        out[c] = out[c].astype(str).str.strip()

    # numéricos
    for c in ["ProteinStart","ProteinEnd","CDNAStart","CDNAEnd","CDSStart","CDSEnd","Start","Stop","Strand"]:
        out[c] = _coerce_numeric(out[c])

    # IsCoding para bool “limpo”
    out["IsCoding"] = out["IsCoding"].astype(str).str.lower().isin(["true"])

    # consequence -> lista
    out["ConsequenceTerms_raw"] = out["ConsequenceTerms"].fillna("").astype(str)
    out["ConsequenceTerms_list"] = out["ConsequenceTerms_raw"].apply(lambda x: [t for t in x.split("|") if t])

    return out

def explode_consequence_terms(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        cs = r.get("ClinicalSignificance", "NA")
        terms = r.get("ConsequenceTerms_list") or []
        if not terms:
            rows.append({"ConsequenceTerm": "NA", "ClinicalSignificance": cs})
        else:
            for t in terms:
                rows.append({"ConsequenceTerm": str(t), "ClinicalSignificance": cs})
    out = pd.DataFrame(rows)
    if "ConsequenceTerm" not in out.columns and "ConsequenceTerms" in out.columns:
        out = out.rename(columns={"ConsequenceTerms": "ConsequenceTerm"})
    if "ConsequenceTerm" not in out.columns:
        out["ConsequenceTerm"] = pd.Series(dtype="string")
    out["ConsequenceTerm"] = out["ConsequenceTerm"].astype("string")
    out["ClinicalSignificance"] = out["ClinicalSignificance"].astype("string")
    return out

def plot_significance_counts(st, df: pd.DataFrame):
    counts = df["ClinicalSignificance"].fillna("NA").value_counts().reset_index()
    counts.columns = ["ClinicalSignificance", "count"]
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("ClinicalSignificance:N", sort="-y", title="Clinical Significance"),
        y=alt.Y("count:Q", title="Contagem"),
        tooltip=["ClinicalSignificance","count"]
    ).properties(title="Distribuição de Clinical Significance", height=300)
    st.altair_chart(chart, use_container_width=True)

def plot_type_by_significance(st, df: pd.DataFrame):
    d = df.copy()
    d["Type"] = d["Type"].fillna("NA")
    agg = d.groupby(["Type","ClinicalSignificance"]).size().reset_index(name="count")
    chart = alt.Chart(agg).mark_bar().encode(
        x=alt.X("Type:N", sort="-y", title="Tipo de mutação"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        tooltip=["Type","ClinicalSignificance","count"]
    ).properties(title="Clinical Significance × Tipo de mutação", height=320).interactive()
    st.altair_chart(chart, use_container_width=True)

def plot_coding_vs_noncoding_by_gene_and_sig(st, df: pd.DataFrame):
    d = df.copy()
    d["CodingClass"] = d["IsCoding"].map({True:"Coding", False:"Non-coding"})
    agg = d.groupby(["GeneSymbol","CodingClass","ClinicalSignificance"]).size().reset_index(name="count")
    chart = alt.Chart(agg).mark_bar().encode(
        x=alt.X("GeneSymbol:N", title="Gene"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        column=alt.Column("CodingClass:N", title=None),
        tooltip=["GeneSymbol","CodingClass","ClinicalSignificance","count"]
    ).properties(title="Codificante vs. Não-codificante por gene (com Clinical Significance)", height=320)
    st.altair_chart(chart, use_container_width=True)


def _only_coding(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["IsCoding"] == True].copy()

def plot_hotspots_by_protein(st, df: pd.DataFrame, bin_size: int):
    # filtra apenas variantes codificantes com ProteinStart/ProteinName válidos
    data = df.copy()
    data = data[data["IsCoding"] == True]
    data = data.dropna(subset=["ProteinStart", "ProteinName", "GeneSymbol"])
    if data.empty:
        st.info("Sem variantes codificantes com ProteinStart/ProteinName.")
        return

    charts = []
    for gene in ["BRCA1", "BRCA2"]:
        g = data[data["GeneSymbol"] == gene].copy()
        if g.empty:
            continue

        # proteína mais prevalente dentro do gene
        top_proteins = g["ProteinName"].value_counts()
        if top_proteins.empty:
            continue
        top_protein = top_proteins.idxmax()

        g = g[g["ProteinName"] == top_protein].copy()
        if g.empty:
            continue

        # binning por posição na proteína
        g["bin"] = (g["ProteinStart"] // bin_size) * bin_size
        agg = g.groupby(["bin"]).size().reset_index(name="count")

        title = f"{gene} – hotspots (codificantes) – proteína mais prevalente: {top_protein} (bin={bin_size} aa)"
        chart = (
            alt.Chart(agg)
            .mark_bar()
            .encode(
                x=alt.X("bin:Q", title="Posição inicial do bin (aa)"),
                y=alt.Y("count:Q", title="Contagem"),
                tooltip=["bin","count"]
            )
            .properties(title=title, height=320)
        )
        charts.append(chart)

    if not charts:
        st.info("Não há dados codificantes suficientes para BRCA1/BRCA2.")
        return

    # empilha BRCA1 e BRCA2 (quando ambos existirem)
    combined = charts[0] if len(charts) == 1 else alt.vconcat(*charts)
    st.altair_chart(combined, use_container_width=True)

def plot_protein_isoform_mix(st, df: pd.DataFrame):
    data = _only_coding(df)
    if data.empty:
        st.info("Sem variantes codificantes para compor o mix de isoformas.")
        return
    agg = data.groupby(["GeneSymbol","ProteinName"]).size().reset_index(name="count")
    # pizza por gene
    chart = alt.Chart(agg).mark_arc(innerRadius=60).encode(
        theta="count:Q",
        color=alt.Color("ProteinName:N", title="Proteína"),
        tooltip=["GeneSymbol","ProteinName","count"]
    ).properties(title="Mix de isoformas (somente codificantes)", height=280)
    chart = chart.facet(column=alt.Column("GeneSymbol:N", title="Gene"))
    st.altair_chart(chart, use_container_width=True)

def plot_consequence_by_significance(st, df: pd.DataFrame, top_n: int = 20):
    d = _only_coding(df)
    exploded = explode_consequence_terms(d)
    # remove NA explícito
    exploded = exploded[exploded["ConsequenceTerm"] != "NA"]
    if exploded.empty:
        st.info("Sem termos de consequência codificantes para exibir.")
        return

    top_terms = (
        exploded.groupby("ConsequenceTerm")
        .size()
        .sort_values(ascending=False)
        .head(top_n)
        .index
    )
    d2 = exploded[exploded["ConsequenceTerm"].isin(top_terms)]
    agg = d2.groupby(["ConsequenceTerm", "ClinicalSignificance"]).size().reset_index(name="count")

    chart = alt.Chart(agg).mark_bar().encode(
        x=alt.X("ConsequenceTerm:N", sort="-y", title="Consequence term"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        tooltip=["ConsequenceTerm", "ClinicalSignificance", "count"]
    ).properties(
        title=f"Top {top_n} consequence terms × Clinical Significance (codificantes)",
        height=350
    ).interactive()

    st.altair_chart(chart, use_container_width=True)

def plot_impact_by_significance(st, df: pd.DataFrame):
    d = _only_coding(df)
    if d.empty:
        st.info("Sem variantes codificantes para ‘Impact’.")
        return
    d["Impact"] = d["Impact"].fillna("NA").astype(str).str.upper()
    agg = d.groupby(["Impact", "ClinicalSignificance"]).size().reset_index(name="count")

    chart = alt.Chart(agg).mark_bar().encode(
        x=alt.X("Impact:N", sort="-y", title="Impact"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        tooltip=["Impact", "ClinicalSignificance", "count"]
    ).properties(title="Impact × Clinical Significance (somente codificantes)", height=300)

    st.altair_chart(chart, use_container_width=True)

def plot_gene_hotspots(st, df: pd.DataFrame, bin_bp: int = 100):
    if bin_bp <= 0:
        bin_bp = 100

    d = df.dropna(subset=["Start", "GeneSymbol"]).copy()
    d["Start"] = _coerce_numeric(d["Start"])
    d = d.dropna(subset=["Start"])
    if d.empty:
        st.info("Sem coordenadas válidas ‘Start’ para hotspots no gene.")
        return

    # Bins em bp
    d["bin"] = (d["Start"] // bin_bp) * bin_bp
    # Para leitura melhor no eixo, converte bin para Mb
    d["bin_Mb"] = d["bin"] / 1e6
    d["bin_range"] = (d["bin"]).astype("Int64").astype(str) + "–" + (d["bin"] + bin_bp - 1).astype("Int64").astype(str)

    agg = (
        d.groupby(["GeneSymbol", "bin", "bin_Mb", "bin_range", "ClinicalSignificance"])
         .size()
         .reset_index(name="count")
    )

    base = alt.Chart(agg).mark_bar().encode(
        x=alt.X("bin_Mb:Q", title=f"Posição inicial do bin no gene (Mb, bin={bin_bp} bp)"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        tooltip=["GeneSymbol","bin_range","ClinicalSignificance","count"]
    ).properties(title="Hotspots por posição no gene (todas as variantes)", height=320)

    chart = base.facet(row=alt.Row("GeneSymbol:N", title=None)).resolve_scale(x='independent')

    st.altair_chart(chart, use_container_width=True)

def plot_top_alleles(st, df: pd.DataFrame, top_n: int = 20):
    d = df[["ReferenceAllele","AlternateAllele","ClinicalSignificance"]].copy()
    d["pair"] = d["ReferenceAllele"].fillna("NA") + "→" + d["AlternateAllele"].fillna("NA")
    agg = d.groupby(["pair","ClinicalSignificance"]).size().reset_index(name="count")
    top_pairs = agg.groupby("pair")["count"].sum().sort_values(ascending=False).head(top_n).index
    agg_top = agg[agg["pair"].isin(top_pairs)]

    chart = alt.Chart(agg_top).mark_bar().encode(
        x=alt.X("pair:N", sort="-y", title="Par Ref→Alt"),
        y=alt.Y("count:Q", title="Contagem"),
        color=alt.Color("ClinicalSignificance:N", title="Clinical Significance"),
        tooltip=["pair","ClinicalSignificance","count"]
    ).properties(
        title=f"Top {top_n} pares de alelos por Clinical Significance", height=350
    ).interactive()

    st.altair_chart(chart, use_container_width=True)
