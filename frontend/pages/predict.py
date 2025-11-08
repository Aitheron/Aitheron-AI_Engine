import streamlit as st
from pathlib import Path
import pandas as pd
from io import BytesIO
from stmol import showmol
import py3Dmol

from _utils import (
    sidebar,
    parse_single_fasta,
    validate_sequence,
    ALLOWED
)
from client.api_client import APIClient

st.set_page_config(page_title="Predição", page_icon="🧬", layout="wide")
sidebar()

CLASS_COL = "Classificação Clinica"
CODING_FLAG_COL = "Ocorreu em região codificante ?"
PROT_START = "Protein_Start"
PROT_END = "Protein_End"
CONSEQUENCE_COL = "Efeito no transcrito"

BENIGN_SET = {"Benigno"}
LIKELY_BENIGN_SET = {"Possivelmente Benigno"}
VUS_SET = {"VUS"}
PATH_SET = {"Patogênico"}

COLOR_BY_CLASS = {
    "Benigno": "green",
    "Provavelmente benigno": "lime",
    "VUS": "yellow",
    "Provavelmente patogênico": "orange",
    "Patogênico": "red"
}

FAINT_BG = {
    "Benigno": "#E6F4EA",
    "Provavelmente benigno": "#F0FAE6",
    "VUS": "#FFF9C4",
    "Provavelmente patogênico": "#FFF3E0",
    "Patogênico": "#FFEBEE",
    "default": "#F5F5F5"
}

GENE_TO_CIF = {
    "BRCA1": "../3d_models/AF_AFP38398F1.cif",
    "BRCA2": ""
}

@st.cache_data(show_spinner=False, ttl=24*3600)
def load_text_from_path(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return f.read()

def classify_from_label(row):
    raw = str(row.get(CLASS_COL, "")).strip()
    if raw in PATH_SET:
        return "Patogênico"
    if raw in BENIGN_SET:
        return "Benigno"
    if raw in LIKELY_BENIGN_SET:
        return "Provavelmente benigno"
    if raw in VUS_SET:
        return "VUS"
    return "VUS"

def prepare_dataframe(df):
    df = df.copy()
    if CODING_FLAG_COL in df.columns:
        coding_mask = df[CODING_FLAG_COL].astype(str).str.lower().str.startswith("s")
    else:
        coding_mask = df[[PROT_START, PROT_END]].notna().all(axis=1)
    df["_ClasseRender"] = df.apply(lambda r: classify_from_label(r), axis=1)
    df["_CorRender"] = df["_ClasseRender"].map(lambda x: COLOR_BY_CLASS.get(x, "magenta"))
    df_coding = df[coding_mask & df[[PROT_START, PROT_END]].notna().all(axis=1)].copy().reset_index(drop=True)
    df_non = df[~(coding_mask & df[[PROT_START, PROT_END]].notna().all(axis=1))].copy().reset_index(drop=True)
    return df, df_coding, df_non

def short_tag(row):
    cons = str(row.get(CONSEQUENCE_COL, "")).strip()
    cls = str(row.get("_ClasseRender", "")).strip()
    if cons and cls:
        return f"{cons} · {cls}"
    if cons:
        return cons
    return cls if cls else "var"

def build_records(df_coding):
    recs = []
    for i, r in df_coding.iterrows():
        try:
            ps, pe = int(r[PROT_START]), int(r[PROT_END])
        except:
            continue
        start, end = min(ps, pe), max(ps, pe)
        recs.append({
            "idx": i,
            "start": start,
            "end": end,
            "class": r["_ClasseRender"],
            "color": r["_CorRender"],
            "tag": short_tag(r),
        })
    return recs

def stagger_offset(i):
    k = i % 6
    offsets = [
        {"x": 12, "y": -10},
        {"x": 16, "y": 10},
        {"x": 22, "y": -16},
        {"x": -14, "y": 8},
        {"x": -18, "y": -12},
        {"x": 10, "y": 14},
    ]
    return offsets[k]

def expand_range(start, end, pad=2):
    a = max(1, start - pad)
    b = end + pad
    return list(range(a, b+1))

def render_structure(coords_str, fmt, df_coding, records, focus_idx=None, width=1000, height=860):
    v = py3Dmol.view(width=width, height=height)
    v.addModel(coords_str, fmt)
    v.setStyle({"cartoon": {"color": "spectrum"}})
    for i, p in enumerate(records):
        if focus_idx is not None and p["idx"] == focus_idx:
            continue
        sel_rng = {"resi": list(range(p["start"], p["end"]+1))}
        v.addStyle(sel_rng, {"stick": {"radius": 0.22, "color": p["color"]}})
        v.addStyle(sel_rng, {"sphere": {"radius": 0.55, "color": p["color"]}})
        off = stagger_offset(i)
        v.addLabel(p["tag"], {
            "backgroundOpacity": 0.65,
            "fontSize": 10,
            "alignment": "left",
            "showBackground": True,
            "borderThickness": 0.5,
            "screenOffset": off
        }, {"resi": [p["start"]]})
    if focus_idx is not None and 0 <= focus_idx < len(df_coding):
        row = df_coding.iloc[focus_idx]
        color = row["_CorRender"]
        start, end = int(row[PROT_START]), int(row[PROT_END])
        start, end = min(start, end), max(start, end)
        sel_focus = {"resi": expand_range(start, end, pad=2)}
        v.addStyle(sel_focus, {"stick": {"radius": 0.45, "color": color}})
        v.addStyle(sel_focus, {"sphere": {"radius": 1.05, "color": color}})
        v.addLabel(short_tag(row), {
            "backgroundOpacity": 0.85,
            "fontSize": 12,
            "alignment": "left",
            "showBackground": True,
            "borderThickness": 1.2,
            "screenOffset": {"x": 2, "y": -18}
        }, {"resi": [start]})
        v.zoomTo(sel_focus)
        v.zoom(1.08)
    else:
        v.zoomTo()
    showmol(v, height=height, width=width)

def reset_state():
    for k in (
        "header",
        "seq",
        "gene",
        "payload_preview",
        "api_result",
        "result_df",
        "result_bytes",
        "result_filename",
        "df_raw",
        "df_coding",
        "df_non",
        "records",
        "focus_idx",
        "cif_text"
    ):
        st.session_state.pop(k, None)

@st.cache_data(show_spinner=False)
def load_fasta_bytes():
    fasta_path = Path(__file__).resolve().parents[2] / "patient.fasta"
    return fasta_path.read_bytes()

def on_select_change():
    choice = st.session_state.get("variant_choice", "(nenhuma)")
    if choice == "(nenhuma)":
        st.session_state["focus_idx"] = None
    else:
        try:
            st.session_state["focus_idx"] = int(choice.split("|")[0].strip().replace("#",""))
        except:
            st.session_state["focus_idx"] = None

st.title("Upload de FASTA para Predição")
top_left, top_right = st.columns([7,3])

with top_right:
    st.info(
        """
        **Observações:**\n\n
        - Este fluxo **não faz autodetecção** de gene.\n
        - Envie **apenas um** registro FASTA por vez.\n
        - Se enviar FASTA de BRCA1 com gene **BRCA2** selecionado (ou vice-versa), **não há correção automática**.
        """
    )
    st.download_button(
        label="📥 Baixar .fasta (Exemplo)",
        data=load_fasta_bytes(),
        file_name="example.fasta",
        mime="text/plain",
        type="secondary",
        use_container_width=False,
    )

with top_left:
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
        btn = st.button("📤 Enviar para processamento", type="secondary", use_container_width=True)
        if btn:
            client = APIClient()
            try:
                with st.spinner("Processando..."):
                    xlsx_bytes = client.predict_from_fasta_file(
                        gene=st.session_state["gene"],
                        file_name=uploaded.name,
                        file_bytes=file_bytes,
                    )
                st.session_state["result_bytes"] = xlsx_bytes
                st.success("Processamento concluído.")
            except Exception as e:
                st.error(f"Falha na chamada de API: {e}")

if "result_bytes" in st.session_state:
    df = pd.read_excel(BytesIO(st.session_state["result_bytes"]), dtype=str)
    st.session_state["result_df"] = df
    st.session_state["df_raw"], st.session_state["df_coding"], st.session_state["df_non"] = prepare_dataframe(df)
    st.session_state["records"] = build_records(st.session_state["df_coding"])
    try:
        cif_path = GENE_TO_CIF.get(st.session_state.get("gene"), "")
        st.session_state["cif_text"] = load_text_from_path(cif_path) if cif_path else None
    except Exception as e:
        st.error(f"Falha ao carregar modelo 3D: {e}")
    viz_left, viz_right = st.columns([4,2], gap="large")
    with viz_left:
        if st.session_state.get("cif_text") and st.session_state["df_coding"] is not None:
            render_structure(st.session_state["cif_text"], "cif",
                             st.session_state["df_coding"],
                             st.session_state["records"],
                             focus_idx=st.session_state.get("focus_idx"))
        else:
            st.info("Modelo 3D indisponível.")
    with viz_right:
        st.markdown("#### Seleção de variante")
        if st.session_state["records"]:
            options = ["(nenhuma)"] + [f'#{r["idx"]} | {r["tag"]}' for r in st.session_state["records"]]
            st.selectbox("Escolha uma variante", options, index=0, key="variant_choice", on_change=on_select_change)
        st.markdown("#### Detalhes")
        if st.session_state.get("focus_idx") is not None and st.session_state.get("df_coding") is not None:
            idx = st.session_state["focus_idx"]
            dfc = st.session_state["df_coding"]
            if 0 <= idx < len(dfc):
                row = dfc.iloc[idx].copy()
                for c in ["_ClasseRender", "_CorRender"]:
                    if c in row.index:
                        row.drop(labels=c, inplace=True)
                gene_keys = ["Gene", "GENE", "Gene Symbol", "Symbol"]
                chrom_keys = ["Chromossomo", "Cromossomo", "Chromosome", "Chr", "Chrom"]
                gcol = next((k for k in gene_keys if k in row.index), None)
                ccol = next((k for k in chrom_keys if k in row.index), None)
                order = list(row.index)
                if gcol and ccol and gcol in order and ccol in order:
                    order.remove(gcol)
                    gi = order.index(ccol)
                    order.insert(gi, gcol)
                row = row[order]
                det_df = row.to_frame(name="Valor")
                if CLASS_COL in det_df.index:
                    cls_value = str(det_df.loc[CLASS_COL, "Valor"])
                    bg = FAINT_BG.get(cls_value, FAINT_BG["default"])
                    styler = det_df.style.set_properties(subset=pd.IndexSlice[[CLASS_COL], :], **{"background-color": bg})
                    st.table(styler)
                else:
                    st.dataframe(det_df, use_container_width=True)
            else:
                st.write("Seleção inválida.")
        else:
            st.caption("Selecione uma variante para ver detalhes e focar no 3D.")
    st.markdown("---")
    st.markdown("#### Variantes não codificantes (não mapeiam para a proteína)")
    df_non_final = st.session_state["df_non"].copy()
    for c in ["_ClasseRender", "_CorRender"]:
        if c in df_non_final.columns:
            df_non_final.drop(columns=c, inplace=True, errors="ignore")
    if len(df_non_final):
        st.dataframe(df_non_final, use_container_width=True)
    else:
        st.write("Sem variantes não codificantes no arquivo.")
    btn_col, _ = st.columns([2, 4])
    with btn_col:
        file_name = (
            st.session_state.get("result_filename")
            or f"predicoes_{st.session_state.get('gene')}.xlsx"
        )
        st.download_button(
            label="📥 Baixar Arquivo com predições",
            data=st.session_state["result_bytes"],
            file_name=file_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            type="secondary",
            use_container_width=False
        )
else:
    st.info("Aguardando o upload de um arquivo FASTA.")
