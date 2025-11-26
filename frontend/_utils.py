import streamlit as st

MIN_LEN = 1000
ALLOWED = set("ACGTN")

def sidebar():
    with st.sidebar:
        st.markdown("## 🧬 Análise BRCA")
        st.caption("Navegue pelas páginas do app:")
        st.page_link("app.py", label="Início", icon="🏠")
        st.page_link("pages/datasets.py", label="Datasets", icon="📊")
        st.page_link("pages/predict.py", label="Predição", icon="🔮")
        st.page_link("pages/results.py", label="Resultados", icon="📑")
        st.divider()

def parse_single_fasta(file_bytes: bytes):
    try:
        text = file_bytes.decode(errors="ignore")
    except Exception:
        return None, None, "Falha ao decodificar o arquivo (não é texto)."
    header = None
    seq_parts = []
    records = 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            records += 1
            if records > 1:
                return None, None, "Foi encontrado mais de um registro FASTA; envie apenas um por vez."
            header = line[1:].strip()
        else:
            seq_parts.append(line.replace(" ", ""))
    if records == 0:
        return None, None, "Formato inválido: não foi encontrado cabeçalho FASTA iniciando com '>'."
    seq = "".join(seq_parts).upper()
    return header, seq, None

def validate_sequence(seq: str):
    if len(seq) < MIN_LEN:
        return False, f"Sequência muito curta ({len(seq)} bp). Mínimo exigido: {MIN_LEN} bp."
    invalid = sorted({ch for ch in seq if ch not in ALLOWED})
    if invalid:
        inv = "".join(invalid)
        return False, f"Caracteres inválidos na sequência: {inv}. Permitidos: {''.join(sorted(ALLOWED))}."
    return True, None