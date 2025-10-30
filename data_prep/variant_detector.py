import requests
import mappy as mp
import pandas as pd

from Bio import SeqIO
from typing import List, Dict, Tuple, Optional

ENSEMBL_SEQ_REGION = "https://rest.ensembl.org/sequence/region/human/{region}"
HEADERS_FASTA = {"Accept": "text/plain"}

GENE_LOCI = {
    "BRCA1": {"chrom": "17", "start": 43044295, "end": 43125483, "strand": 1},
    "BRCA2": {"chrom": "13", "start": 32315474, "end": 32400266, "strand": +1},
}

MINIMAP2_PRESET = "map-ont"
MIN_MAPQ = 0
FLANK_BP = 500
MIN_SOFTCLIP_INS = 1

def fetch_reference(chrom: str, start: int, end: int, strand: int) -> str:
    region = f"{chrom}:{start}..{end}:{strand}"
    url = ENSEMBL_SEQ_REGION.format(region=region)
    r = requests.get(url, headers=HEADERS_FASTA, timeout=60)
    r.raise_for_status()
    return r.text.strip().replace("\n", "")

def read_fasta_string(path: str) -> List[Tuple[str, str]]:
    recs = []
    for rec in SeqIO.parse(path, "fasta"):
        recs.append((rec.id, str(rec.seq).upper()))
    if not recs:
        raise ValueError("FASTA do paciente vazio.")
    return recs

def coalesce_microindels(events: List[Dict]) -> List[Dict]:
    ev = sorted(events, key=lambda e: (e["Start"], e["End"]))
    out = []
    i = 0
    while i < len(ev):
        a = ev[i]
        merged = False
        if i + 1 < len(ev):
            b = ev[i + 1]
            if a["Chromosome"] == b["Chromosome"]:
                gap = b["Start"] - a["End"]
                near = gap <= 1
                pair = {a["Type"], b["Type"]}
                if near and pair.issubset({"Deleção", "Inserção"}):
                    ref = a["Ref"] if a["Type"] == "Deleção" else b["Ref"]
                    alt = a["Alt"] if a["Type"] == "Inserção" else b["Alt"]
                    start = min(a["Start"], b["Start"])
                    end = max(a["End"], b["End"])
                    out.append({"Chromosome": a["Chromosome"], "Start": start, "End": end, "Ref": ref, "Alt": alt, "Type": "Indel", "Origin": "Coalesced"})
                    i += 2
                    merged = True
                elif near and (("SNV" in pair) and (("Inserção" in pair) or ("Deleção" in pair))):
                    if a["Type"] == "SNV" and b["Type"] == "Inserção":
                        ref = a["Ref"]
                        alt = a["Alt"] + b["Alt"]
                        start = a["Start"]
                        end = max(a["End"], b["End"])
                    elif a["Type"] == "SNV" and b["Type"] == "Deleção":
                        ref = a["Ref"] + b["Ref"]
                        alt = a["Alt"]
                        start = min(a["Start"], b["Start"])
                        end = max(a["End"], b["End"])
                    elif a["Type"] == "Inserção" and b["Type"] == "SNV":
                        ref = b["Ref"]
                        alt = a["Alt"] + b["Alt"]
                        start = min(a["Start"], b["Start"])
                        end = b["End"]
                    elif a["Type"] == "Deleção" and b["Type"] == "SNV":
                        ref = a["Ref"] + b["Ref"]
                        alt = b["Alt"]
                        start = a["Start"]
                        end = max(a["End"], b["End"])
                    else:
                        ref = a["Ref"]
                        alt = a["Alt"]
                        start = a["Start"]
                        end = a["End"]
                    out.append({"Chromosome": a["Chromosome"], "Start": start, "End": end, "Ref": ref, "Alt": alt, "Type": "Indel", "Origin": "Coalesced"})
                    i += 2
                    merged = True
        if not merged:
            out.append(a)
            i += 1
    return out

def parse_cs_to_variants(cs: str, ref_start_1based: int, ref_name: str, ref_seq: str, locus_start_1b: int) -> List[Dict]:
    i = 0
    pos = ref_start_1based
    variants = []
    mnp_ref = []
    mnp_alt = []
    mnp_pos = []
    def flush_mnp(current_pos):
        nonlocal mnp_ref, mnp_alt, mnp_pos
        if mnp_ref and mnp_alt and len(mnp_ref) == len(mnp_alt):
            for k in range(len(mnp_ref)):
                if mnp_ref[k] != mnp_alt[k]:
                    p = mnp_pos[k]
                    variants.append({"Chromosome": ref_name, "Start": p, "End": p, "Ref": mnp_ref[k], "Alt": mnp_alt[k], "Type": "SNV", "Origin": "CS"})
        mnp_ref, mnp_alt, mnp_pos = [], [], []
    while i < len(cs):
        c = cs[i]
        if c == ':':
            i += 1
            j = i
            while j < len(cs) and cs[j].isdigit():
                j += 1
            run = int(cs[i:j]) if j > i else 0
            if mnp_ref:
                flush_mnp(pos)
            pos += run
            i = j
            continue
        if c == '*':
            ref_base = cs[i+1]
            alt_base = cs[i+2]
            mnp_ref.append(ref_base.upper())
            mnp_alt.append(alt_base.upper())
            mnp_pos.append(pos)
            pos += 1
            i += 3
            continue
        if c == '-':
            if mnp_ref:
                flush_mnp(pos)
            i += 1
            j = i
            while j < len(cs) and cs[j] not in ':*+-':
                j += 1
            deleted = cs[i:j].upper()
            start = pos
            end = pos + len(deleted) - 1
            variants.append({"Chromosome": ref_name, "Start": start, "End": end, "Ref": deleted, "Alt": "-", "Type": "Deleção", "Origin": "CS"})
            pos += len(deleted)
            i = j
            continue
        if c == '+':
            if mnp_ref:
                flush_mnp(pos)
            i += 1
            j = i
            while j < len(cs) and cs[j] not in ':*+-':
                j += 1
            inserted = cs[i:j].upper()
            anchor = pos - 1 if pos > ref_start_1based else pos
            typ = "Inserção"
            end_idx = anchor - locus_start_1b
            start_idx = end_idx - len(inserted) + 1
            if start_idx >= 0 and end_idx >= 0:
                left = ref_seq[start_idx:end_idx+1].upper()
                right_start = end_idx + 1
                right_end = right_start + len(inserted)
                right = ref_seq[right_start:right_end].upper()

                if left == inserted or right == inserted:
                    typ = "Duplicação"
            variants.append({"Chromosome": ref_name, "Start": anchor, "End": anchor, "Ref": "-", "Alt": inserted, "Type": typ, "Origin": "CS"})
            i = j
            continue
        i += 1
    if mnp_ref:
        flush_mnp(pos)
    variants = coalesce_microindels(variants)
    return variants

def call_structurals_from_splits(hits) -> List[Dict]:
    svs = []
    if not hits:
        return svs
    hits = sorted(hits, key=lambda h: (h.q_st, h.q_en))
    for a, b in zip(hits, hits[1:]):
        if a.ctg != b.ctg:
            continue
        if a.strand == b.strand:
            gap = max(0, b.r_st - a.r_en)
            if gap <= 50:
                start = min(a.r_st, b.r_st) + 1
                end = max(a.r_en, b.r_en)
                svs.append({"Chromosome": a.ctg, "Start": start, "End": end, "Ref": "-", "Alt": "-", "Type": "Duplicação", "Origin": "SPLIT"})
    return svs

def emit_softclip_variants(h, qseq: str, ref_seq: str, chrom: str, ref_global_start_1b: int) -> List[Dict]:
    vars_sc = []
    left_len = h.q_st
    right_len = len(qseq) - h.q_en
    if left_len >= MIN_SOFTCLIP_INS:
        ins = qseq[:left_len].upper()
        if left_len == 1 and h.r_st > 0:
            rb = ref_seq[h.r_st - 1].upper()
            qb = ins
            if rb != qb:
                pos = ref_global_start_1b + h.r_st
                vars_sc.append({"Chromosome": chrom, "Start": pos, "End": pos, "Ref": rb, "Alt": qb, "Type": "SNV", "Origin": "left_clip_rescue"})
        else:
            anchor = ref_global_start_1b + h.r_st
            vars_sc.append({"Chromosome": chrom, "Start": anchor, "End": anchor, "Ref": "-", "Alt": ins, "Type": "Inserção", "Origin": "left_clip"})
    if right_len >= MIN_SOFTCLIP_INS:
        ins = qseq[-right_len:].upper()
        if right_len == 1 and h.r_en < len(ref_seq):
            rb = ref_seq[h.r_en].upper()
            qb = ins
            if rb != qb:
                pos = ref_global_start_1b + h.r_en + 1
                vars_sc.append({"Chromosome": chrom, "Start": pos, "End": pos, "Ref": rb, "Alt": qb, "Type": "SNV", "Origin": "right_clip_rescue"})
        else:
            anchor = ref_global_start_1b + h.r_en
            vars_sc.append({"Chromosome": chrom, "Start": anchor, "End": anchor, "Ref": "-", "Alt": ins, "Type": "Inserção", "Origin": "right_clip"})
    if h.q_st < len(qseq) and h.r_st < len(ref_seq):
        rb = ref_seq[h.r_st].upper()
        qb = qseq[h.q_st].upper()
        if rb != qb:
            pos = ref_global_start_1b + h.r_st + 1
            vars_sc.append({"Chromosome": chrom, "Start": pos, "End": pos, "Ref": rb, "Alt": qb, "Type": "SNV", "Origin": "left_edge_snv"})
    if h.q_en > 0 and h.r_en - 1 < len(ref_seq) and h.q_en - 1 < len(qseq):
        rb = ref_seq[h.r_en - 1].upper()
        qb = qseq[h.q_en - 1].upper()
        if rb != qb:
            pos = ref_global_start_1b + h.r_en
            vars_sc.append({"Chromosome": chrom, "Start": pos, "End": pos, "Ref": rb, "Alt": qb, "Type": "SNV", "Origin": "right_edge_snv"})
    return vars_sc

def analyze_patient_fasta(
    ref_seq: str,
    chrom: str,
    ref_global_start_1b: int,
    patient_fasta_path: str,
    preset: str,
    min_mapq: int,
    out_csv : Optional[str] = None
):
    al = mp.Aligner(seq=ref_seq, preset=preset)
    if not al:
        raise RuntimeError("Falha ao inicializar o minimap2 (mappy).")
    patient_seqs = read_fasta_string(patient_fasta_path)
    rows = []
    for _, pacient_seq in patient_seqs:
        hits = [h for h in al.map(pacient_seq, cs=True) if h.mapq >= min_mapq]
        for hit in hits:
            if hit.cs:
                seg_ref_start_1b = ref_global_start_1b + hit.r_st
                vars_seg = parse_cs_to_variants(hit.cs, seg_ref_start_1b, chrom, ref_seq, ref_global_start_1b)
                for v in vars_seg:
                    rows.append(v)
            vars_sc = emit_softclip_variants(hit, pacient_seq, ref_seq, chrom, ref_global_start_1b)
            for v in vars_sc:
                rows.append(v)
        svs = call_structurals_from_splits(hits)
        for s in svs:
            rows.append(s)
    rows.sort(key=lambda r: (r["Chromosome"], r["Start"], r["End"], r.get("Type","")))
    cols = ["Chromosome", "Start", "End", "Ref", "Alt", "Type", "Origin"]
    
    cols = ["Chromosome", "Start", "End", "Ref", "Alt", "Type", "Origin"]
    rows.sort(key=lambda r: (r["Chromosome"], r["Start"], r["End"], r.get("Type","")))
    df = pd.DataFrame(rows, columns=cols)

    if not df.empty:
        df["Start"] = pd.to_numeric(df["Start"], errors="coerce").astype("Int64")
        df["End"]   = pd.to_numeric(df["End"], errors="coerce").astype("Int64")
    
    if out_csv:
        df.to_csv(out_csv, index=False)

    return df

def detect_variants_from_fasta(
    target_gene : str,
    fasta_file_path : str
):
    locus = GENE_LOCI[target_gene]
    chrom = locus["chrom"]
    start = max(1, locus["start"] - FLANK_BP)
    end = locus["end"] + FLANK_BP
    strand = locus["strand"]
    ref_seq = fetch_reference(chrom, start, end, strand)
    variants = analyze_patient_fasta(
        ref_seq=ref_seq,
        chrom=chrom,
        ref_global_start_1b=start,
        patient_fasta_path=fasta_file_path,
        preset=MINIMAP2_PRESET,
        min_mapq=MIN_MAPQ
    )
    print(target_gene)
    print(f"{chrom}:{start}-{end} strand {strand} len={len(ref_seq)}")
    print(f"Detected {len(variants)} variants")

    return variants

if __name__ == "__main__":
    detect_variants_from_fasta(
        target_gene="BRCA1",
        fasta_file_path="./patient.fasta"
    )
