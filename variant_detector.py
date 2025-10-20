import csv
from typing import List, Dict, Tuple
import requests
import mappy as mp
from Bio import SeqIO

ENSEMBL_SEQ_REGION = "https://rest.ensembl.org/sequence/region/human/{region}"
HEADERS_FASTA = {"Accept": "text/plain"}

GENE_LOCI = {
    "BRCA1": {"chrom": "17", "start": 43044295, "end": 43125483, "strand": -1},
    "BRCA2": {"chrom": "13", "start": 32315474, "end": 32400266, "strand": +1},
}

TARGET_GENE = "BRCA1"
PATIENT_FASTA_PATH = "./patient.fasta"
OUT_CSV = "patient_variants.csv"
MINIMAP2_PRESET = "asm5"
MIN_MAPQ = 0
FLANK_BP = 500 # Adicona algumas bases adiconais para verificar existencia de mutações no começo e fim do gene
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

def parse_cs_to_variants(cs: str, ref_start_1based: int, ref_name: str) -> List[Dict]:
    i = 0
    pos = ref_start_1based
    variants = []
    mnp_ref = []
    mnp_alt = []
    mnp_start = None
    def flush_mnp(current_pos):
        nonlocal mnp_ref, mnp_alt, mnp_start
        if mnp_ref and mnp_ref != mnp_alt:
            ref = "".join(mnp_ref)
            alt = "".join(mnp_alt)
            end = current_pos - 1
            start = mnp_start
            variants.append({
                "Chromosome": ref_name,
                "Start": start,
                "End": end,
                "Ref": ref,
                "Alt": alt,
                "Type": "MNP" if len(ref) > 1 else "SNV",
                "Origin": "CS"
            })
        mnp_ref, mnp_alt, mnp_start = [], [], None
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
            if mnp_start is None:
                mnp_start = pos
            ref_base = cs[i+1]
            alt_base = cs[i+2]
            mnp_ref.append(ref_base.upper())
            mnp_alt.append(alt_base.upper())
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
            variants.append({
                "Chromosome": ref_name,
                "Start": start,
                "End": end,
                "Ref": deleted,
                "Alt": "-",
                "Type": "DEL",
                "Origin": "CS"
            })
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
            variants.append({
                "Chromosome": ref_name,
                "Start": anchor,
                "End": anchor,
                "Ref": "-",
                "Alt": inserted,
                "Type": "INS",
                "Origin": "CS"
            })
            i = j
            continue
        i += 1
    if mnp_ref:
        flush_mnp(pos)
    return variants

def call_structurals_from_splits(hits) -> List[Dict]:
    svs = []
    if not hits:
        return svs
    hits = sorted(hits, key=lambda h: (h.q_st, h.q_en))
    for a, b in zip(hits, hits[1:]):
        if a.ctg != b.ctg:
            continue
        if a.strand != b.strand:
            start = min(a.r_st, b.r_st) + 1
            end = max(a.r_en, b.r_en)
            svs.append({
                "Chromosome": a.ctg,
                "Start": start,
                "End": end,
                "Ref": "-",
                "Alt": "-",
                "Type": "INV",
                "Origin": "SPLIT"
            })
        else:
            gap = max(0, b.r_st - a.r_en)
            if gap <= 50:
                start = min(a.r_st, b.r_st) + 1
                end = max(a.r_en, b.r_en)
                svs.append({
                    "Chromosome": a.ctg,
                    "Start": start,
                    "End": end,
                    "Ref": "-",
                    "Alt": "-",
                    "Type": "DUP",
                    "Origin": "SPLIT"
                })
    return svs

def emit_softclip_variants(h, qseq: str, ref_seq: str, chrom: str, ref_global_start_1b: int) -> List[Dict]:
    vars_sc = []
    left_len = h.q_st
    right_len = len(qseq) - h.q_en

    if left_len >= 1:
        ins = qseq[:left_len].upper()
        if left_len == 1 and h.r_st > 0:
            rb = ref_seq[h.r_st - 1].upper()
            qb = ins
            if rb != qb:
                pos = ref_global_start_1b + h.r_st  # posição 1-based do rb
                vars_sc.append({
                    "Chromosome": chrom,
                    "Start": pos,
                    "End": pos,
                    "Ref": rb,
                    "Alt": qb,
                    "Type": "SNV",
                    "Origin": "left_clip_rescue"
                })
            # se rb == qb, ignora
        else:
            anchor = ref_global_start_1b + h.r_st
            vars_sc.append({
                "Chromosome": chrom,
                "Start": anchor,
                "End": anchor,
                "Ref": "-",
                "Alt": ins,
                "Type": "INS",
                "Origin": "left_clip"
            })

    if right_len >= 1:
        ins = qseq[-right_len:].upper()
        if right_len == 1 and h.r_en < len(ref_seq):
            rb = ref_seq[h.r_en].upper()
            qb = ins
            if rb != qb:
                pos = ref_global_start_1b + h.r_en + 1  # posição do rb
                vars_sc.append({
                    "Chromosome": chrom,
                    "Start": pos,
                    "End": pos,
                    "Ref": rb,
                    "Alt": qb,
                    "Type": "SNV",
                    "Origin": "right_clip_rescue"
                })
            # se rb == qb, ignora
        else:
            anchor = ref_global_start_1b + h.r_en
            vars_sc.append({
                "Chromosome": chrom,
                "Start": anchor,
                "End": anchor,
                "Ref": "-",
                "Alt": ins,
                "Type": "INS",
                "Origin": "right_clip"
            })

    if h.q_st < len(qseq) and h.r_st < len(ref_seq):
        rb = ref_seq[h.r_st].upper()
        qb = qseq[h.q_st].upper()
        if rb != qb:
            pos = ref_global_start_1b + h.r_st + 1
            vars_sc.append({
                "Chromosome": chrom,
                "Start": pos,
                "End": pos,
                "Ref": rb,
                "Alt": qb,
                "Type": "SNV",
                "Origin": "left_edge_snv"
            })

    if h.q_en > 0 and h.r_en - 1 < len(ref_seq) and h.q_en - 1 < len(qseq):
        rb = ref_seq[h.r_en - 1].upper()
        qb = qseq[h.q_en - 1].upper()
        if rb != qb:
            pos = ref_global_start_1b + h.r_en
            vars_sc.append({
                "Chromosome": chrom,
                "Start": pos,
                "End": pos,
                "Ref": rb,
                "Alt": qb,
                "Type": "SNV",
                "Origin": "right_edge_snv"
            })

    return vars_sc

def analyze_patient_to_csv(ref_seq: str, chrom: str, ref_global_start_1b: int, patient_fasta_path: str, out_csv: str, preset: str, min_mapq: int):
    al = mp.Aligner(seq=ref_seq, preset=preset)
    if not al:
        raise RuntimeError("Falha ao inicializar o minimap2 (mappy).")
    patient_seqs = read_fasta_string(patient_fasta_path)
    rows = []
    for _, pacient_seq in patient_seqs:
        hits = [h for h in al.map(pacient_seq, cs=True) if h.mapq >= min_mapq]
        for hit in hits:
            if hit.cs:
                ref_start_1b_seg = ref_global_start_1b + hit.r_st
                vars_seg = parse_cs_to_variants(hit.cs, ref_start_1b_seg + 1, chrom)
                for v in vars_seg:
                    rows.append(v)
            vars_sc = emit_softclip_variants(hit, pacient_seq, ref_seq, chrom, ref_global_start_1b + 0)
            for v in vars_sc:
                rows.append(v)
        svs = call_structurals_from_splits(hits)
        for s in svs:
            rows.append(s)
    rows.sort(key=lambda r: (r["Chromosome"], r["Start"], r["End"], r.get("Type","")))
    cols = ["Chromosome", "Start", "End", "Ref", "Alt", "Type"]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    return rows

def detect_variants_from_fasta():
    locus = GENE_LOCI[TARGET_GENE]
    chrom = locus["chrom"]
    start = max(1, locus["start"] - FLANK_BP)
    end = locus["end"] + FLANK_BP
    strand = locus["strand"]
    ref_seq = fetch_reference(chrom, start, end, strand)
    variants = analyze_patient_to_csv(
        ref_seq=ref_seq,
        chrom=chrom,
        ref_global_start_1b=start,
        patient_fasta_path=PATIENT_FASTA_PATH,
        out_csv=OUT_CSV,
        preset=MINIMAP2_PRESET,
        min_mapq=MIN_MAPQ
    )
    print(TARGET_GENE)
    print(f"{chrom}:{start}-{end} strand {strand} len={len(ref_seq)}")
    print(PATIENT_FASTA_PATH)
    print(OUT_CSV)
    print(f"Detected {len(variants)} variants")

if __name__ == "__main__":
    detect_variants_from_fasta()
