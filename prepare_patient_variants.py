import json, time
import requests
import pandas as pd

from core.settings import IMPACT_RANK, TERM_TO_FLAG

VEP_URL = "https://rest.ensembl.org/vep/human/region"
HEADERS_JSON = {"Content-Type":"application/json","Accept":"application/json"}

# termos que caracterizam região codificante (para setar is_coding=1)
CODING_SET = {
    "missense_variant","synonymous_variant","frameshift_variant",
    "inframe_insertion","inframe_deletion","start_lost","stop_gained","stop_lost",
    "protein_altering_variant","coding_sequence_variant","stop_retained_variant"
}

# termos que obrigam is_coding=0 quando são os mais severos
NONCODING_OVERRIDE = {
    "intron_variant","upstream_gene_variant","downstream_gene_variant",
    "5_prime_UTR_variant","3_prime_UTR_variant","regulatory_region_variant",
    "intergenic_variant"
}

SEQ_URL_TMPL = "https://rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1"
HEADERS_SEQ = {"Accept":"text/plain"}

def fetch_base(chrom,pos,cache):
    key=(chrom,pos)
    if key in cache: return cache[key]
    url=SEQ_URL_TMPL.format(chrom=chrom,start=pos,end=pos)
    r=requests.get(url,headers=HEADERS_SEQ,timeout=60)
    r.raise_for_status()
    base=r.text.strip().upper()
    cache[key]=base
    return base

def to_vcf_record(row, base_cache):
    chrom=str(row["Chromosome"])
    t=str(row["Type"]).strip()
    ref=str(row["Ref"])
    alt=str(row["Alt"])
    start=int(row["Start"])

    if t in ("SNV","Substituição"):
        pos=start
        vref=ref
        valt=alt
    elif t in ("Deleção","Deletion"):
        pos=max(1,start-1)
        anchor=fetch_base(chrom,pos,base_cache)
        vref=anchor + ("" if ref=="-" else ref)
        valt=anchor
    elif t in ("Inserção","Insertion","Duplicação","Duplication"):
        pos=start
        anchor=fetch_base(chrom,pos,base_cache)
        vref=anchor
        valt=anchor + ("" if alt=="-" else alt)
    elif t=="Indel":
        pos=max(1,start-1)
        anchor=fetch_base(chrom,pos,base_cache)
        vref=anchor + ("" if ref=="-" else ref)
        valt=anchor + ("" if alt=="-" else alt)
    else:
        pos=start
        vref=ref if ref!="-" else "N"
        valt=alt if alt!="-" else "N"

    return chrom,vref,valt,pos

def vep_call(vcf_like):
    payload={"variants":vcf_like,"canonical":1,"merged":1,"protein":1}
    r=requests.post(VEP_URL,headers=HEADERS_JSON,data=json.dumps(payload),timeout=180)
    if r.status_code==429:
        time.sleep(1.0)
        r=requests.post(VEP_URL,headers=HEADERS_JSON,data=json.dumps(payload),timeout=180)
    r.raise_for_status()
    return r.json()

def pick_best_transcript(rec, most):
    tcs=rec.get("transcript_consequences") or []
    c1=[tc for tc in tcs if (tc.get("biotype")=="protein_coding" and most in (tc.get("consequence_terms") or []))]
    if c1: return c1[0]
    c2=[tc for tc in tcs if (tc.get("canonical")=="YES" and tc.get("biotype")=="protein_coding")]
    if c2: return c2[0]
    c3=[tc for tc in tcs if tc.get("biotype")=="protein_coding"]
    if c3: return c3[0]
    c4=[tc for tc in tcs if tc.get("canonical")=="YES"]
    if c4: return c4[0]
    return tcs[0] if tcs else None

def extract_term_flags(rec, chosen_tc):
    # marca flags com base no transcrito escolhido; se não houver, usa união de todos
    flags = { col: 0 for col in TERM_TO_FLAG.values() }
    if chosen_tc and isinstance(chosen_tc.get("consequence_terms"), list):
        observed = set(chosen_tc["consequence_terms"])
    else:
        observed = set()
        for tc in rec.get("transcript_consequences", []) or []:
            for term in tc.get("consequence_terms", []) or []:
                if isinstance(term, str):
                    observed.add(term)
    for term, col in TERM_TO_FLAG.items():
        if term in observed:
            flags[col] = 1
    return flags, observed

def decide_is_coding(rec, most, chosen_tc, observed_terms):
    protein_coding = bool(chosen_tc and chosen_tc.get("biotype")=="protein_coding")
    has_coding_term = any((t in CODING_SET) for t in (observed_terms or []))
    is_coding = 1 if (protein_coding and has_coding_term) else 0
    if most in NONCODING_OVERRIDE:
        is_coding = 0
    return is_coding

def enrich_with_vep_to_csv(in_csv, out_csv):
    df=pd.read_csv(in_csv)
    base_cache={}
    vcf_variants=[]; mapping=[]
    for i,row in df.iterrows():
        chrom,vref,valt,pos=to_vcf_record(row,base_cache)
        vcf_variants.append(f"{chrom} {pos} . {vref} {valt}")
        mapping.append(i)

    def batched(seq,n=200):
        for i in range(0,len(seq),n):
            yield i,seq[i:i+n]

    annos=[None]*len(df)
    for i0,chunk in batched(vcf_variants,200):
        resp=vep_call(chunk)
        for j,rec in enumerate(resp):
            annos[mapping[i0+j]]=rec
        time.sleep(0.15)

    # garantir que colunas Is* existam mesmo se não aparecerem (preenchemos ao final)
    is_cols = list(TERM_TO_FLAG.values())

    records=[]
    for i,row in df.iterrows():
        rec=annos[i] or {}
        most=rec.get("most_severe_consequence")

        # escolhe o transcrito antes de extrair flags
        tc=pick_best_transcript(rec,most) if most else None

        # flags Is* a partir do transcrito escolhido (fallback: união)
        term_flags, observed_terms = extract_term_flags(rec, tc)

        # impacto
        impact_str=None; impact_num=None
        if tc and tc.get("impact"):
            impact_str=tc["impact"]
        if not impact_str:
            for _tc in rec.get("transcript_consequences",[]):
                if _tc.get("impact"):
                    impact_str=_tc["impact"]; break
        if impact_str:
            impact_num=IMPACT_RANK.get(impact_str,0)

        is_coding = decide_is_coding(rec, most, tc, observed_terms)

        # comprimentos
        len_ref=len(str(row["Ref"]).replace("-",""))
        len_alt=len(str(row["Alt"]).replace("-",""))
        length_change=len_alt-len_ref

        rec_out=dict(row)
        rec_out.update({
            "ConsequenceTerms":most,
            "impact":impact_str,
            "impact_rank":impact_num,
            "IsCoding": is_coding,
            "Protein_Start":tc.get("protein_start") if tc else None,
            "Protein_End":tc.get("protein_end") if tc else None,
            "CDNAStart":tc.get("cdna_start") if tc else None,
            "CDNAEnd":tc.get("cdna_end") if tc else None,
            "Codons":tc.get("codons") if tc else None,
            "len_ref":len_ref,
            "len_alt":len_alt,
            "length_change":length_change,
        })
        rec_out.update(term_flags)

        records.append(rec_out)

    out=pd.DataFrame(records)

    for c in is_cols:
        if c not in out.columns:
            out[c]=0

    out.to_csv(out_csv,index=False)
    return out

if __name__=="__main__":
    in_csv="patient_variants.csv"
    out_csv="patient_variants_annotated.csv"
    df=enrich_with_vep_to_csv(in_csv,out_csv)
    print(f"Anotadas {len(df)} variantes -> {out_csv}")
