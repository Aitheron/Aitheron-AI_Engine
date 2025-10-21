from pathlib import Path
import re

input = Path("patient.txt")      # ajuste o nome
output = Path("patient.fasta")
raw = input.read_text().strip()

# mantém apenas letras ACGTURYKMSWBDHVN (qualquer case); remove espaços/números
seq = re.sub(r'[^ACGTURYKMSWBDHVNacgturykmswbdhvn]', '', raw).upper()
assert seq, "Sequência vazia após sanitização."

with output.open('w') as f:
    f.write(">sample_BRCA?\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")

print(f"FASTA escrito: {output} | comprimento: {len(seq)} nt")
