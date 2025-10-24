import pandas as pd

def clean_data_per_submit_confidence(
    input_files : list[str],
    valid_confidence : list[str]
) -> list[str]:

    new_files_path = []
    for file_path in input_files:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
        
        confidence_columns = [c for c in df.columns 
                            if "review" in c.lower() or "criteria" in c.lower() or "assertion" in c.lower()]
        
        if not confidence_columns:
            raise ValueError(f"No confidence/review column found in: {file_path}")
        
        confidence_col = confidence_columns[0]  # take the first matching column
        df_filtered = df[df[confidence_col].str.lower().isin(valid_confidence)]
        
        output_path = file_path.replace(".tsv", "_filtered.tsv")
        new_files_path.append(output_path)
        df_filtered.to_csv(output_path, sep="\t", index=False)
        print(f"Original df: {len(df)} - Processed df: {len(df_filtered)}")
        print(f"Filtered file saved to: {output_path} ({len(df_filtered)} variants kept)")
    return new_files_path

if __name__ == "__main__":
    input_files = [
        "../../files/clinvar_BRCA1_GRCh38.tsv",
        "../../files/clinvar_BRCA2_GRCh38.tsv"
    ]
    # Terms that represent high-confidence variants (≥ 3 stars)
    
    valid_confidence = ["reviewed by expert panel", "practice guideline"]
    clean_data_per_submit_confidence(
        input_files=input_files,
        valid_confidence=valid_confidence
    )
