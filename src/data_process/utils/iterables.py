import pandas as pd

def read_data(file_path):
    data = pd.read_csv(file_path, sep='\t')
    
    return data

def tsv_to_xlsx(input_path: str, output_path: str = None):

    if output_path is None:
        output_path = input_path.rsplit(".", 1)[0] + ".xlsx"

    df = read_data(input_path)

    df.to_excel(output_path, index=False)

    return output_path

if __name__ == "__main__":
    tsv_to_xlsx('../files/clinvar_BRCA1_GRCh38_filtered.tsv')
    tsv_to_xlsx('../files/clinvar_BRCA2_GRCh38_filtered.tsv')