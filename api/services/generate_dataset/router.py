from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pathlib import Path
import os
import pandas as pd

from services.generate_dataset.schemas import DatasetItem
from services.generate_dataset.generate_dataset_service import run_generate_dataset_service

router = APIRouter(prefix="/api/datasets", tags=["Gerar Datasets"])

@router.post("/download", response_class=FileResponse)
def download_excel(item: DatasetItem):

    is_both = item.both or (len(set(item.genes)) == 2)

    try:
        res = run_generate_dataset_service(
            genes=item.genes,
            force=item.force,
            both=is_both
        )

        excel_mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"

        if is_both:
            path = Path(res["merged_excel"])
            if not path or not path.exists():
                raise HTTPException(status_code=500, detail="Merged Excel not found.")
            # nome amigável e determinístico
            filename = "dataset_BRCA1_BRCA2.xlsx"
            return FileResponse(path, filename=filename, media_type=excel_mime)

        gene = item.genes[0]
        csv_path = Path(res["by_gene"][gene]["with_uniprot"])
        if not csv_path.exists():
            raise HTTPException(status_code=500, detail=f"CSV for {gene} not found.")

        excel_path = csv_path.with_suffix(".xlsx")
        if item.force or not excel_path.exists():
            df = pd.read_csv(csv_path)
            tmp = excel_path.with_suffix(".xlsx.tmp")
            with pd.ExcelWriter(tmp, engine="openpyxl") as w:
                df.to_excel(w, index=False, sheet_name="data")
            os.replace(tmp, excel_path)

        return FileResponse(excel_path, filename=excel_path.name, media_type=excel_mime)

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
