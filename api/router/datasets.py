from fastapi import APIRouter, HTTPException
from starlette.responses import StreamingResponse
from io import BytesIO
from pathlib import Path
import pandas as pd

from schemas.datasets import DatasetItem
from services.generate_dataset_service import run_generate_dataset_service
from logger import get_logger

logger = get_logger(__name__)

router = APIRouter(prefix="/api/datasets", tags=["Gerar Datasets"])

@router.post("/download")
def download_excel(item: DatasetItem):
    logger.info("Preparing dataset to download !")
    if not item.genes:
        logger.error("At least one gene must be provided !")
        raise HTTPException(status_code=400, detail="At least one gene must be provided.")
    genes_sorted = sorted(list(dict.fromkeys(item.genes)))
    kind_suffix = "training" if item.is_training_dt else "merged"
    filename = f"dataset_{'_'.join(genes_sorted)}_{kind_suffix}.xlsx"
    try:
        res = run_generate_dataset_service(
            genes=genes_sorted,
            force=False,
            is_training_dt=item.is_training_dt
        )
        csv_path = Path(res.get("dataset") or "")
        if not csv_path.exists():
            raise HTTPException(status_code=500, detail="CSV not found.")
        df = pd.read_csv(csv_path, dtype=str)
        output = BytesIO()
        with pd.ExcelWriter(output, engine="openpyxl") as writer:
            df.to_excel(writer, index=False, sheet_name="data")
        output.seek(0)
        logger.info("Dataset is ready for download !")
        return StreamingResponse(
            output,
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
