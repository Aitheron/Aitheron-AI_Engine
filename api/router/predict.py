import os
import numpy as np
import pandas as pd
from io import BytesIO
from tempfile import NamedTemporaryFile
from typing import Literal
from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from starlette.responses import StreamingResponse

from core.settings import ALLOWED_GENES
from services.predict_variants_service import run_predict_variants_service
from logger import get_logger

logger = get_logger(__name__)

router = APIRouter(prefix="/api/ml", tags=["Classificar Variantes Genéticas"])

@router.post("/predict")
async def predict(
    gene: Literal["BRCA1", "BRCA2"] = Form(...),
    file: UploadFile = File(...),
):

    logger.info("Fasta receveid, starting inference process !")

    if gene not in ALLOWED_GENES:
        logger.error(f"Invalid gene: {gene}. Options: {ALLOWED_GENES}")
        raise HTTPException(
            status_code=422,
            detail=f"Invalid gene: {gene}. Options: {ALLOWED_GENES}"
        )

    data = await file.read()
    if len(data) == 0:
        logger.error("File should not be empty !")
        raise HTTPException(422, detail="File should not be empty.")

    if not file.filename.endswith((".fasta", ".fa", ".fna")):
        logger.error("Invalid file extension. Please upload a FASTA file (.fasta, .fa, .fna) !")
        raise HTTPException(
            422,
            detail="Invalid file extension. Please upload a FASTA file (.fasta, .fa, .fna)"
        )

    with NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
        tmp.write(data)
        tmp_path = tmp.name

    try:
        df = run_predict_variants_service(
            gene=gene,
            file_path=tmp_path
        )
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    if not isinstance(df, pd.DataFrame):
        logger.error("Internal error: prediction did not return a DataFrame !")
        raise HTTPException(500, detail="Internal error: prediction did not return a DataFrame.")

    df = df.replace([np.inf, -np.inf], np.nan)

    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="data")
    output.seek(0)

    filename = f"prediction_{gene}.xlsx"
    return StreamingResponse(
        output,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
