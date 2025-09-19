from pydantic import BaseModel, Field
from typing import List, Literal

Gene = Literal["BRCA1", "BRCA2"]

class DatasetItem(BaseModel):
    genes: List[Gene] = Field(..., min_length=1, max_length=2)
    force: bool = False
    both: bool = False
