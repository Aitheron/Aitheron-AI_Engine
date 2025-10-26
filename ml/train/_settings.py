from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class ModelFlags:
    use_mlp: bool = True
    use_xgb: bool = False
    use_rf: bool = False

@dataclass
class Paths:
    train_csv_path: str = "./filess/clinvar_BRCA1_BRCA2_GRCh38_training.csv"
    artifacts_dir: str = "./artifacts"

@dataclass
class Columns:
    label_col: str = "Label"
    gene_col: str = "GeneSymbol"

    numeric_cols: List[str] = field(default_factory=lambda: [
        "Start", "End",
        "len_ref", "len_alt", "length_change",
        "ImpactRank",
        "Protein_Start", "Protein_End",
        "CDNAStart", "CDNAEnd",
    ])

    cat_cols: List[str] = field(default_factory=list)

    # prefixo das colunas binárias (todas as Is* geradas + IsCoding)
    is_prefix: str = "Is"

@dataclass
class CV:
    n_splits: int = 10
    shuffle: bool = True
    random_state: int = 42

@dataclass
class TrainParams:
    hidden_dims: List[int] = field(default_factory=lambda: [256, 128])
    dropout: float = 0.2
    batch_size: int = 256
    epochs: int = 60
    lr: float = 1e-3
    weight_decay: float = 1e-5
    early_stopping_patience: int = 8
    device: Optional[str] = None  # "cuda" | "cpu" | None (auto)

@dataclass
class Config:
    # IMPORTANT: todos que são objetos devem usar default_factory
    flags: ModelFlags = field(default_factory=ModelFlags)
    paths: Paths = field(default_factory=Paths)
    cols: Columns = field(default_factory=Columns)
    cv: CV = field(default_factory=CV)
    train: TrainParams = field(default_factory=TrainParams)
