import pytest
from pydantic import ValidationError
from schemas.datasets import DatasetItem

def test_datasetitem_valid_defaults():
    item = DatasetItem(genes=["BRCA1"])
    assert item.genes == ["BRCA1"]
    assert item.force is False
    assert item.is_training_dt is False

def test_datasetitem_two_genes_ok():
    item = DatasetItem(genes=["BRCA1", "BRCA2"], is_training_dt=True)
    assert item.genes == ["BRCA1", "BRCA2"]
    assert item.is_training_dt is True

def test_datasetitem_min_length_enforced():
    with pytest.raises(ValidationError):
        DatasetItem(genes=[])

def test_datasetitem_max_length_enforced():
    with pytest.raises(ValidationError):
        DatasetItem(genes=["BRCA1", "BRCA2", "BRCA1"])

def test_datasetitem_literal_enforced():
    with pytest.raises(ValidationError):
        DatasetItem(genes=["TP53"])
