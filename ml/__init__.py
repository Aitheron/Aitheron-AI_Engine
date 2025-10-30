from functools import lru_cache

@lru_cache(maxsize=1)
def get_predictor():
    from ml.train._settings import Config
    from ml.predict import Predictor
    cfg = Config()
    return Predictor("../artifacts/final_model/model.pth", cfg)

__all__ = ["get_predictor"]
