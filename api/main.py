from fastapi import FastAPI

from router import (
    datasets_router,
    predict_router
)

API_PREFIX = "/api"

def create_app() -> FastAPI:
    app = FastAPI(title="BRCA API", version="0.1.0")

    @app.get(f"{API_PREFIX}/health")
    def health():
        return {"status": "ok"}

    app.include_router(datasets_router)
    app.include_router(predict_router)

    return app

app = create_app()
