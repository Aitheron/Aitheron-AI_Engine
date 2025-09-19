from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from services.generate_dataset.router import router as datasets_router

API_PREFIX = "/api"

def create_app() -> FastAPI:
    app = FastAPI(title="BRCA API", version="0.1.0")

    @app.get(f"{API_PREFIX}/health")
    def health():
        return {"status": "ok"}

    app.include_router(datasets_router)

    return app

app = create_app()
