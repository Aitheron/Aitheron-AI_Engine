#!/usr/bin/env sh
cd api
PYTHONPATH="$(pwd)/.." uvicorn main:app --reload --host 0.0.0.0 --port 8000 &

cd ..
cd frontend
streamlit run app.py --server.address 0.0.0.0 --server.port 8501