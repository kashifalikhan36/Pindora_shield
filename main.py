from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from routes.drugs import router as text_router

app = FastAPI(title="Pindora Shield API",description="Drug discovery and molecule generation API",version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(text_router)