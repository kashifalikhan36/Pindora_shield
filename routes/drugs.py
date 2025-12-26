from fastapi import APIRouter, HTTPException
from models.schemas import TextInput, TextResponse, HealthCheck
from pindora import Pindora
import json

router = APIRouter(
    prefix="/api",
    tags=["pindora"],
    responses={404: {"description": "Not found"}},
)

@router.post("/drug_konsi_ohundni_h", response_model=TextResponse)
async def process_text(request: TextInput):
    if not request.text or len(request.text.strip()) == 0:
        raise HTTPException(status_code=400, detail="Text cannot be empty")

    pindora = Pindora()
    results = pindora.drug_discovery_pipeline(request.text)
      
    return {
        "results": results,
        "status": "success",
        }
