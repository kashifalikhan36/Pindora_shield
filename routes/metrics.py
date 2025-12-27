from fastapi import APIRouter, HTTPException
from models.schemas import Generate3DInput
from utils.copilot import AzureOpenAIChatClient
from utils.matrix_file import MatrixPredictor

router = APIRouter(
    prefix="/metrics",
    tags=["pindora"],
    responses={404: {"description": "Not found"}},
)

@router.post("/metrics_data")
async def metrics_data(request: Generate3DInput):

    client = AzureOpenAIChatClient()
    matrix = MatrixPredictor()

    results = matrix.predict_all(request.input_smile)

    report = client.generate_report_from_smiles_ic50_value_association_score_target_symbol_max_phase(request.input_smile, results["IC50"], results["Association_Score"], results["Predicted_Target"], results["Max_Clinical_Phase"])
    
    if not request.input_smile or len(request.input_smile.strip()) == 0:
        raise HTTPException(status_code=400, detail="SMILES string is required")
    return {
        "report": report,
        "status": "success",
    }
