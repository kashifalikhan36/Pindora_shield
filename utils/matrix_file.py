import joblib
import numpy as np
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
except Exception as e:
    raise ImportError("RDKit is required for matrix_file. Install rdkit and try again.") from e

from typing import Any, Dict, Optional


class MatrixPredictor:

    def __init__(self, model_dir: Optional[str | Path] = None):
        self.model_dir = Path(model_dir) if model_dir else Path(__file__).parent
        self.ic50_bundle = self._load_bundle("catboost_ic50.pkl")
        self.assoc_bundle = self._load_bundle("catboost_association.pkl")
        self.phase_bundle = self._load_bundle("catboost_max_phase.pkl")
        self.target_bundle = self._load_bundle("catboost_target.pkl")

        self.radius = self.ic50_bundle["radius"]
        self.fp_size = self.ic50_bundle["fp_size"]
        self.morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=self.radius,
            fpSize=self.fp_size
        )

    def _load_bundle(self, filename: str) -> Dict[str, Any]:
        path = self.model_dir / filename
        if not path.exists():
            raise FileNotFoundError(
                f"Required model bundle not found: {path}. Pass the correct model_dir to MatrixPredictor() or place the bundle in this directory."
            )
        return joblib.load(path)

    def smiles_to_fp(self, smiles: str) -> np.ndarray:
        """Convert SMILES to fingerprint array using the configured Morgan generator."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        fp = self.morgan_gen.GetFingerprint(mol)
        return np.array(fp).reshape(1, -1)

    def predict_ic50(self, smiles: str) -> Dict[str, float]:
        X = self.smiles_to_fp(smiles)
        log_ic50 = self.ic50_bundle["model"].predict(X)[0]
        ic50 = 10 ** log_ic50
        return {"IC50": float(ic50), "log10_IC50": float(log_ic50)}

    def predict_association(self, smiles: str) -> Dict[str, float]:
        X = self.smiles_to_fp(smiles)
        log_assoc = self.assoc_bundle["model"].predict(X)[0]
        assoc = 10 ** log_assoc
        return {"Association_Score": float(assoc), "log10_Association": float(log_assoc)}

    def predict_phase(self, smiles: str) -> int:
        X = self.smiles_to_fp(smiles)
        max_phase = self.phase_bundle["model"].predict(X)[0]
        return int(max_phase)

    def predict_target(self, smiles: str) -> str:
        X = self.smiles_to_fp(smiles)
        target_idx = self.target_bundle["model"].predict(X)[0]
        target_symbol = self.target_bundle["label_encoder"].inverse_transform([int(target_idx)])[0]
        return str(target_symbol)

    def predict_all(self, smiles: str) -> Dict[str, Any]:
        """Run all predictions and return a consolidated dictionary (same shape as previous module-level version)."""
        X = self.smiles_to_fp(smiles)

        log_ic50 = self.ic50_bundle["model"].predict(X)[0]
        ic50 = 10 ** log_ic50

        log_assoc = self.assoc_bundle["model"].predict(X)[0]
        assoc = 10 ** log_assoc

        max_phase = self.phase_bundle["model"].predict(X)[0]

        target_idx = self.target_bundle["model"].predict(X)[0]
        target_symbol = self.target_bundle["label_encoder"].inverse_transform([int(target_idx)])[0]

        return {
            "IC50": float(ic50),
            "Association_Score": float(assoc),
            "Max_Clinical_Phase": int(max_phase),
            "Predicted_Target": str(target_symbol)
        }


if __name__ == "__main__":
    sample_smiles = "CC(=O)NC1=CC=C(O)C=C1"

    try:
        predictor = MatrixPredictor()
    except FileNotFoundError as exc:
        print("Model files not found. Skipping predictions.")
        print(str(exc))
    else:
        print("Running sample prediction for:", sample_smiles)
        results = predictor.predict_all(sample_smiles)
        for k, v in results.items():
            print(f"{k}: {v}")

        assert "IC50" in results and results["IC50"] >= 0
        assert "Association_Score" in results
        assert "Predicted_Target" in results
        print("Sample prediction completed successfully.")
