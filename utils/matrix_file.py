import joblib
import numpy as np
from pathlib import Path
from typing import Any, Dict, Optional, Union

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


class MatrixPredictor:
    def __init__(self, model_dir: Optional[Union[str, Path]] = None):
        # Resolve model directory
        self.model_dir = Path(model_dir) if model_dir else Path(__file__).parent.parent

        # Load model bundles
        self.ic50_bundle = self._load_bundle("matriX_model/ic50.pkl")
        self.assoc_bundle = self._load_bundle("matriX_model/association.pkl")
        self.phase_bundle = self._load_bundle("matriX_model/max_phase.pkl")
        self.target_bundle = self._load_bundle("matriX_model/target.pkl")

        # Fingerprint config (single source of truth)
        self.radius = self.ic50_bundle["radius"]
        self.fp_size = self.ic50_bundle["fp_size"]

        self.morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=self.radius,
            fpSize=self.fp_size
        )

    def _load_bundle(self, relative_path: str) -> Dict[str, Any]:
        path = self.model_dir / relative_path
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")

        return joblib.load(path)

    def smiles_to_fp(self, smiles: str) -> np.ndarray:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")

        fp = self.morgan_gen.GetFingerprint(mol)
        return np.array(fp).reshape(1, -1)

    # ---------------- Predictions ---------------- #

    def predict_ic50(self, smiles: str) -> Dict[str, float]:
        X = self.smiles_to_fp(smiles)
        log_ic50 = self.ic50_bundle["model"].predict(X)[0]
        return {
            "IC50": float(10 ** log_ic50),
            "log10_IC50": float(log_ic50)
        }

    def predict_association(self, smiles: str) -> Dict[str, float]:
        X = self.smiles_to_fp(smiles)
        log_assoc = self.assoc_bundle["model"].predict(X)[0]
        return {
            "Association_Score": float(10 ** log_assoc),
            "log10_Association": float(log_assoc)
        }

    def predict_phase(self, smiles: str) -> int:
        X = self.smiles_to_fp(smiles)
        return int(self.phase_bundle["model"].predict(X)[0])

    def predict_target(self, smiles: str) -> str:
        X = self.smiles_to_fp(smiles)
        idx = int(self.target_bundle["model"].predict(X)[0])
        return self.target_bundle["label_encoder"].inverse_transform([idx])[0]

    def predict_all(self, smiles: str) -> Dict[str, Any]:
        X = self.smiles_to_fp(smiles)

        log_ic50 = self.ic50_bundle["model"].predict(X)[0]
        log_assoc = self.assoc_bundle["model"].predict(X)[0]
        phase = self.phase_bundle["model"].predict(X)[0]
        target_idx = self.target_bundle["model"].predict(X)[0]

        target = self.target_bundle["label_encoder"].inverse_transform(
            [int(target_idx)]
        )[0]

        return {
            "IC50": float(10 ** log_ic50),
            "Association_Score": float(10 ** log_assoc),
            "Max_Clinical_Phase": int(phase),
            "Predicted_Target": str(target)
        }

print(MatrixPredictor().predict_all("CC(=O)NC1=CC=C(O)C=C1"))