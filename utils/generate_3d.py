from rdkit import Chem
from rdkit.Chem import AllChem
import os

class Molecule3DGenerator:
    def __init__(self):
        self.output_dir = "3d_models"
        os.makedirs(self.output_dir, exist_ok=True)

    def _generate_3d(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        mol = Chem.AddHs(mol)

        status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if status != 0:
            raise RuntimeError("3D embedding failed")

        AllChem.UFFOptimizeMolecule(mol)

        output_path = os.path.join(self.output_dir, "molecule.sdf")
        writer = Chem.SDWriter(output_path)
        writer.write(mol)
        writer.close()

        return output_path
