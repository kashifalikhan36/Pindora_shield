from Tengan.generate_from_smiles import MoleculeGenerator
from fetch_data import FetchData
from copilot import AzureOpenAIChatClient
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem
from rdkit import DataStructs

data_processor = FetchData()
copilot=AzureOpenAIChatClient()
model_path = "Tengan/res/save_models/ZINC/TenGAN_0.5/rollout_8/batch_64/druglikeness/g_pretrained.pkl"
generator = MoleculeGenerator(
    model_path=model_path,
    batch_size=8,
    max_len=120
)

disease_c=copilot.generate_desease_name_from_prompt("i have breast cancer and diabetes type 2")
diseases=json.loads(disease_c)["desease"]
all_data = []

for disease_name in diseases:
    efo_ids = data_processor.map_disease_to_efo(disease_name)
    print(f"Disease: {disease_name} -> EFO IDs: {efo_ids}")
    if not efo_ids:
        print(f"No EFO ID found for disease: {disease_name}")
    for i in efo_ids:
        targets = data_processor.get_associated_targets(i, max_targets=50)
        all_data = []
        for i, target in enumerate(targets):
            drugs = data_processor.get_known_drugs_for_target(target["target_id"], max_drugs=10)
            if len(drugs)==0:
                continue
            for j, drug in enumerate(drugs):
                ic50_data = data_processor.get_ic50_data_for_molecule(drug["drug_id"], limit=100)
                features = data_processor.get_molecule_properties(drug["drug_id"])
                for ic50 in ic50_data:
                    row = {
                        "disease_name": disease_name,
                        "efo_id": i,
                        
                        # Target information
                        "target_id": target["target_id"],
                        "target_symbol": target["approved_symbol"],
                        "association_score": target["association_score"],
                        
                        # Drug information
                        "drug_id": drug["drug_id"],
                        "drug_name": drug["pref_name"],
                        "clinical_phase": drug["phase"],
                        
                        # IC50 bioactivity data
                        "ic50_value": ic50["standard_value"],
                        "ic50_units": ic50["standard_units"],
                        "target_chembl_id": ic50["target_chembl_id"],
                        "assay_chembl_id": ic50["assay_chembl_id"],
                        "pchembl_value": ic50["pchembl_value"]
                    }
                    for key, value in features.items():
                        if key != "molecule_chembl_id" and value is not None:
                            row[key] = value

                    all_data.append(row)
                    break
                break

print(f"Total records collected: {len(all_data)}")

def split_smiles_components(smiles: str) -> list:
    if not smiles or not isinstance(smiles, str):
        return []
    
    if '$' in smiles:
        components = smiles.split('$')
    else:
        components = [smiles]
    
    cleaned = [c.strip() for c in components if c.strip() and len(c.strip()) > 1]
    return cleaned if cleaned else [smiles.strip()]

def calculate_similarity(smiles1: str, smiles2: str) -> float:
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        if mol1 and mol2:
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
            return DataStructs.TanimotoSimilarity(fp1, fp2)
    except:
        pass
    return 0.0

def get_molecular_properties(smiles: str) -> dict:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                "valid": True,
                "molecular_weight": Descriptors.MolWt(mol),
                "logp": Crippen.MolLogP(mol),
                "hbd": Descriptors.NumHDonors(mol),
                "hba": Descriptors.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "aromatic_rings": Descriptors.NumAromaticRings(mol)
            }
    except:
        pass
    return {"valid": False}

gen_mol = []
for i in all_data:
    input_smiles = i["smiles"]
    results = generator.generate_from_smiles(input_smiles, 5)
    
    generated_with_props = []
    for result_smiles in results:
        components = split_smiles_components(result_smiles)
        component_data = []
        
        for component in components:
            props = get_molecular_properties(component)
            similarity = calculate_similarity(input_smiles, component)
            component_data.append({
                "smiles": component,
                "similarity": round(similarity, 3),
                "properties": props
            })
        
        generated_with_props.append(components)
    
    gen_mol.append({
        "input_smile": input_smiles,
        "disease_name": i["disease_name"],
        "target_symbol": i["target_symbol"],
        "drug_name": i["drug_name"],
        "generated_molecules": generated_with_props
    })
    
    print(f"\nInput: {input_smiles} (Drug: {i['drug_name']})")
    for j, components in enumerate(generated_with_props, 1):
        print(f"  Generated {j}: {components}")

with open("generated_molecules_new.json", "w", encoding="utf-8") as f:
    json.dump(gen_mol, f, indent=2)