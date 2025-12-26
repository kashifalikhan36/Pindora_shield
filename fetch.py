import requests
import json
import pandas as pd
from typing import List, Dict, Any, Optional
import time

# ============================
# API Endpoints
# ============================
OPEN_TARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"

def query_open_targets(query: str, variables: Dict[str, Any] = None, max_retries: int = 3) -> Dict[str, Any]:
    payload = {"query": query}
    if variables:
        payload["variables"] = variables
    for attempt in range(max_retries):
        try:
            response = requests.post(OPEN_TARGETS_URL, json=payload, timeout=30)
            response.raise_for_status()
            data = response.json()
            return data
            
        except requests.exceptions.HTTPError as e:
            if response.status_code in [502, 503, 504] and attempt < max_retries - 1:
                wait_time = 2 ** attempt
                time.sleep(wait_time)
                continue
            else:
                raise
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                time.sleep(wait_time)
                continue
            else:
                raise
    

def query_chembl(endpoint: str, params: Dict[str, Any] = None) -> Dict[str, Any]:
    if endpoint.endswith(".json"):
        url = f"{CHEMBL_URL}/{endpoint}"
    elif "/" in endpoint and "molecule/" in endpoint:
        url = f"{CHEMBL_URL}/{endpoint}.json"
    else:
        url = f"{CHEMBL_URL}/{endpoint}.json"
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()



# ============================
# Step 1: Disease → EFO ID
# ============================

def map_disease_to_efo(disease_name: str) -> List[str]:

    query = """
        query MapIds($terms: [String!]!, $entityNames: [String!]) {
            mapIds(queryTerms: $terms, entityNames: $entityNames) {
                mappings {
                    term
                    hits {
                        id
                        entity
                    }
                }
            }
        }
    """
    
    variables = {
        "terms": [disease_name],
        "entityNames": ["disease"]
    }
    
    data = query_open_targets(query, variables)
    efo_ids = []
    
    for mapping in data["data"]["mapIds"]["mappings"]:
        for hit in mapping["hits"]:
            if hit["entity"] == "disease":
                efo_ids.append(hit["id"])
    
    return efo_ids

# ============================
# Step 2: Disease → Target Proteins
# ============================

def get_associated_targets(efo_id: str, max_targets: int = 50) -> List[Dict[str, Any]]:
    query = """
        query AssociatedTargets($efoId: String!, $pageIndex: Int!, $pageSize: Int!) {
            disease(efoId: $efoId) {
                associatedTargets(page: { index: $pageIndex, size: $pageSize }) {
                    count
                    rows {
                        score
                        target {
                            id
                            approvedSymbol
                        }
                    }
                }
            }
        }
    """
    
    targets = []
    page_index = 0
    page_size = 50

    while len(targets) < max_targets:
        variables = {
            "efoId": efo_id,
            "pageIndex": page_index,
            "pageSize": page_size
        }
        
        data = query_open_targets(query, variables)
        rows = data["data"]["disease"]["associatedTargets"]["rows"]
        
        if not rows:
            break
        for row in rows:
            targets.append({
                "target_id": row["target"]["id"],
                "approved_symbol": row["target"]["approvedSymbol"],
                "association_score": row["score"]
            })
        
        page_index += 1
        
        if len(rows) < page_size:
            break
    return targets[:max_targets]

# ============================
# Step 3: Target → Known Drugs
# ============================

def get_known_drugs_for_target(target_id: str, max_drugs: int = 50) -> List[Dict[str, Any]]:

    query = """
        query KnownDrugs($targetId: String!, $size: Int!, $cursor: String) {
            target(ensemblId: $targetId) {
                knownDrugs(size: $size, cursor: $cursor) {
                    count
                    cursor
                    rows {
                        drugId
                        prefName
                        phase
                    }
                }
            }
        }
    """

    drugs = []
    cursor = None
    size = 50

    while len(drugs) < max_drugs:
        variables = {
            "targetId": target_id,
            "size": size,
            "cursor": cursor
        }

        data = query_open_targets(query, variables)
        known_drugs = data["data"]["target"]["knownDrugs"]

        for row in known_drugs["rows"]:
            drugs.append({
                "drug_id": row["drugId"],
                "pref_name": row["prefName"],
                "phase": row["phase"]
            })
        cursor = known_drugs.get("cursor")
        if not cursor:
            break
    return drugs[:max_drugs]

# ============================
# Step 4: Drug → IC50 Data (ChEMBL)
# ============================

def get_ic50_data_for_molecule(molecule_chembl_id: str, limit: int = 1000) -> List[Dict[str, Any]]:
    params = {
        "molecule_chembl_id": molecule_chembl_id,
        "standard_type__exact": "IC50",
        "limit": limit,
        "offset": 0
    }
    
    data = query_chembl("activity", params)

    ic50_records = []
    
    for activity in data.get("activities", []):
        if activity.get("standard_value") is None:
            continue
        
        ic50_records.append({
            "molecule_chembl_id": molecule_chembl_id,
            "standard_value": activity["standard_value"],
            "standard_units": activity.get("standard_units", "nM"),
            "target_chembl_id": activity.get("target_chembl_id"),
            "target_pref_name": activity.get("target_pref_name"),
            "assay_chembl_id": activity.get("assay_chembl_id"),
            "pchembl_value": activity.get("pchembl_value"),
            "document_chembl_id": activity.get("document_chembl_id")
        })
    
    return ic50_records

# ============================
# Step 5: Drug → Molecular Features (ChEMBL)
# ============================

def get_molecule_properties(molecule_chembl_id: str) -> Dict[str, Any]:
    data = query_chembl(f"molecule/{molecule_chembl_id}")
    molecule = data.get("molecule", {})
    
    props = molecule.get("molecule_properties", {})
    return {
        "molecule_chembl_id": molecule_chembl_id,
        "pref_name": molecule.get("pref_name"),
        "full_mwt": props.get("full_mwt"),
        "alogp": props.get("alogp"),
        "hba": props.get("hba"),
        "hbd": props.get("hbd"),
        "psa": props.get("psa"),
        "rtb": props.get("rtb"),
        "ro5_violations": props.get("num_ro5_violations"),
        "qed_weighted": props.get("qed_weighted")
    }


# ============================
# Main Pipeline
# ============================

def pipeline(disease_name: str, max_targets: int = 5, max_drugs_per_target: int = 5) -> pd.DataFrame:

    # ================================================================================
    # STEP 1: Map disease name to EFO identifier
    # ================================================================================
    efo_ids = map_disease_to_efo(disease_name)
    
    if not efo_ids:
        return pd.DataFrame()
    
    efo_id = efo_ids[0]
    # ================================================================================
    # STEP 2: Retrieve associated targets (genes/proteins) for the disease
    # ================================================================================
    targets = get_associated_targets(efo_id, max_targets=max_targets)
    # ================================================================================
    # STEP 3-5: For each target, get drugs, IC50 data, and molecular features
    # ================================================================================
    all_data = []

    for i, target in enumerate(targets, 1):
        # ----------------------------------------------------------------------------
        # STEP 3: Fetch known drugs for this target
        # ----------------------------------------------------------------------------
        drugs = get_known_drugs_for_target(target["target_id"], max_drugs=max_drugs_per_target)
        
        if not drugs:
            print(f"  ⚠ No drugs found for this target, skipping...")
            continue

        for j, drug in enumerate(drugs, 1):
            # ------------------------------------------------------------------------
            # STEP 4: Retrieve IC50 data from ChEMBL
            # ------------------------------------------------------------------------
            ic50_data = get_ic50_data_for_molecule(drug["drug_id"], limit=100)
            # ------------------------------------------------------------------------
            # STEP 5: Fetch molecular features/properties
            # ------------------------------------------------------------------------
            features = get_molecule_properties(drug["drug_id"])

            for ic50 in ic50_data:
                row = {
                    # Disease information
                    "disease_name": disease_name,
                    "efo_id": efo_id,
                    
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
                    if value is not None and key != "molecule_chembl_id":
                        row[key] = value
                
                all_data.append(row)
            
            if ic50_data:
                print(f"      ✓ Added {len(ic50_data)} rows to dataset")
            else:
                print(f"      ⚠ No IC50 data available for this drug")
    
    return pd.DataFrame(all_data)

# ============================
# Example Usage & Discussion
# ============================

def collect_large_dataset(diseases: List[str], max_targets: int = 50, max_drugs_per_target: int = 10, 
                         output_file: str = "ic50_large_dataset.csv") -> pd.DataFrame:
    all_dataframes = []
    total_records = 0
    
    for i, disease in enumerate(diseases, 1):
        
        try:
            df = pipeline(disease, max_targets=max_targets, max_drugs_per_target=max_drugs_per_target)
            
            if not df.empty:
                all_dataframes.append(df)
                total_records += len(df)
                if all_dataframes:
                    temp_df = pd.concat(all_dataframes, ignore_index=True)
                    temp_df.to_csv(output_file, index=False)
            else:
                print(f"\n⚠ No IC50 data collected for {disease}")
                
        except Exception as e:
            if all_dataframes:
                temp_df = pd.concat(all_dataframes, ignore_index=True)
                temp_df.to_csv(output_file, index=False)
            continue
    
    if all_dataframes:
        final_df = pd.concat(all_dataframes, ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()
    


def main():
    diseases = [
        "breast cancer",
        "prostate cancer",
        "lung cancer",
        "colorectal cancer",
        "diabetes type 2",
        "hypertension",
        "alzheimer disease"
    ]
    
    df = collect_large_dataset(diseases, max_targets=50, max_drugs_per_target=10, 
                              output_file="ic50_large_dataset.csv")

    if df.empty:
        return

    display_cols = ['disease_name', 'target_symbol', 'drug_name', 'ic50_value', 'ic50_units', 'pchembl_value']
    available_cols = [col for col in display_cols if col in df.columns]

    # Save to CSV
    output_file = "ic50_large_dataset.csv"
    df.to_csv(output_file, index=False)

    target_counts = df.groupby(['target_symbol', 'target_id']).size().sort_values(ascending=False).head(10)

    drug_counts = df.groupby(['drug_name', 'drug_id']).size().sort_values(ascending=False).head(10)
if __name__ == "__main__":
    main()