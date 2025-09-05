"""
This script demos searching the PDB archive for the best ligand structure resolved 
by MX method (macromolecular crystallography, predominantly by X-ray) in terms of 
goodness-of-fit between the ligand model coordinates and the experimental data.

This script makes use of the `rcsb-api` Python package for performing this task.

The process includes 3 steps:
Step 1: Search for CCD IDs by name, then choose the first matched CCD ID
Step 2: Search for PDB IDs with the CCD ID from step 1
Step 3: Retrieve ligand quality metrics for each PDB entry, returns the PDB ID
with the best fitted ligand

Requirements:
    pip install "rcsb-api>=1.4.0"

Usage:
    python find_best_ligand_with_rcsbapi.py
Output:
    PDB ID with best ligand displayed on the terminal

"""

from rcsbapi.data import DataQuery
from rcsbapi.search import AttributeQuery, NestedAttributeQuery


## Step 1: Search for CCD IDs by name, then choose the first matched CCD ID

# provide example query payload for CCD ID search
ccd_search_query = AttributeQuery("chem_comp.name", operator="exact_match", value="IBUPROFEN", service="text_chem")

# execute query and retrieve CCD IDs from response
ccd_search_results = list(ccd_search_query.exec(return_type="mol_definition"))
print(f"found {len(ccd_search_results)} matching CCD IDs")

# use the first matched CCD ID for the subsequent PDB query
ccd_id = ccd_search_results[0]
print(f"searching PDB entries with {ccd_id}")


## Step 2: Search for PDB IDs with the CCD ID from step 1

# provide example query payload for PDB ID search
# ccd_id = "IBP"
sub_q1 = AttributeQuery("rcsb_nonpolymer_instance_annotation.comp_id", operator="exact_match", value=ccd_id, service="text")
sub_q2 = AttributeQuery("rcsb_nonpolymer_instance_annotation.type", operator="exact_match", value="HAS_NO_COVALENT_LINKAGE", service="text")
pdb_search_query = NestedAttributeQuery(sub_q1, sub_q2)

# execute query and retrieve the PDB IDs from the response
pdb_search_results = list(pdb_search_query.exec())
print(f"found {len(pdb_search_results)} PDB entries with {ccd_id}")


## Step 3: Retrieve ligand quality metrics for each PDB entry, returns the PDB ID with the best fitted ligand

# prepare data API query
data_query = DataQuery(
    input_type="entries",
    input_ids=pdb_search_results,
    return_data_list=[
        "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score",
        "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_annotation.comp_id",
    ]
)

# execute the query and retrieve the requested ligand quality metrics data
data_results = data_query.exec(progress_bar=True)
entry_data = data_results["data"]["entries"]

# review each PDB entry for best fitted ligand
pdb_id_best = None
best_score = 0

for entry in entry_data:
    pdb_id = entry["rcsb_id"]
    if "nonpolymer_entities" in entry:
        for each_entity in entry["nonpolymer_entities"]:
            for each_instance in each_entity["nonpolymer_entity_instances"]:
                ligand_id = each_instance["rcsb_nonpolymer_instance_annotation"][0]["comp_id"]
                if ligand_id.upper() == ccd_id.upper():
                    scores = each_instance["rcsb_nonpolymer_instance_validation_score"]
                    if scores:
                        score_1 = scores[0]
                        ranking_model_fit = score_1["ranking_model_fit"]
                        if ranking_model_fit and ranking_model_fit > best_score:
                            best_score = ranking_model_fit
                            pdb_id_best = pdb_id


# print result at the terminal
print(f"PDB entry {pdb_id_best} has the best fitted MX ligand structure for {ccd_id} at {best_score*100}%")
