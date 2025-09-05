"""
This script demos searching the PDB archive for the best ligand structure resolved
by MX method (macromolecular crystallography, predominantly by X-ray) in terms of
goodness-of-fit between the ligand model coordinates and the experimental data.

This script allows use of "copy-and-paste" from RCSB.org API query editors for performing this task.

The process includes 3 steps:
Step 1: Search for CCD IDs by name, then choose the first matched CCD ID
Step 2: Search for PDB IDs with the CCD ID from step 1
Step 3: Retrieve ligand quality metrics for each PDB entry, returns the PDB ID
with the best fitted ligand

Requirements:
    pip install requests
    pip install python_graphql_client

Usage:
    python find_best_ligand.py
Output:
    PDB ID with best ligand displayed on the terminal

"""

import json
import requests
import sys
from python_graphql_client import GraphqlClient

## Step 1: Search for CCD IDs by name, then choose the first matched CCD ID

# define search API end point
search_api_url = "https://search.rcsb.org/rcsbsearch/v2/query"

# provide example query payload for CCD ID search
attribute = "chem_comp.name"
value = "IBUPROFEN"
payload_ccd_search = {
    "query": {
        "type": "terminal",
        "label": "text_chem",
        "service": "text_chem",
        "parameters": {
            "attribute": attribute,
            "operator": "exact_match",
            "value": value,
        },
    },
    "return_type": "mol_definition",
    "request_options": {
        "paginate": {"start": 0, "rows": 25},
        "results_content_type": ["experimental"],
        "sort": [{"sort_by": "score", "direction": "desc"}],
        "scoring_strategy": "combined",
    },
}

# send POST request
response_ccd = requests.post(search_api_url, json=payload_ccd_search)

# check response
if response_ccd.status_code == 200:
    results_ccd = response_ccd.json()
else:
    print(f"Error {response_ccd.status_code}: {response_ccd.text}")
    sys.exit()

# retrieve CCD IDs from response
l_ccd_id = []
for each in results_ccd["result_set"]:
    l_ccd_id.append(each["identifier"])
print(f"found {len(l_ccd_id)} CCD IDs with example match to {attribute} : {value}")

# use the first matched CCD ID for the subsequent PDB query
ccd_id = l_ccd_id[0]
print(f"searching PDB entries with {ccd_id}")


## Step 2: Search for PDB IDs with the CCD ID from step 1

# provide example query payload for PDB ID search
# ccd_id = "IBP"
payload_pdb_search = {
    "query": {
        "type": "group",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                    "operator": "exact_match",
                    "value": ccd_id,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_instance_annotation.type",
                    "operator": "exact_match",
                    "value": "HAS_NO_COVALENT_LINKAGE",
                },
            },
        ],
        "logical_operator": "and",
        "label": "nested-attribute",
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {"start": 0, "rows": 25},
        "results_content_type": ["experimental"],
        "sort": [{"sort_by": "score", "direction": "desc"}],
        "scoring_strategy": "combined",
    },
}

# send POST request
response_pdb = requests.post(search_api_url, json=payload_pdb_search)

# check response
if response_pdb.status_code == 200:
    results_pdb = response_pdb.json()
else:
    print(f"Error {response_pdb.status_code}: {response_pdb.text}")
    sys.exit()

# retrieve PDB IDs from response
l_pdb_id = []
for each in results_pdb["result_set"]:
    l_pdb_id.append(each["identifier"])
print(f"found {len(l_pdb_id)} PDB entries with {ccd_id}")


## Step 3: Retrieve ligand quality metrics for each PDB entry, returns the PDB ID with the best fitted ligand

# define data API end point and initate query client
data_api_url = "https://data.rcsb.org/graphql"
client = GraphqlClient(endpoint = data_api_url)

# provide query to retrieve ligand quality metrics
# (copy and paste from RCSB.org Data API query editor)
ligand_query = """
query test($pdb_id: String!) {
  entry(entry_id:$pdb_id){
    nonpolymer_entities {
      nonpolymer_entity_instances {
        rcsb_nonpolymer_instance_validation_score {
          average_occupancy,
          mogul_bonds_RMSZ,
          mogul_angles_RMSZ,
          mogul_bond_outliers,
          mogul_angle_outliers,
          RSR, 
          RSCC,
          completeness,
          ranking_model_fit,
          ranking_model_geometry
        }
        rcsb_nonpolymer_instance_annotation {
          comp_id
        }
      }
    }
  }
}
"""

# review each PDB entry for best fitted ligand
pdb_id_best = None
best_score = 0
for pdb_id in l_pdb_id:
    print(f"checking ligand quality in {pdb_id}")
    ligand_query_var = {"pdb_id": pdb_id}
    ligand_data = client.execute(query=ligand_query, variables=ligand_query_var)
    if ligand_data["data"]["entry"]["nonpolymer_entities"]:
        for each_entity in ligand_data["data"]["entry"]["nonpolymer_entities"]:
            for each_instance in each_entity["nonpolymer_entity_instances"]:
                ligand_id = each_instance["rcsb_nonpolymer_instance_annotation"][0][
                    "comp_id"
                ]
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
