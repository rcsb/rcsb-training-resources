"""
Retrieve non-redundant set of for proteins that interact with ligands of interest.

General workflow:
    1. First, use the search API to get a list of chemical components of interest
    2. Then, use the search API to find PDB sequences that have components of interest as ligands in proximity
    3. Next, remove sequences that share at least 95% identity and return sequences from structures with the
       highest resolution

To run:
    python3 non_redundant_targets.py

"""


import requests
import time


def exec_search(query):
    base_url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    response = requests.post(base_url, json=query)

    if response.status_code == 200:
        return response.json()['result_set']
    else:
        print(f"Request failed with status code: {response.status_code}")


def search_ccd():

    ccd_search_query = {
        "query": {
            "type": "terminal",
            "label": "text_chem",
            "service": "text_chem",
            "parameters": {
                "attribute": "chem_comp.formula_weight",
                "operator": "greater",
                "value": 180
            }
        },
        "return_type": "mol_definition",
        "request_options": {
            "results_verbosity": "compact",
            "return_all_hits": True
        }
    }

    return exec_search(ccd_search_query)


def search_sequence_targets(ligands):

    targets_search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "greater",
                        "value": 20
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_ligand_neighbors.ligand_comp_id",
                        "operator": "in",
                        "value": ligands
                    }
                }
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "results_verbosity": "compact",
            "return_all_hits": True,
            "group_by": {
                "aggregation_method": "sequence_identity",
                "similarity_cutoff": 95,
                "ranking_criteria_type": {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc"
                }
            },
            "group_by_return_type": "representatives"
        }
    }

    return exec_search(targets_search_query)


if __name__ == "__main__":
    s_time = time.time()
    ccd_ids = search_ccd()
    print(f"Total number of CCD IDs: {len(ccd_ids)}")
    final_sequences = search_sequence_targets(ccd_ids)
    print(f"Total number of results: {len(final_sequences)}")
    print(f"List of results: {final_sequences}")
    e_time = time.time()
    print(f"Total execution time: {e_time - s_time:.4f} seconds")
