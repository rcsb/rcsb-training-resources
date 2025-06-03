import requests
from python_graphql_client import GraphqlClient
import json
from urllib.parse import quote
from rcsbapi.search import AttributeQuery
from rcsbapi.data import DataQuery as Query
import time


def exec_search_library():
    q1 = AttributeQuery(
        attribute="rcsb_polymer_instance_feature_summary.type",
        operator="exact_match",
        value="UNOBSERVED_RESIDUE_XYZ"
    )

    q2 = AttributeQuery(
        attribute="rcsb_polymer_instance_feature_summary.coverage",
        operator="greater",
        value=0
    )

    q3 = AttributeQuery(
        attribute="entity_poly.rcsb_entity_polymer_type",
        operator="exact_match",
        value="Protein"
    )

    q4 = AttributeQuery(
        attribute="entity_poly.rcsb_sample_sequence_length",
        operator="greater",
        value=20
    )

    query = (q1 & q2) & q3 & q4

    result_set = []
    i = 0
    for result in query(return_type="polymer_instance", results_verbosity="compact", rows=1000):
        result_set.append(result)
        i += 1
        if i == 1000:
            break
    # result_set = list(query(
    #     return_type="polymer_instance",
    #     results_verbosity="compact"
    # ))
    return result_set


def exec_search():
    search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_instance_feature_summary.coverage",
                                "operator": "greater",
                                "value": 0
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_instance_feature_summary.type",
                                "operator": "exact_match",
                                "value": "UNOBSERVED_RESIDUE_XYZ"
                            }
                        }
                    ]
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
                    "label": "text",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "greater",
                        "value": 20
                    }
                }
            ]
        },
        "return_type": "polymer_instance",
        "request_options": {
            # Every search hit is returned as a simple string, e.g. "4HHB.A", with no additional metadata
            "results_verbosity": "compact",
            "return_all_hits": True,
            # "paginate": {
            #     "start": 0,
            #     "rows": 1000
            # }
        }
    }

    json_string = json.dumps(search_query)
    encoded_query = quote(json_string)

    base_url = 'https://search.rcsb.org/rcsbsearch/v2/query'
    url = base_url + '?json=' + encoded_query
    response = requests.get(url)

    if response.status_code == 200:
        return response.json()['result_set']
    else:
        print(f"Request failed with status code: {response.status_code}")


def exec_data_library(id_batches):
    selected_chain_ids = []
    for batch in id_batches:
        data_query = Query(
            input_type="polymer_entity_instances",
            input_ids=batch,
            return_data_list=[
                "rcsb_id",
                "polymer_entity.entity_poly.rcsb_sample_sequence_length",
                "rcsb_polymer_instance_feature.type",
                "rcsb_polymer_instance_feature.feature_positions.beg_seq_id",
                "rcsb_polymer_instance_feature.feature_positions.end_seq_id"
            ]
        )
        data = data_query.exec()
        parse_data(data, selected_chain_ids)
    return selected_chain_ids


def exec_data(id_batches):
    base_url = 'https://data.rcsb.org/graphql'
    client = GraphqlClient(endpoint=base_url)
    query_method = """
    query chains_with_missing_coords ($ids: [String!]!) {
      polymer_entity_instances(instance_ids:$ids){
            rcsb_id
            polymer_entity {
              entity_poly {
                rcsb_sample_sequence_length
              }
            }
            rcsb_polymer_instance_feature {
              type
              feature_positions {
                beg_seq_id
                end_seq_id
              }
            }
        }
    }
    """

    selected_chain_ids = []
    for batch in id_batches:
        query_variables = {"ids": batch}
        # Execute the query with a timeout (in seconds)
        data = client.execute(query=query_method, variables=query_variables, timeout=60)
        parse_data(data, selected_chain_ids)
    return selected_chain_ids


def parse_data(data, ids):
    for d in data['data']['polymer_entity_instances']:
        pdb_id = d['rcsb_id']
        sequence_length = d['polymer_entity']['entity_poly']['rcsb_sample_sequence_length']
        for f in d['rcsb_polymer_instance_feature']:
            if f['type'] == 'UNOBSERVED_RESIDUE_XYZ':
                unobserved_range_count = 0
                for r in f['feature_positions']:
                    # include only non-terminal regions
                    if r['beg_seq_id'] != 1 and r['end_seq_id'] != sequence_length:
                        unobserved_range_count += 1
                # include only chains with 2 or more non-terminal regions
                if unobserved_range_count >= 2:
                    ids.append(pdb_id)


if __name__ == "__main__":
    s_time = time.time()

    # Find all protein chains with missing coordinates
    search_s_time = time.time()
    chain_ids = exec_search()
    # chain_ids = exec_search_library()
    search_e_time = time.time()
    print(f"Search execution time: {search_e_time - search_s_time:.4f} seconds. Retrieved {len(chain_ids)} IDs")

    # Fetch data for missing coordinates sequence ranges and select chains with non-terminal, non-contiguous regions
    data_s_time = time.time()
    size = 5_000
    batches = [chain_ids[i:i + size] for i in range(0, len(chain_ids), size)]
    final_chains = exec_data(batches)
    # final_chains = exec_data_library(batches)
    final_structures = list(set([s.split('.')[0] for s in final_chains]))
    data_e_time = time.time()
    print(f"Data execution time: {data_e_time - data_s_time:.4f} seconds")

    e_time = time.time()
    print(f"Total execution time: {e_time - s_time:.4f} seconds")
    if len(final_chains) > 0:
        print(f"Total number of chains with missing coordinates in non-terminal, non-contiguous regions: "
              f"{len(final_chains)}. Example: {final_chains[0]}")
        print(f"Total number of structures with missing coordinates in non-terminal, non-contiguous regions: "
              f"{len(final_structures)}")
    else:
        print("No chains matching criteria is found")
