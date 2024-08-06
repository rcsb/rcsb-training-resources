# Any code that is provided should include an explicit list of all libraries that must be installed/imported
# along with the command needed to do so. In this case
# !pip install rcsbsearchapi
# !pip install python_graphql_client

# Minimal Python script to fetch all chemical component IDs and then fetch the SMILES, InChI, etc. strings associated with them
#
# Requires the following modules to be installed:
import json
from rcsbsearchapi.search import AttributeQuery
from python_graphql_client import GraphqlClient

# Step 1: Run search to retrieve all chemical component IDs
q1 = AttributeQuery("rcsb_chem_comp_container_identifiers.comp_id", "exists", service="text_chem")
results = [mol for mol in q1("mol_definition")]

subListSize = 300
resultsSubListL = [results[i:i+subListSize] for i in range(0, len(results), subListSize)]

# Step 2: Run data API query to retrieve data on all chemical components
url_data_api = 'https://data.rcsb.org/graphql'
client = GraphqlClient(endpoint=url_data_api)  # instantiate client with the RCSB Data API endpoint
query_method = """
query structure ($comp_ids: [String!]!) {
  chem_comps(comp_ids:$comp_ids){
        chem_comp {
            id
            name
            formula
            pdbx_formal_charge
            formula_weight
            type
        }
        rcsb_chem_comp_descriptor {
            InChI
            InChIKey
            SMILES
            SMILES_stereo
        }
    }
}
"""

# Iterate over each sublist and perform the data API query
for subList in resultsSubListL:
    query_variables = {"comp_ids": subList}
    dataResult = client.execute(query=query_method, variables=query_variables)
    data = dataResult['data']  # This will contain your data API query results--process/rewrangle this as needed or desired
    # Print out the first result
    print(json.dumps(data["chem_comps"][0], indent=2))  # Probably don't want to print out everything, but this is one way to do it while testing
