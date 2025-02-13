"""  
Python script for fetching the chemical descriptors (SMILES, InChI, InChIKey) for all released chemical components (via RCSB PDB APIs).  
The results are saved into a CSV file, `rcsb_chemical_descriptors.csv`.  

This script requires the following packages, which can be installed with:   
pip install requests  
pip install rcsb-api  
""" 

import requests
import csv
from rcsbapi.data import DataQuery as Query

# Step 1: Retrieve all CCD IDs from Data API
url = 'https://data.rcsb.org/rest/v1/holdings/current/ccd_ids'
response = requests.get(url)
ids = eval(response.text)

# Step 2: Split full list of IDs into batches
batchSize = 5_000
idBatches = [ids[i:i+batchSize] for i in range(0, len(ids), batchSize)]

#Step 3: Query chemical descriptors
chemical_descriptors = []
spec = [
    ["SMILES (OpenEye)", "SMILES_CANONICAL", "OpenEye OEToolkits"],
    ["SMILES (OpenEye with stereo)", "SMILES", "OpenEye OEToolkits"],
    ["SMILES (CACTVS)", "SMILES_CANONICAL", "CACTVS"],
    ["SMILES (CACTVS with stereo)", "SMILES", "CACTVS"],
    ["InChI", "InChI", "InChI"],
    ["InChIKey", "InChIKey", "InChI"]
]
for batch in idBatches:
    query = Query(
        input_type="chem_comps",
        input_ids=batch,
        return_data_list=[
            "pdbx_chem_comp_descriptor.type", 
            "pdbx_chem_comp_descriptor.descriptor", 
            "pdbx_chem_comp_descriptor.program"
            ]
    )
    data = query.exec()
    for d in data['data']['chem_comps']:
        if (d["pdbx_chem_comp_descriptor"]):
            comp_id = d['rcsb_id']
            row = {"CCD ID": comp_id}
            descriptors = d["pdbx_chem_comp_descriptor"]
            for s in spec:
                l = [i for i in descriptors if i["type"]==s[1] and i["program"]==s[2]]
                descriptor = None
                if len(l) > 0:
                    descriptor = l[0]["descriptor"]
                row[s[0]]= descriptor
            chemical_descriptors.append(row)

with open("rcsb_chemical_descriptors.csv", "w") as handle:
    headers = list(chemical_descriptors[-1].keys())
    writer = csv.DictWriter(handle, fieldnames=headers)
    writer.writeheader()
    writer.writerows(chemical_descriptors)
print("Wrote chemical descriptors to rcsb_chemical_descriptors.csv")