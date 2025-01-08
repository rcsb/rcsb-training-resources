import requests
import csv
from dateutil import parser
from rcsbapi.data import DataQuery as Query

# Step 1: Retrieve all PDB IDs from Data API
url = 'https://data.rcsb.org/rest/v1/holdings/current/entry_ids'
response = requests.get(url)
ids = eval(response.text)

# Step 2: Split full list of IDs into batches
batchSize = 5_000
idBatches = [ids[i:i+batchSize] for i in range(0, len(ids), batchSize)]

#Step 3: Query release date
release_dates = []
for batch in idBatches:
    query = Query(
        input_type="entries",
        input_ids=batch,
        return_data_list=["rcsb_accession_info.initial_release_date"]
    )
    data = query.exec()
    for d in data['data']['entries']:
        entry_id = d['rcsb_id']
        isodate = d["rcsb_accession_info"]["initial_release_date"]
        date = parser.parse(isodate).strftime('%Y-%m-%d')
        release_dates.append({
            "pdb_id": entry_id,
            "release_date": date
        })

with open("rcsb_release_dates.csv", "w") as handle:
    headers = list(release_dates[-1].keys())
    writer = csv.DictWriter(handle, fieldnames=headers)
    writer.writeheader()
    writer.writerows(release_dates)
print("Wrote release dates to rcsb_release_dates.csv")