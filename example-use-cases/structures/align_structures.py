"""
This script is designed to simplify the use of the RCSB.org Alignment API and make it more accessible (https://alignment.rcsb.org).

This service allows for the alignment of two individual structural chains and computes the RMSD and TM score (when available)
using one of multiple alignment algorithm options. Support for aligning multimeric assemblies is not yet available.

For help and a list of all possible options/flags, run:

    python alignment_api.py --help

EXAMPLE USAGE:

1) Align chain A of PDB ID 5QU3 with chain A of PDB ID 10GS:

    python alignment_api.py 5QU3 A 10GS A

2) Align a local mmCIF structure file (chain B) with PDB ID 10GS (chain A):

    python alignment_api.py local_file.cif B 10GS A

3) Specify the alignment method to use (by default, it uses fatcat-rigid):

    python alignment_api.py 5QU3 A 10GS A --method ce

4) If you select a method that has parameters, you can specify them as a list of parameter=value pairs:

    python alignment_api.py 5QU3 A 10GS A --method ce --params gap_max_size=20

5) If you want the result to be saved as a CSV file:

    python alignment_api.py 5QU3 A 10GS A --csv results.csv

"""

import requests
import json
import time
import argparse
import os
import csv

SUBMIT_URL = "https://alignment.rcsb.org/api/v1/structures/submit"

METHOD_DEFAULTS = {
    "fatcat-rigid": {
        "rmsd_cutoff": 3,
        "afp_dist_cutoff": 1600,
        "fragment_length": 8
    },

    "fatcat-flexible": {
        "rmsd_cutoff": 3,
        "afp_dist_cutoff": 1600,
        "fragment_length": 8,
        "max_num_twists": 5
    },

    "ce": {
        "gap_max_size": 30,
        "gap_opening_penalty": 5,
        "gap_extension_penalty": 0.5,
        "fragment_size": 8,
        "rmsd_threshold": 3,
        "max_opt_rmsd": 99
    },

    "ce-cp": {
        "gap_max_size": 30,
        "gap_opening_penalty": 5,
        "gap_extension_penalty": 0.5,
        "fragment_size": 8,
        "rmsd_threshold": 3,
        "max_opt_rmsd": 99,
        "min_cp_length": 5
    },

    "smith-waterman-3d": {
        "gap_opening_penalty": 8,
        "gap_extension_penalty": 1
    },

    "tm-align": {} 

}


def resolve_input(user_input):
    if os.path.isfile(user_input):
        return {"type": "file", "value": user_input}
    else:
        return {"type": "pdb_id", "value": user_input}


def parse_arguments():
    p = argparse.ArgumentParser(description="Run RCSB pairwise alignment")
    p.add_argument("input1", help="PDB ID or path to local PDB file")
    p.add_argument("chain1", help="Chian ID for first structure")
    p.add_argument("input2", help="PDB ID or path to local PDB file")
    p.add_argument("chain2", help="Chain ID for second structure")
    p.add_argument("--csv", help="Optional output CSV file to save summary results")

    p.add_argument(
        "--method",
        choices=[
            "tm-align",
            "fatcat-rigid",
            "fatcat-flexible",
            "ce",
            "ce-cp",
            "smith-waterman-3d"
        ],
        default="fatcat-rigid",
        help="Alignment method (default: fatcat-rigid)"
    )

    p.add_argument(
        "--params",
        nargs="*",
        help="Optional method parameters as key=value pairs"
    )

    return p.parse_args()


def build_method_parameters(method, user_params):
    params = METHOD_DEFAULTS.get(method, {}).copy()

    if not user_params:
        return params

    for item in user_params:
        if "=" not in item:
            raise ValueError(f"Invalid parameter format: {item}. Use key=value.")

        key, value = item.split("=", 1)

        if key not in params:
            raise ValueError(f"Parameter '{key}' not valid for method '{method}'")

        try:
            if "." in value:
                value = float(value)
            else:
                value = int(value)
        except ValueError:
            pass

        params[key] = value

    return params


def build_structure_json(resolved, chain):
    if resolved["type"] == "pdb_id":
        return {
            "entry_id": resolved["value"],
            "selection": {"asym_id": chain}
        }

    if resolved["type"] == "file":
        return {
            "format": "mmcif",
            "selection": {"asym_id": chain}
        }


def submit_alignment_job(pdb1, pdb2, chain1, chain2, method, user_params):

    structure1 = build_structure_json(pdb1, chain1)
    structure2 = build_structure_json(pdb2, chain2)
    method_parameters = build_method_parameters(method, user_params)

    query = {
        "context": {
            "mode": "pairwise",
            "method": {"name": method,
                       "parameters": method_parameters},
            "structures": [structure1, structure2]
        }        
    }

    data = {"query": json.dumps(query)}

    files = []


    if pdb1["type"] == "file":
        files.append(
            (f"files", ("filename", open(pdb1["value"], "r"))),
        )

    if pdb2["type"] == "file":
        files.append(
            (f"files", ("filename", open(pdb2["value"], "r"))),
        )

    response = requests.post(SUBMIT_URL, data=data, files=files)

    if response.status_code != 200:
        raise Exception(f"Job submission failed with {response.status_code}: {response.text}")

    ticket = response.text.strip()
    print(f"Submitted job → Ticket: {ticket}")

    return ticket


def get_alignment_results(ticket):

    RESULTS_URL = f"https://alignment.rcsb.org/api/v1/structures/results?uuid={ticket}"

    while True:
        response = requests.get(RESULTS_URL)

        if response.status_code != 200:
            print("Waiting for job to finish...")
            time.sleep(2)
            continue

        try:
            js = response.json()
        except:
            print("Not ready yet...")
            time.sleep(1)
            continue

        status = None

        if "info" in js:
            status = js["info"].get("status", "").upper()

        print(f"Job status: {status}")

        if status == "COMPLETE":
            return js

        time.sleep(1)
        print("Results:")
        print(response.text)


def extract_summary(results, method):
    print(json.dumps(results, indent=2))

    try:
        if isinstance(results.get("results"), list) and len(results["results"]) > 0:
            result_obj = results["results"][0]
        else:
            print("Results array is missing or empty.")
            return

        summary = result_obj.get("summary", {})
        scores = summary.get("scores", [])

        rmsd = "N/A"
        tmscore = "N/A"

        for s in scores:
            if not isinstance(s, dict):
                continue
            t = s.get("type", "").lower()
            if t == "rmsd":
                rmsd = s.get("value", "N/A")
            elif t in ("tm-score", "tmscore", "tm_score"):
                tmscore = s.get("value", "N/A")

        print("\n===== SUMMARY =====")
        print(f"You used {method} for alignment")
        print(f"RMSD: {rmsd}")
        print(f"TM-score: {tmscore}")
        print("===================")

    except Exception as e:
        print("Could not extract summary:", e)


def write_csv_summary(csv_path, ticket, method, results):
    try:
        if isinstance(results.get("results"), list) and len(results["results"]) > 0:
            result_obj = results["results"][0]
        else:
            print("No results to write to CSV.")
            return

        summary = result_obj.get("summary", {})
        scores = summary.get("scores", [])

        rmsd = "N/A"
        tmscore = "N/A"

        for s in scores:
            if not isinstance(s, dict):
                continue
            t = s.get("type", "").lower()
            if t == "rmsd":
                rmsd = s.get("value", "N/A")
            elif t in ("tm-score", "tmscore", "tm_score"):
                tmscore = s.get("value", "N/A")

        file_exists = os.path.isfile(csv_path)

        with open(csv_path, mode="a", newline="") as f:
            writer = csv.writer(f)

            if not file_exists:
                writer.writerow(["ticket", "method", "rmsd", "tm_score"])

            writer.writerow([ticket, method, rmsd, tmscore])

        print(f"CSV written to {csv_path}")

    except Exception as e:
        print("Failed to write CSV:", e)


def main():

    args = parse_arguments()

    pdb1 = resolve_input(args.input1)
    pdb2 = resolve_input(args.input2)

    chain1 = args.chain1
    chain2 = args.chain2

    ticket = submit_alignment_job(pdb1, pdb2, chain1, chain2, args.method, args.params) 
    results = get_alignment_results(ticket)

    extract_summary(results, args.method)

    if args.csv:
        write_csv_summary(args.csv, ticket, args.method, results)


if __name__ == "__main__":
    main()
