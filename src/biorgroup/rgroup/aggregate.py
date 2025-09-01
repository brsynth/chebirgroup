import argparse
import ast
import json
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from biorgroup.utils.atomic import is_candidate
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-chebi-csv", required=True, help="Rhea csv file")
    parser.add_argument(
        "--input-refine-json", nargs="*", required=True, help="Input refine json files"
    )
    parser.add_argument(
        "--output-chebi-csv", required=True, help="Output chebi csv file"
    )
    args = parser.parse_args()

    # Init
    df = pd.read_csv(args.input_chebi_csv)
    df["chebi"] = df["chebi"].apply(ast.literal_eval)

    df["core_superstructure_smiles"] = [[] for _ in range(df.shape[0])]
    df["core_superstructure_pubchem_cid"] = [[] for _ in range(df.shape[0])]
    df["rgroup_extended_smiles"] = [[] for _ in range(df.shape[0])]
    df["rgroup_extended_pubchem_cid"] = [[] for _ in range(df.shape[0])]

    for ix, row in df.iterrows():
        if not is_candidate(row["smiles"]):
            continue
        # print("Deal with:", row["chebi"])

        chebi_id = row["chebi"][0]
        chebi_id = chebi_id.lower().replace(":", "_")

        path_refines = [
            x
            for x in args.input_refine_json
            if os.path.basename(x) == chebi_id + ".json"
        ]
        assert len(path_refines) == 1, "Retrieving refine file"

        with open(path_refines[0]) as fd:
            data_refine = json.load(fd)
        smiles, cids = [], []
        for smi, cid in data_refine["base"].items():
            smiles.append(smi)
            cids.append(cid)
        df.at[ix, "core_superstructure_smiles"] = smiles
        df.at[ix, "core_superstructure_pubchem_cid"] = cids

        smiles, cids = [], []
        for smi, cid in data_refine["onlyrgroup"].items():
            smiles.append(smi)
            cids.append(cid)
        df.at[ix, "rgroup_extended_smiles"] = smiles
        df.at[ix, "rgroup_extended_pubchem_cid"] = cids

    # Re-order columns
    cols = [
        "smiles",
        "chebi",
        "num_heavy_atoms",
        "exact_mol_wt",
        "core_superstructure_smiles",
        "core_superstructure_pubchem_cid",
        "rgroup_extended_smiles",
        "rgroup_extended_pubchem_cid",
    ]
    df = df[cols]

    # Save
    df.to_csv(args.output_chebi_csv, index=False)
