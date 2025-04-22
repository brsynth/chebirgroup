import argparse
import ast
import json
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from chebirgroup.utils.atomic import is_candidate
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-chebi-csv", required=True, help="Rhea csv file")
    parser.add_argument(
        "--input-refine-json", nargs="*", required=True, help="Input refine json files"
    )
    parser.add_argument("--output-chebi-csv", required=True, help="Output chebi csv file")
    args = parser.parse_args()

    # Init
    df = pd.read_csv(args.input_chebi_csv)
    df["chebi"] = df["chebi"].apply(ast.literal_eval)

    df["additional_substituents_smiles"] = [[] for _ in range(df.shape[0])]
    df["additional_substituents_pubchem_cid"] = [[] for _ in range(df.shape[0])]
    df["only_match_rgroup_smiles"] = [[] for _ in range(df.shape[0])]
    df["only_match_rgroup_pubchem_cid"] = [[] for _ in range(df.shape[0])]

    for ix, row in df.iterrows():
        if not is_candidate(row["smiles_rhea"]):
            continue
        print("Deal with:", row["chebi"])

        chebi_id = row["chebi"][0]
        chebi_id = chebi_id.lower().replace(":", "_")

        path_refines = [x for x in args.input_refine_json if os.path.basename(x) == chebi_id + ".json"]
        assert len(path_refines) == 1, "Retrieving refine file"

        with open(path_refines[0]) as fd:
            data_refine = json.load(fd)
        data = defaultdict(list)
        for smi, cids in data_refine["base"].items():
            cid = ",".join(cids)
            data[smi].append(cid)
        smiles, cids = [], []
        for smi, cid in data.items():
            smiles.append(smi)
            cids.append((",".join(cid)).split(","))
        df.at[ix, "additional_substituents_smiles"] = smiles
        df.at[ix, "additional_substituents_pubchem_cid"] = cids

        data = defaultdict(list)
        for smi, cids in data_refine["onlyrgroup"].items():
            cid = ",".join(cids)
            data[smi].append(cid)
        smiles, cids = [], []
        for smi, cid in data.items():
            smiles.append(smi)
            cids.append((",".join(cid)).split(","))
        df.at[ix, "only_match_rgroup_smiles"] = smiles
        df.at[ix, "only_match_rgroup_pubchem_cid"] = cids

    df.to_csv(args.output_chebi_csv, index=False)
    print("End")
