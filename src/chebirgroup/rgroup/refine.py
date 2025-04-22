import argparse
import json
import multiprocessing
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from chebirgroup.pubchem.compound import Compound
from chebirgroup.utils.multiprocessing import chunk_iterable
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem.rdRGroupDecomposition import RGroupDecompositionParameters
from rdkit import RDLogger
from sqlalchemy import and_, create_engine
from sqlalchemy.orm import sessionmaker
from tqdm import tqdm


RDLogger.DisableLog("rdApp.*")


def check_rgroup_base(data):
    ps = RGroupDecompositionParameters()
    ps.allowNonTerminalRGroups = True

    res = None
    mol = Chem.MolFromSmiles(data["smiles"])

    if mol:
        matched, unmatched = rdRGD.RGroupDecompose(
            data["targets"],
            [mol],
            asSmiles=True,
            asRows=True,
            options=ps,
        )
        if len(unmatched) < 1:
            res = dict(cid=data["cid"], smiles=data["smiles"])
    return res


def check_rgroup_only(data):
    ps = RGroupDecompositionParameters()
    ps.allowNonTerminalRGroups = True
    ps.onlyMatchAtRGroups = True

    res = None
    mol = Chem.MolFromSmiles(data["smiles"])

    if mol:
        matched, unmatched = rdRGD.RGroupDecompose(
            data["targets"],
            [mol],
            asSmiles=True,
            asRows=True,
            options=ps,
        )
        if len(unmatched) < 1:
            res = dict(cid=data["cid"], smiles=data["smiles"])
    return res


def fragment_mol(smiles, cids, original_wt):
    res = []
    for smiles_query in smiles.split("."):
        fragment = Chem.MolFromSmiles(smiles_query)
        if "*" in smiles_query:
            continue
        try:
            wt_target = Descriptors.ExactMolWt(fragment)
        except:
            continue
        if wt_target < original_wt:
            continue
        res.append((Chem.MolToSmiles(fragment), cids))
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-chebi-csv", required=True, help="Rhea csv file")
    parser.add_argument(
        "--parameter-chebi-int", required=True, type=int, help="ChEBI id int"
    )
    parser.add_argument(
        "--parameter-chunk-size-int", default=int(5e5), type=int, help="Manage chunk size for memory purpose"
    )
    parser.add_argument("--input-search-txt", required=True, help="Input txt file")
    parser.add_argument("--output-refine-json", required=True, help="Output json file")
    parser.add_argument("-t", "--input-thread-int", type=int, default=1, help="Threads")
    args = parser.parse_args()

    # Init
    threads = args.input_thread_int
    chunk_size = args.parameter_chunk_size_int
    finput = args.input_search_txt
    foutput = args.output_refine_json

    results = {}

    # Get Mol
    df = pd.read_csv(args.input_chebi_csv)
    original_smiles = None
    for _, row in df.iterrows():
        if f"CHEBI:{args.parameter_chebi_int}" in row["chebi"]:
            original_smiles = row["smiles_rhea"]
            break
    assert (
        original_smiles is not None
    ), f"Smiles not retrieved for ChEBI: {args.parameter_chebi_int}"

    original_mol = Chem.MolFromSmiles(original_smiles)
    original_wt = Descriptors.ExactMolWt(original_mol)

    smiles_sample = {}
    with open(finput) as fd:
        for line in fd:
            line = line.replace("\n", "")
            if line == "" or line == "|":
                continue
            smiles, cids = line.split("|")
            smiles_sample[smiles] = [cids]

    smiles_to_process = {}
    smiles_not_process = {}
    for smiles, cids in tqdm(smiles_sample.items(), total=len(smiles_sample)):
        if "." in smiles:
            smiles_to_process[smiles] = cids
            continue
        smiles_not_process[smiles] = cids

    # Check rgroup_base
    batchs = []
    for smiles, cids in smiles_to_process.items():
        batchs.append((smiles, cids))
    rgroup_bases = []
    for batch in chunk_iterable(batchs, chunk_size=chunk_size):
        with multiprocessing.Pool(processes=threads) as pool:
            async_results = [
                pool.apply_async(fragment_mol, (query[0], query[1], original_wt))
                for query in batch
            ]
            for async_result in async_results:
                try:
                    rgroup_bases.extend(async_result.get(timeout=300))
                except multiprocessing.TimeoutError:
                    continue
    del smiles_to_process

    new_smiles = []
    for batch in chunk_iterable(rgroup_bases, chunk_size=chunk_size):
        with multiprocessing.Pool(processes=threads) as pool:
            async_results = [
                pool.apply_async(
                    check_rgroup_base,
                    (
                        dict(
                            targets=[Chem.MolFromSmiles(original_smiles)],
                            smiles=query[0],
                            cid=query[1],
                        ),
                    ),
                )
                for query in batch
            ]
            for async_result in async_results:
                try:
                    new_smiles.append(async_result.get(timeout=10))
                except Exception:
                    continue
    del batchs

    results["base"] = {}
    for smiles, cids in smiles_not_process.items():
        results["base"][smiles] = cids
    for data in new_smiles:
        if data:
            if data["smiles"] in results["base"]:
                results["base"][data["smiles"]].extend(data["cid"])
            else:
                results["base"][data["smiles"]] = data["cid"]

    # Check rgroup_only
    batchs = []
    for smiles, cids in results["base"].items():
        batchs.append((smiles, cids))
    new_smiles = []
    for batch in chunk_iterable(batchs, chunk_size=chunk_size):
        with multiprocessing.Pool(processes=threads) as pool:
            async_results = [
                pool.apply_async(
                    check_rgroup_only,
                    (
                        dict(
                            targets=[Chem.MolFromSmiles(original_smiles)],
                            smiles=query[0],
                            cid=query[1],
                        ),
                    ),
                )
                for query in batch
            ]
            for async_result in async_results:
                try:
                    new_smiles.append(async_result.get(timeout=10))
                except Exception:
                    continue
    results["onlyrgroup"] = {}
    for data in new_smiles:
        if data:
            if data["smiles"] in results["onlyrgroup"]:
                results["onlyrgroup"][data["smiles"]].extend(data["cid"])
            else:
                results["onlyrgroup"][data["smiles"]] = data["cid"]

    # Write output
    with open(foutput, "w") as fd:
        json.dump(results, fd)
