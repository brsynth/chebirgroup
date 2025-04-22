import argparse
import multiprocessing
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from chebirgroup.pubchem.compound import Compound
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem.rdRGroupDecomposition import RGroupDecompositionParameters
from rdkit import RDLogger
from sqlalchemy import and_, create_engine
from sqlalchemy.orm import sessionmaker


RDLogger.DisableLog("rdApp.*")


def check_rgroup(data):
    ps = RGroupDecompositionParameters()
    ps.allowNonTerminalRGroups = True

    res = None
    mol = Chem.MolFromInchi(data["inchi"])

    if mol:
        matched, unmatched = rdRGD.RGroupDecompose(
            data["targets"],
            [mol],
            asSmiles=True,
            asRows=True,
            options=ps,
        )
        if len(unmatched) < 1:
            res = dict(cid=data[cid], inchi=data["inchi"])
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-chebi-csv", required=True, help="Rhea csv file")
    parser.add_argument(
        "--input-pubchem-db", required=True, help="Pubchem sql database"
    )
    parser.add_argument(
        "--parameter-chebi-int", required=True, type=int, help="ChEBI id int"
    )
    parser.add_argument("--output-search-txt", required=True, help="Output txt file")
    parser.add_argument("-t", "--input-thread-int", type=int, default=1, help="Threads")
    args = parser.parse_args()

    # Init
    threads = args.intput_thread_int

    engine = create_engine(f"sqlite:///{args.input_pubchem_db}")
    Session = sessionmaker(bind=engine)
    session = Session()

    # Get Mol
    df = pd.read_csv(args.input_chebi_csv)
    smiles = None
    for _, row in df.iterrows():
        if f"CHEBI:{args.parameter_chebi_int}" in row["chebi"]:
            smiles = row["smiles"]
            break
    assert (
        smiles is not None
    ), f"Smiles not retrieved for ChEBI: {args.parameter_chebi_int}"

    # Find R-group
    mol = Chem.MolFromSmiles(smiles)
    wt_query = Descriptors.ExactMolWt(mol)

    queries = (
        session.query(Compound.cid, Compound.inchi)
        .filter(
            and_(Compound.is_biochemical == 1, Compound.exact_mol_wt > exact_mol_wt)
        )
        .yield_per(int(1e7))
    )
    results = []
    batchs = set()
    targets = [mol]
    sub_start = time.time()
    for count, query in enumerate(queries):
        batchs.add(query)
        if len(batchs) > 1e6:
            with multiprocessing.Pool(processes=threads) as pool:
                async_results = [
                    pool.apply_async(
                        check_rgroup,
                        (dict(targets=targets, cid=query[0], inchi=query[1]),),
                    )
                    for query in batchs
                ]
                for async_result in async_results:
                    try:
                        results.append(async_result.get(timeout=10))
                    except multiprocessing.TimeoutError:
                        print("Timeout error")
            batchs = set()
    with multiprocessing.Pool(processes=threads) as pool:
        async_results = [
            pool.apply_async(
                check_rgroup, (dict(targets=targets, cid=query[0], inchi=query[1]),)
            )
            for query in batchs
        ]
        for async_result in async_results:
            try:
                results.append(async_result.get(timeout=10))
            except multiprocessing.TimeoutError:
                print("Timeout error")

    session.close()
    # Dereplicate

    datas = defaultdict(list)
    for result in results:
        mol = Chem.MolFromInchi(result["inchi"])
        if mol:
            smiles = Chem.MolToSmiles(mol)
            datas[smiles].append(result["cid"])

    # Write output
    with open(args.output_search_txt, "w") as fd:
        for smiles, cids in datas.items():
            cid = ",".join(cids)
            fd.write(f"{smiles}|{cid}\n")
