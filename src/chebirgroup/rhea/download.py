import argparse
import os
import shutil
import tarfile
import tempfile
from typing import Any, Dict, List, Optional, Tuple

from chebirgroup.utils import cmd
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm


class Rhea(object):
    URL = "https://ftp.expasy.org/databases/rhea/old_releases"

    def __init__(
        self,
        path: str,
        release: int,
        *args,
        **kwargs,
    ) -> None:
        super(Rhea, self).__init__(*args, **kwargs)
        self.path = path
        self.release = release

    def download(self) -> None:
        url = self.URL + f"/{self.release}.tar.bz2"
        tmp = tempfile.NamedTemporaryFile(suffix=".tar.bz2")
        cmd.url_download(url=url, path=tmp.name)
        target_file = f"{self.release}/tsv/rhea-chebi-smiles.tsv"
        with tarfile.open(tmp.name, "r:bz2") as tar:
            tar.extract(target_file, path=self.path)
        shutil.move(
            os.path.join(self.path, target_file),
            os.path.join(self.path, "rhea-chebi-smiles.tsv"),
        )
        shutil.rmtree(os.path.join(self.path, f"{self.release}"))

    def format(self) -> None:
        df = pd.read_csv(
            os.path.join(self.path, "rhea-chebi-smiles.tsv"),
            sep="\t",
            names=["chebi", "smiles"],
        )
        df = df.groupby("smiles")["chebi"].apply(list).reset_index()

        # Check Rhea
        df["num_heavy_atoms"] = 0
        for ix, row in tqdm(df.iterrows(), total=df.shape[0]):
            rhea_smiles = row["smiles"]
            mol = Chem.MolFromSmiles(rhea_smiles)
            rhea_smiles = Chem.MolToSmiles(mol)

            df.at[ix, "smiles"] = rhea_smiles
            df.at[ix, "num_heavy_atoms"] = mol.GetNumHeavyAtoms()
            df.at[ix, "exact_mol_wt"] = round(Descriptors.ExactMolWt(mol), 2)
            df.at[ix, "chebi"] = sorted(list(set(row["chebi"])))
        df.to_csv(os.path.join(self.path, "rhea-chebi-smiles.csv"), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir-str", required=True, help="Output for the Rhea directory"
    )
    parser.add_argument(
        "--parameter-release-int", type=int, default=134, help="Release version"
    )
    args = parser.parse_args()

    # Create database
    output_dir_str = args.output_dir_str
    os.makedirs(output_dir_str, exist_ok=True)
    rhea = Rhea(path=output_dir_str, release=args.parameter_release_int)
    rhea.download()
    rhea.format()
