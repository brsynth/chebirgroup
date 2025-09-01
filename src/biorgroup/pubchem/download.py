import argparse
import glob
import gzip
import json
import os
import xml.etree.ElementTree as ET
from multiprocessing.pool import ThreadPool as Pool
from typing import Any, Dict, List, Optional, Tuple

from bs4 import BeautifulSoup
from biorgroup.pubchem.compound import Compound
from biorgroup.utils import cmd, atomic
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

RDLogger.DisableLog("rdApp.*")


class Prop(object):
    def __init__(
        self,
        label: Optional[str] = None,
        name: Optional[str] = None,
        *args,
        **kwargs,
    ) -> None:
        super(Prop, self).__init__(*args, **kwargs)
        self.label = label
        self.name = name

    def to_dict(self) -> Dict[str, Optional[str]]:
        return dict(
            label=self.label,
            name=self.name,
        )

    def __repr__(self) -> str:
        return str(self.to_dict())


class Pubchem(object):
    URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/"
    XML_NS = {
        "ncbi": "http://www.ncbi.nlm.nih.gov",
        "ncbi_complete": "{http://www.ncbi.nlm.nih.gov}",
    }

    def __init__(
        self,
        path: str,
        *args,
        **kwargs,
    ) -> None:
        super(Pubchem, self).__init__(*args, **kwargs)
        self.path = path

    def download(self) -> None:
        os.makedirs(self.path, exist_ok=True)
        # Grep root file
        index_path = os.path.join(self.path, "index.html")
        cmd.url_download(self.URL, index_path)
        with open(index_path) as fid:
            soup = BeautifulSoup(fid, features="html.parser")
        regular_files = []
        md5_files = []
        for item in soup.find_all("a"):
            value = item.getText()
            if value.startswith("Compound"):
                if value.endswith("md5"):
                    md5_files.append(value)
                else:
                    regular_files.append(value)
        # Download md5 files
        for filename in md5_files:
            if not os.path.isfile(os.path.join(self.path, filename)):
                cmd.url_download(self.URL + filename, os.path.join(self.path, filename))
        # Download regular file
        count = 0
        for filename in regular_files:
            md5_file = [x for x in md5_files if x.startswith(filename)][0]
            with open(os.path.join(self.path, md5_file)) as fid:
                data = fid.read()
                target_md5 = data.split()[0]
            current_md5 = ""
            if os.path.isfile(os.path.join(self.path, filename)):
                current_md5 = cmd.md5(path=os.path.join(self.path, filename))
            current_try = 0
            while current_md5 != target_md5:
                cmd.url_download(self.URL + filename, os.path.join(self.path, filename))
                current_md5 = cmd.md5(path=os.path.join(self.path, filename))
                current_try += 1
                if current_try > 3:
                    break
            count += 1
            if count > 1:
                break

    def create_csv(self, threads: int = 1) -> None:
        def _parse_xml(filename: str) -> None:
            fcsv = filename.replace(".xml.gz", ".csv.gz")
            if os.path.isfile(fcsv):
                return
            fid = gzip.open(filename, "r")
            compounds = []
            prop = Prop()
            data = {}
            for event, elem in ET.iterparse(fid, events=("start", "end")):
                tag = elem.tag.replace(self.XML_NS["ncbi_complete"], "")
                if event == "end" and tag == "PC-Compound" and len(data) > 0:
                    compounds.append(data)
                    data = {}
                    elem.clear()
                elif event == "start" and tag == "PC-Compound":
                    data = {}
                elif event == "end" and tag == "PC-CompoundType_id_cid":
                    data["cid"] = int(elem.text)
                elif event == "end" and tag == "PC-InfoData":
                    prop = Prop()
                elif event == "end" and tag == "PC-Urn_label":
                    prop.label = elem.text
                elif event == "end" and tag == "PC-InfoData_value_sval":
                    if prop.label == "InChI":
                        data["inchi"] = elem.text
            df = pd.DataFrame.from_dict(compounds)
            df.to_csv(fcsv, index=False)
            fid.close()

        pool = Pool(threads)
        pool.map(_parse_xml, glob.glob(os.path.join(self.path, "*.xml.gz")))

    def create_sql(self, path: str) -> None:
        engine = create_engine(f"sqlite:///{path}")
        Session = sessionmaker(bind=engine)
        session = Session()
        Compound.__table__.create(engine)

        for filename in glob.glob(os.path.join(self.path, "*.csv.gz")):
            df = pd.read_csv(filename)
            for ix, row in df.iterrows():
                mol = Chem.MolFromInchi(row["inchi"])
                if mol:
                    df.at[ix, "is_biochemical"] = atomic.is_biochemical(mol=mol)
                    df.at[ix, "exact_mol_wt"] = Descriptors.ExactMolWt(mol)
            df = df[~pd.isna(df["is_biochemical"])]
            session.add_all([Compound(**record) for record in df.to_dict("records")])
            session.commit()

        session.close()

    def to_dict(self) -> Dict[str, str]:
        return dict(path=self.path)

    def __repr__(self) -> str:
        return str(self.to_dict())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-pubchem-dir", required=True, help="Pubchem download directory"
    )
    parser.add_argument(
        "--output-pubchem-db", required=True, help="Pubchem sqlite database"
    )
    parser.add_argument(
        "-t", "--parameter-thread-int", type=int, default=1, help="Threads"
    )
    args = parser.parse_args()

    # Create database
    threads = args.parameter_thread_int
    output_pubchem_dir = args.output_pubchem_dir
    os.makedirs(output_pubchem_dir, exist_ok=True)

    pubchem = Pubchem(path=output_pubchem_dir)
    pubchem.download()
    pubchem.create_csv(threads=threads)
    pubchem.create_sql(path=args.output_pubchem_db)
