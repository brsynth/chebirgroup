# ChEBI R-group

## Installation

```sh
conda env create --file recipes/worklow.yaml --name chebirgroup
pip install --no-deps -e .
```

## Build dataset

### 1 - Download PubChem
```sh
python -m chebirgroup.pubchem.download \
    --output-pubchem-dir <output dir> \
    --output-pubchem-db <output sql database>
```

### 2 - Download Rhea
```sh
python -m chebirgroup.rhea.download \
    --output-rhea-dir <output dir> \
    --parameter-release-int <int>
```

### 3 - R-group search
```sh
snakemake \
    -p \
    -j 48 \
    -c 48 \
    --workflow-profile template/chebirgroup \
    -s ./src/chebirgroup/rgroup/Snakefile \
    --use-conda \
    --latency-wait 5 \
    --rerun-incomplete \
    --config depot_dir=./src/chebirgroup/rgroup input_chebi_csv=rhea-chebi-smiles.csv input_pubchem_db=pubchem.db output_dir_str=chebi
```

## Dataset overview
The Snakemake workflow produces a `csv.gz` file containing:  
| smiles_rhea | chebi | num_heavy_atoms | exact_mol_wt | additional_substituents_smiles | additional_substituents_pubchem_cid | only_match_rgroup_smiles | only_match_rgroup_pubchem_cid |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `str` | `List[str]` | `int` | `int` | `List[str]` | List[List[str]] | List[str] | List[List[str]] |
