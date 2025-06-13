# ChEBI R-group

## Data location

[https://doi.org/10.57745/V3URYA](https://doi.org/10.57745/V3URYA)

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

| column name | type |
| --- | --- |
| smiles_rhea | `str`|
| chebi | `List[str]` |
| num_heavy_atoms | `int` |
| exact_mol_wt | `int` |
| additional_substituents_smiles | `List[str]` |
| additional_substituents_pubchem_cid | `List[List[str]]` |
| only_match_rgroup_smiles | `List[str]` |
| only_match_rgroup_pubchem_cid | `List[List[str]]` |
