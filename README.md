# BioRGroup dataset

## Data location

[https://doi.org/10.57745/V3URYA](https://doi.org/10.57745/V3URYA)

## Installation

```sh
conda env create --file recipes/worklow.yaml --name biorgroup
pip install --no-deps -e .
```

## Build dataset

### 1 - Download PubChem
```sh
python -m biorgroup.pubchem.download \
    --output-pubchem-dir <output dir> \
    --output-pubchem-db <output sql database>
```

### 2 - Download Rhea
```sh
python -m biorgroup.rhea.download \
    --output-rhea-dir <output dir> \
    --parameter-release-int <int>
```

### 3 - R-group search
```sh
snakemake \
    -p \
    -j 48 \
    -c 48 \
    --workflow-profile template/biorgroup \
    -s ./src/biorgroup/rgroup/Snakefile \
    --use-conda \
    --latency-wait 5 \
    --rerun-incomplete \
    --config input_depot_str=./src/biorgroup/rgroup input_chebi_csv=rhea-chebi-smiles.csv input_pubchem_db=pubchem.db output_dir_str=chebi parameter_search_timeout_int=10
```

## Dataset overview
The Snakemake workflow produces a `csv.gz` file containing:  

| column name | type |
| --- | --- |
| smiles_rhea | `str`|
| chebi | `List[str]` |
| num_heavy_atoms | `int` |
| exact_mol_wt | `float` |
| core_superstructure_smiles | `List[str]` |
| core_superstructure_pubchem_cid | `List[List[str]]` |
| rgroup_extended_smiles | `List[str]` |
| rgroup_extended_pubchem_cid | `List[List[str]]` |
