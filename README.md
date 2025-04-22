# ChEBI R-group

## Installation

```sh
conda env create --file recipes/worklow.yaml --name chebi-rgroup
pip install --no-deps -e .
```

## Build dataset

### Download PubChem
```sh
python -m chebirgroup.pubchem.download \
    --output-pubchem-dir <output dir> \
    --output-pubchem-db <output sql database>
```

### Download Rhea
```sh
python -m chebirgroup.rhea.download \
    --output-rhea-dir <output dir> \
    --parameter-release-int <int>
```

### R-group search
```sh
snakemake -n --config input_chebi_csv=$work/rgroup/rhea/rhea-chebi-smiles.csv input_pubchem_db=$work/rgroup/pubchem.db output_dir_str=$work/rgroup/chebi
```
