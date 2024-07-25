# reAlldock
reAlldock downloads all AlphaFold structures from the proteome of any species in the AlphaFold Database and docks ligands against all of them with Vina-GPU. After that, performs Enrichment analysis using String API for the top N docking score proteins.

## Requirement
- Currently, it supports only Linux.
- GPU is required and make sure the version of GPU driver is up to date.
- Sufficient free storage space over 1TB is recommended.

## Installation
1. Install [Vina-GPU](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1).

2. Clone current repo.
```
git clone https://github.com/toxtoxcat/reAlldock
cd reAlldock
```

3. Create Conda environment.
```
conda create -n realldock
conda activate realldock
conda install -c conda-forge biopython openbabel fpocket pandas matplotlib seaborn rdkit
conda install -c anaconda requests
```

## Usage
1. Download AlphaFold structures and prepare docking.
Running the code will display a list of species available for download from [AlphaFold Database](https://alphafold.ebi.ac.uk/). Enter the number of the species in the list, then download will begin. After download, Creating pdbqt files and Running fpocket will begin.
```
python main_prepare.py \
        --output_path {Absolute path of the output directory}
```

2. Docking with Vina-GPU.

```
python main_dock.py \
        --output_path {Absolute path of the output directory}
        --ligand_path {path of the ligand for docking to the AlphaFold structures}
        --vina_executable_path {Absolute path of the Vina-GPU executable file}
```

3. Perform Enrichment analysis using String API for the top N docking score proteins.

```
python main_string.py \
        --output_path {Absolute path of the output directory}
        --vina_csv_path {path of the Vina-GPU result csv}
```
