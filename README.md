# Path-Explorer

Path-Explorer is a new suite of user-friendly computational tools that enable reconstruction of co-occurring strains, as well as provide other fundamental analytical capabilities that will enhance metagenomic-based WBS of various pathogen and antimicrobial resistance targets.

<div align="center">
	<img width="558" alt="new" src="https://github.com/haoqiusong/Path-Explorer/assets/106828678/c8ff02ff-e917-4a92-896d-386fe066bd78">
</div>

# Installation & Usage

## Requirements

1. Linux Operating System
2. Conda

## Installation

```
git clone https://github.com/haoqiusong/Path-Explorer.git
conda env create -f environment.yml
conda activate path-explorer
```

## (A). Read mapping to quantify pathogen and ARG markers of interest

```
scp install.sh read_mapping_and_quantification
cd read_mapping_and_quantification
bash install.sh
```

### Download the compressed standard 16 GB kraken2 DB and uncompress it

```
mkdir k2_DB
cd k2_DB
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240112.tar.gz
tar -zxvf k2_standard_16gb_20240112.tar.gz
rm k2_standard_16gb_20240112.tar.gz
```

## (B). De novo assembly to assess potential new pathogen variants and to determine “resistome risk” (i.e., degree to which ARGs have propensity to spread to pathogens)

## (C). Strain-resolved pathogen detection
