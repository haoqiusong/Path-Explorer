# Path-Explorer

Path-Explorer is a new suite of user-friendly computational tools that enable reconstruction of co-occurring strains, as well as provide other fundamental analytical capabilities that will enhance metagenomic-based WBS of various pathogen and antimicrobial resistance targets.

<div align="center">
	<img width="1006" alt="Screenshot 2025-02-14 at 2 10 03 PM" src="https://github.com/user-attachments/assets/b26ebd44-3433-426b-83dd-7dc7d60a6e47" />

</div>

# Installation

## Environment installation

```
git clone https://github.com/haoqiusong/Path-Explorer.git
conda env create -f environment.yml
conda activate path-explorer
```

Please make sure you have a Conda installed on your machine.

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

```
scp install.sh de_novo_assembly_and_resistome_risk
cd de_novo_assembly_and_resistome_risk
bash install.sh
```

### Download the compressed Blast Database file from Zenodo (25 GB) to run MetaCompare and uncompress it

```
wget https://zenodo.org/records/10471551/files/BlastDB.tar.gz
tar -zxvf BlastDB.tar.gz
```

### Download the compressed DeepARG-DB and mobileOG database (DB.tar.gz), put it inside "de_novo_assembly_and_resistome_risk" directory and uncompress it

Use this [link](https://drive.google.com/file/d/10YuSxmre1bIg8V6-oAHjBkMcUWGJUmE4/view?usp=drive_link) to download **DB.tar.gz**. Put it inside the **de_novo_assembly_and_resistome_risk** directory and uncompress it using the following command:

```
tar -zxvf DB.tar.gz
```

## (C). Strain-resolved pathogen detection

```
scp install.sh strain_resolved_detection_and_recovery
cd strain_resolved_detection_and_recovery
bash install.sh
```

### Users w/ their own reference genome list

Make sure you have a sequence data ( * .fasta/ * .fa/ * .fasta.gz/ * .fa.gz) containing all your reference genome sequences.

### Users w/o their own reference genome list

A default list of reference genomes will be used, which contains multiple pathogens at different species and strain level. In this list, three genomes of each ESKAPEE (Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, Acinetobacgter baumannii, Pseudomonas aeruginosa, Enterobacter spp., and Escherichia coli) pathogens are included. Besides, there are also several viruses including influenza virus (H1N1, H3N2, and H5N1), Norovirus (GI, GII, and GV), Coronavirus (OC43, NL63, 229E, and HKU1), and SARS-CoV-2 strains.

# Usage

## (A). Read mapping to quantify pathogen and ARG markers of interest

To run the assembly pipeline on metagenomic paired-end short read data ( * .fastq/ * .fq/ * .fastq.gz/ * .fq.gz), use the following command:
```
nextflow run short-read-pipeline.nf --R1 <absolute/path/to/forward/read/file> --R2 <absolute/path/to/reverse/read/file> --out_fname <prefix of output file name>
rm -r work
```

The command line options for this script (**short-read-pipeline.nf**) are:

**--R1**: The absolute path of the fastq file containing forward read sequences

**--R2**: The absolute path of the fastq file containing reverse read sequences

**--out_fname**: The prefix of the output file name

### Output

With **--out_fname S1**, output files include: **S1_rpoB_ARG_norm.tsv** and **S1_drug_wise_rpoB_norm.tsv**.

## (B). De novo assembly to assess potential new pathogen variants and to determine “resistome risk” (i.e., degree to which ARGs have propensity to spread to pathogens)

To run the assembly pipeline on metagenomic paired-end short read data ( * .fastq/ * .fq/ * .fastq.gz/ * .fq.gz), use the following command:

```
nextflow run assembly_pipeline.nf --R1 <absolute/path/to/forward/read/file> --R2 <absolute/path/to/reverse/read/file> --out_fname <prefix of output file name>
rm -r work
```

The command line options for this script (**assembly_pipeline.nf**) are:

**--R1**: The absolute path of the fastq file containing forward read sequences

**--R2**: The absolute path of the fastq file containing reverse read sequences

**--out_fname**: The prefix of the output file name

### Output

With **--out_fname S1**, output files include: **S1_resistome_risk.txt**, **S1_ARGs.faa**, and **S1_ARGs_and_mobility.tsv**.

## (C). Strain-resolved pathogen detection

To run the strain-resolved detection and recovery pipeline on metagenomic paired-end short read data ( * .fastq/ * .fq/ * .fastq.gz/ * .fq.gz), use the following command:

```
nextflow run strain_detection.nf --R1 <absolute/path/to/forward/read/file> --R2 <absolute/path/to/reverse/read/file> (--Ref <absolute/path/to/reference/genome/file>) --out_fname <prefix of output file name>
rm -r work
```

The command line options for this script (**strain_detection.nf**) are:

**--R1**: The absolute path of the fastq file containing forward read sequences

**--R2**: The absolute path of the fastq file containing reverse read sequences

**--Ref**: The absolute path of the fasta file containing all reference genome sequences (If not provided, a default reference genome list will be used.)

**--out_fname**: The prefix of the output file name

### Output

With **--out_fname S1**, output files include: **S1_recovered_genomes.fasta**, **S1_relative_abundance.png**, and **S1_genome_recovery_rate.tsv**.

# Citation

If you use Path-Explorer, please cite our work:
TBD...
