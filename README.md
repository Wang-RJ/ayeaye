# Unprecedented Female Mutation Bias in the Aye-Aye: Analysis Repository

This repository contains the code and resources used for the analysis described in the manuscript **"Unprecedented female mutation bias in the aye-aye."** The analysis focuses on downstream processing of sequencing data (VCF and BAM files) for two species: the aye-aye and the baboon.

## Repository Structure

The repository is organized as follows:

- **Top-level files**:
  - Three `.txt` files containing identified de novo mutations from:
    - Aye-aye analyses
    - Baboon analyses
    - Reanalyses of baboon samples published in Wu et al.

- **Subfolders**:
  - `ayeaye`: Contains data and scripts for the aye-aye analyses, including sex bias studies.
  - `baboon`: Contains data and scripts for the baboon analyses.
  - `wu_baboon`: Contains data and scripts for reanalysis of the Wu et al. data.

### Subfolder Contents

Each subfolder includes:

- **Shell Scripts**:
  - Located in the `shell` subfolder, these scripts perform callability and depth analyses.
  - Scripts in the aye-aye folder also include analyses related to:
    - Phasing with **POOHA** and **Unfazed**.
    - Sex bias, as detailed in the manuscript.

- **R Scripts**:
  - Located in the main folder of each analysis:
    - `callability.R`: Calculates callability from tables created in shell scripts.
    - `denovo_<samples>.R`: Main script for identifying de novo mutations.
    - `denovo_utilities.R`, `processBAM.R`, `processRelated.R`: Scripts with helper functions.

## Dependencies

### Shell Code

The shell scripts depend on the following tools and versions:

- `gcc/9.3.0`
- `python/3.9.8`
- `bcftools/1.17`

Additionally, the shell scripts call two bioinformatics phasing programs:

- **POOHA**
- **Unfazed**

### R Code

The R scripts require **R/4.0.2** and the following R packages:

- `dplyr`
- `ggplot2`
- `memoise`
