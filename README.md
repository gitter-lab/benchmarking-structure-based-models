# Exploring zero-shot structure-based protein fitness prediction
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13821572.svg)](https://doi.org/10.5281/zenodo.13821572)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13819823.svg)](https://doi.org/10.5281/zenodo.13819823)

This is a repository of code needed to replicate benchmarks of structure-based protein fitness prediction models in ProteinGym.

[Exploring zero-shot structure-based protein fitness prediction](https://openreview.net/forum?id=DtpbhlaNdb).  
Arnav Sharma, Anthony Gitter.  
*Generative and Experimental Perspectives for Biomolecular Design Workshop at the 13th International Conference on Learning Representations*. 2025.

## Table of Contents
  * [Setup](#setup)
    * [Conda environment](#conda)
  * [Accessing Data](#access_data)

## Setup
### Conda environment setup
The conda environment can be setup using 
```
  conda env create -f environment.yml
```
### SSEmb setup
* Download the `test.tar.gz` file from the [Zenodo dataset repository](https://doi.org/10.5281/zenodo.13819823) and place it in the assets directory of this repo.
* Use `scripts/setup_ssemb.sh` to setup the data and code.
* Follow instructions on the [SSEmb Github](https://github.com/KULL-Centre/_2023_Blaabjerg_SSEmb) to run ProteinGym assays.

### ProteinGym setup for running ESM-IF1 with experimental structures

* Download the `experimental_struct_artifacts.tar.gz` file from the [Zenodo dataset repository](https://doi.org/10.5281/zenodo.13819823) and place it in the assets directory of this repo.
* Use `scripts/setup_proteinGym.sh` to setup the data and code.
* Follow instructions on the [ProteinGym Github](https://github.com/OATML-Markslab/ProteinGym) to setup `scipts/zero_shot_config.sh` within ProteinGym.
* cd into the cloned `ProteinGym/scripts/scoring_DMS_zero_shot directory`. Run `sh score_multichain_structs_esm.sh` to generate ESM-IF1 scores for experimental
structures.

## Accessing Data:

### ProteinGym Data:
Get zero-shot scores:
```
wget https://marks.hms.harvard.edu/proteingym/zero_shot_substitutions_scores.zip
```
These scores can be extracted to the location of your choosing and then each script can be pointed to this location.

### Generated Data:
Data to run experiments can be found in our [Zenodo dataset repository](https://doi.org/10.5281/zenodo.13819823).

  * `experimental_struct_artifacts.tar.gz` contains the experimentally determined structures
for ProteinGym assays used in our analysis along with the reference file needed to
generate ESM-IF1 predictions for these structures in ProteinGym. This compressed file contains:
    * Experimental structures for ProteinGym assays. These structures were filtered using logic similar to that used in
    the ["Structure-informed protein engineering with equivariant graph neural networks" repository](https://github.com/semiluna/partIII-amino-acid-prediction).
    * Modified files for running ESM-IF1 on multichain complexes. The function is similar to that in the [ESM repository](https://github.com/facebookresearch/esm/blob/main/examples/inverse_folding/score_log_likelihoods.py)
    * Reference file to map all DMS ids to their pdb files
    * Script to run inference for each DMS id in the reference file.

  * `results.tar.gz` contains the prediction results obtained by running SSEmb on the 216
ProteinGym assays being considered in this study. Files compressed are:
    * `df_total_proteingym_1.csv` which has all the SSEmb predictions for the 216 ProteinGym assays.
    * `experimental_struct_scores.csv` which has all the ESM inverse folding predictions for the 61 assays with experimental
      structures.
  * `test.tar.gz` contains all the structures from ProteinGym as well as MSAs generated
using mmseqs2. To use this directory:
    * Setup SSEmb as directed in its [repository](https://github.com/KULL-Centre/_2023_Blaabjerg_SSEmb)
    * Download this file and extract it in the data folder.
