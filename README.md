# AnnotateSnakeMake

[![DOI](https://zenodo.org/badge/542176125.svg)](https://zenodo.org/doi/10.5281/zenodo.11084022)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.12.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


**AnnotateSnakeMake** is a Snakemake workflow routinely used in the [Marlétaz lab](https://fmarletaz.github.io/) to **annotate genomes** of echinoderm species. It takes as input a **genome sequence** (*fasta*) and **RNA-seq datasets** (*fastq.gz*, paired-end) and outputs the set of annotated **protein-coding genes** (*bed*, *gtf* and *fasta* files).


## Table of content
  - [Installation](#installation)
    - [Conda and Snakemake](#conda-and-snakemake)
    - [AnnotateSnakeMake](#annotatesnakemake)
  - [Usage](#usage)
  	- [Workflow breakdown](#workflow-breakdown)
    - [Test set](#test-set)
    - [User-provided data](#user-provided-data)
  - [Contact](#contact)
  - [Reference](#reference)

## Installation

### Conda and Snakemake

The pipeline uses Conda to deploy all of its dependencies automatically. The recommended way to install Snakemake is thus to follow the Conda/Mamba installation guidelines detailed in the [Snakemake documentaion](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

We recommend creating a conda environment with snakemake 7.12 and tabulate 0.8:

```
mamba create -n snakemake -c bioconda -c conda-forge snakemake=7.12 tabulate=0.8
```

### AnnotateSnakeMake

Clone the repository:

  ```
  git clone https://github.com/eparey/AnnotateSnakeMake.git
  ```

## Usage

### Worflow breakdown

The following steps are ~ sequentially ran (see workflow image in `images/workflow.pdf`):

- repeat annotation and masking (RepeatModeler, RepeatMasker): `workflow/rules/repeats.smk`
- transcriptome assembly, de novo (Trinity, gmap): `workflow/rules/transcriptome_de_novo.smk`
- transcriptome assembly, genome-guided (STAR, STRINGTIE, TACO): `workflow/rules/transcriptome_using_genome.smk`
- consensus transcriptome (Mikado): `workflow/rules/consensus_transcriptome.smk`
- model training and gene prediction with AUGUSTUS (hints from Metaeuk, Mikado and Portcullis): `workflow/rules/gene_pred.smk`
- repeat filtering from gene models: `workflow/rules/filter_repeats.smk`
- gene models refinement with PASA: `workflow/rules/gene_annot.smk`
- evaluation of the annotation (BUSCO, PFAM), gene naming (DIAMOND, SWISSPROT) and file formatting: `workflow/rules/final_models.smk`

### Test set

The test sets serves as an example for input data specification and formatting (see `config/config.yaml` and `resources/`).

Before running, users need to download and format the pfam database, and specify its location in `config/config.yaml` (please see download instructions at https://github.com/aziele/pfam_scan, 'Prepare library of Pfam HMMs').

- The pipeline has to be run from the `AnnotateSnakeMake` folder with the snakemake environment activated (see above to create the 'snakemake' conda env):

	```
	cd AnnotateSnakeMake
	conda activate snakemake
	```

- To print a dry-run:

	```
	snakemake --configfile=config/config.yaml --use-conda --cores 48 -n
	````

- To run on the test set:

	```
	snakemake --configfile=config/config.yaml --use-conda --cores 48
	````


### User-provided data

- Create a specific directory for the run:

	```
	mkdir ../Afil_annotation
	```

- Create a configuration file. We recommend copying and editing the example config (`cp config/config.yaml ../Afil_annotation/config_afil.yaml` ). Note that relative paths should be replaced by absolute paths and paths to the user-provided data.

- Run the pipeline, still from within the `AnnotateSnakeMake/` folder but providing the run directory with `--directory`:

	```
	conda activate snakemake
	snakemake --configfile ../Afil_annotation/config_afil.yaml --directory ../Afil_annotation/ --use-conda --cores 48
	```


## Contributors

The workflow was primarily deposited for reproducibility purposes, but if you are interested in running it but are encountering issues, please do not hesitate to contact us, we are happy to help:

- [Elise Parey](e.parey@ucl.ac.uk)
- [Ferdinand Marlétaz](f.marletaz@ucl.ac.uk)

We also thank Chema Martin-Duran for contributing to the development of the pipeline.

## References

The annotation workflow is described in the methods section of the brittle star *Amphiura filiformis* genome paper:

Parey, E. et al. The brittle star genome illuminates the genetic basis of animal appendage regeneration. Nat Ecol Evol (2024). https://doi.org/10.1038/s41559-024-02456-y
