# README
<!-- toc -->

- [P.knowlesi de novo genome assembly](#Pknowlesi-de-novo-genome-assembly)
  * [Abstract](#Abstract)
  * [How to Cite](#How-to-Cite)
  * [Table of Contents](#Table-of-Contents)
  * [Data Access](#Data-Access)
  * [Installation](#Installation)
    + [Processing](#Processing)
    + [Conda](#Conda)
    + [Other tools not on Conda](#Other-tools-not-on-Conda)
    + [Pipeline](#Pipeline)
      - [Preparation](#Preparation)

<!-- tocstop -->
# P.knowlesi de novo genome assembly
## Abstract
A pipeline to take input data from raw nanopore reads to being ready for manual or automated annotation. 
It is important to note that while this pipeline automates alot of commands, it cannot automate everything. 
User input is needed at various points to confirm which branch of the pipeline to follow as well as when to continue at certain 'break-points'
I'm hoping that this has been fairly prepped in a fairly straight-forward manner that is somewhat approachable
However, if you want to tweak the script, the vast majority is written in Bash with some secondary scripts in python.
Please be aware that the first portion of this script is tuned for barcoded ONT sequencing reads however, it is possible to use non-barcoded ONT sequenced read. 

## How to Cite

## Table of Contents
- Abstract
- How to Cite
- Table of Contents
- Data Access
- Installation
  - Processing
	- Conda YAML files
	- Other Tools not on Conda
- Pipeline
	- Preparation
	- Basecalling
	- Demultiplexing
	- Adapter Removal
	- Checks
	- Alignment
	- De novo assembly
	- Genome size estimation
	- Decontamination
	- Racon Polishing
	- Medaka Correction
	- Pilon Correction
	- Remove API and MIT
	- RepeatMasking
	- De-chrimerisation
	- Quality Assessment
		- QUAST
		- BUSCO
		- Pomoxis
		- Assembly-stats
	- Ready for annotation
	- FAQ
  	- I haven't used barcodes. What do I do?
	- Acknowledgements
	- References
## Data Access
The completed genomes generated from this pipeline was uploaded to ... with accession code...
Requirements:
## Installation
### Processing
- You will need a GPU in order to carry out basecalling! CPU only basecalling while available for Guppy is simply too slow for this pipeline. This pipeline was developed using an NVIDIA GTX1060 with 6GB of VRAM. 
- Other NVIDIA GPUs are likely to work as long as 

### Conda
- First and foremost, you will need conda in order to install the vast majority of tools used in this pipeline
- The yaml file "pipeline.yml" contains a list of tools and the versions which have been tested with this build of the pipeline.
- Install with `conda env update -n test_env --file pipeline.yml`
- In some cases, the conda yaml installation may not work due to unsolvable conflicts. If this occurs, here are the tools and versions you can download from cond:a
  - assembly-stats v1.0.1
  - bedtools v2.30
  - bioawk v1.0
  - blobtools v1.1.1
  - cd-hit v4.8.1
  - flye v2.8.3
  - genometools v1.6.1
  - jellyfish v2.2.10
  - minimap2 v2.17
  - pigz v2.5
  - pilon v1.23
  - pomoxis v0.2.5
  - porechop v0.2.4
  - pv v1.6.6
  - qcat v1.1.0
  - quast v5.0.2
  - racon v1.4.20
  - ragtag v1.1.1
  - repeatmasker v4.1.1
  - repeatmodeler v2.0.1
  - samtools v1.11
  - seqkit v0.15.0
### Other tools not on Conda
- **Guppy**: Down and install the version which suits your preference from the Nanopore Community site. This pipeline has been tested with **version 4.0.15**
- **Busco** docker image:
``` bash
docker pull ezlabgva/busco:v5.0.0_cv1
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.0.0_cv1
```
- Download and install BLAST+ package: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- Ensure that the BLAST bin and BLAST databases (particularly nt) are in your $PATH

- These packages have been tested to work together although more recent version may result in breakage
- There are also tools which cannot be found on Conda
  - CENSOR: https://www.girinst.org/censor/
  - OneCodeToFindThemAll
  - TransposonPSI
  - RepBase: https://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz
    - You will need an academic license for this
  - Companion
### Pipeline
#### Preparation
- Open the `Whole_pipeline.sh` and enter the paths to the different files, folders indicated:
  - **SCPTS**: Enter the full path to where the package was downloaded and unzipped so that secondary scripts can be called e.g. `/path/to/PknowlesiDeNovo/package/Scripts`
  - **RPATH**: If you are starting from raw sequencing reads, enter the full path to the raw reads e.g `/path/to/ONT/raw/reads/pass/`. Leave blank otherwise
  - **BSECL**: If your data has already been basecalled then enter **yes**. Otherwise **no**
  - **DPATH**: **If starting from raw reads, leave blank**. If starting from basecalled reads, enter full path here e.g. `/path/to/basecalled/data/folder/experimentFolder`
  - **ILLUMPATH**: If you have access to Illumina short reads for the ONT sequenced isolates, enter the full path to the folder containing the R1 and R2 fasta/q here e.g `/path/to/the/Illumina/read/folder`
  - **SPATH**: Enter the full path to the save folder that will hold all the output folders and files. Note that, if this folder was not previously created, it will be created here also. `/full/path/to/save/folder`
  - **REFPATH**: Enter the full path to the reference genome you are using to remove contamination. This pipeline is tuned to align your raw reads against the human reference to remove human contaminant sequences. Path would be `/path/to/Human/reference/file.fasta (or fna, fa etc)`
  - **PKREF**: Enter the full path to the desired Plasmodium reference fasta file. It will only be used for comparison **after** *de novo* assembly `/path/to/reference.fasta`
  - **REFGFF**: Enter the full path to the human reference gff `/path/to/human/reference.gff3`
  - **PKGFF**: Enter the full path to the stated Plasmodium reference GFF file `/path/to/reference.gff3 (or gff)`
  - **GENOME**: Enter the estimated genome size for the genome of interest
  - **EXPMT**: Give the name of the experiment or give this analysis run a name e.g. `Pk_April2021_Clinical_Isolate`
  - **BCODES**: Give the barcodes that were used to label the isolates in quotations, separated by space e.g `("barcode05" "barcode06" "barcode07")`
  - **ISOLATES**: What are the isolate names/IDs for the isolates used. Please enter the isolates in order of the barcodes used in quotations, separated by space. i.e. barcode05 is sks047. e.g. `("sks047" "sks047" "sks048" "sks048" "sks058")`
  - **DUPL**: Are there duplicated isolates? i.e. in this run, has a single isolate being sequenced with multiple barcodes. **Yes or No**
  - **DUPPY**: Enter the duplicate isolates in quotations, separated by space. `("sks047" "sks048)"`
  - **ILLUMISO**: If Illumina short reads are available, enter the isolates which have Illumina reads in quotation, separated by space `("sks047" "sks048)"`
  - **KIT**: The sequencing kit ID used e.g. `"SQK-RBK004"`. Check for other kits on the terminal with:  `/path/to/ont-guppy/bin/guppy_basecaller --print_workflows`
  - **KIT2**: The sequencing kit used however using the qcat nomenclature. You can find this using: `qcat --list-kits` and choosing the correct kit which you have used. 
  - **Note: Both KIT and KIT2 are compulsory**
  - **FLOWCELL**: Enter the flowcell version you have used e.g. `"FLO-MIN106"`. This is needed for basecalling
  - **THREADS**: Number of threads to use e.g `"32"`
  - **KMER**: Kmer length to use to process files. However, this is only necessary if you are carrying out genomescope/jellyfish of your Illumina reads. *Note: kmer length should only be odd numbers. This pipeline was tuned with a kmer length of 27* e.g. `"27"`
  - **SPEC**: State the species (or closest species) that has been trained on Augustus e.g. `"pfalciparum"`. List can be found here: https://funannotate.readthedocs.io/en/latest/commands.html#funannotate-species` 
  - **ONT**: Enter the full path the the 'bin' of the guppy package you have downloaded e.g. `/path/to/ont-guppy/bin`
  - **blastDATA**: 
  - **busDock**: Enter full path to the docker image for the busco tool e.g. `"/path/to/busco_v5.0.0_cv1.sif"`
  - **LINEAGE**: Enter the lineage library you wish to use for busco e.g. `"plasmodium_odb10"` To see the full list available, go to: https://busco-data.ezlab.org/v5/data/lineages/  
