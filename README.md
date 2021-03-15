# README
<!-- toc -->

- [P.knowlesi de novo genome assembly](#Pknowlesi-de-novo-genome-assembly)
  * [Abstract](#Abstract)
  * [How to Cite](#How-to-Cite)
  * [Data Access](#Data-Access)
  * [Installation](#Installation)
    + [Processing](#Processing)
    + [Conda](#Conda)
    + [Other tools not on Conda](#Other-tools-not-on-Conda)
    + [Script and Pipeline Preparation](#Script-and-Pipeline-Preparation)
  * [Pipeline](#Pipeline)
    + [Data Preparation](#Data-Preparation)
      - [Basecalling](#Basecalling)
      - [Demultiplexing](#Demultiplexing)
      - [Adapter Removal](#Adapter-Removal)
      - [Checks](#Checks)
    + [De Novo Assembly](#De-Novo-Assembly)
      - [Alignment](#Alignment)
      - [De novo assembly](#De-novo-assembly)
      - [Genome size estimation](#Genome-size-estimation)
      - [Decontamination](#Decontamination)
    + [Correction and Polishing](#Correction-and-Polishing)
      - [Racon Polishing](#Racon-Polishing)
      - [Medaka Correction](#Medaka-Correction)
      - [Pilon Correction](#Pilon-Correction)
    + [[Optional *but recommended*]: Apicoplast and Mitochondrial Separation and Circlarisation](#Optional-but-recommended-Apicoplast-and-Mitochondrial-Separation-and-Circlarisation)
      - [Remove API and MIT](#Remove-API-and-MIT)
    + [RepeatMasking](#RepeatMasking)
      - [RepeatMasking](#RepeatMasking-1)
      - [De-chrimerisation](#De-chrimerisation)
    + [Quality Assessment](#Quality-Assessment)
      - [QUAST](#QUAST)
      - [BUSCO](#BUSCO)
      - [Pomoxis](#Pomoxis)
      - [Assembly-stats](#Assembly-stats)
    + [Ready for annotation](#Ready-for-annotation)
  * [FAQ](#FAQ)
    + [I haven’t used barcodes. What do I do?](#I-haven’t-used-barcodes-What-do-I-do)
  * [Acknowledgements](#Acknowledgements)
  * [References](#References)

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

## Data Access
The completed genomes generated from this pipeline was uploaded to ... with accession code...
## Installation
### Processing
- You will need a GPU in order to carry out basecalling! CPU only basecalling while available for Guppy is simply too slow for this pipeline. This pipeline was developed using an NVIDIA GTX1060 with 6GB of VRAM. 
- Other NVIDIA GPUs are likely to work as long as they have at least **NVIDIA driver version 455 and NVIDIA compute version 6.1**. It is advised that you follow the guppy installation guide on the Oxford Nanopore Community website

### Conda
- First and foremost, you will need conda in order to install the vast majority of tools used in this pipeline
- The yaml files "Main_pipeline.yml", "Medaka_env.yml", "Blobtools_env.yml", "Quast_env.yml", "Pomoxis_env.yml" contain lists of tools and the versions which have been tested with this build of the pipeline.
- Install with `conda env update -n test_env --file pipeline.yml`
- In some cases, the conda yaml installation may not work due to unsolvable conflicts. If this occurs, here are the tools and versions you can download from conda.
  - Note: I would suggest installing the tools in different environments as stated below in the order stated
  - Main pipeline environment
      | **Tools**      | **Version** |
      | -------------- | ----------- |
      | qcat           | 1.1.0       |
      | samtools       | 1.11        |
      | minimap2       | 2.17        |
      | pigz           | 2.5         |
      | assembly-stats | 1.0.1       |
      | porechop       | 0.2.4       |
      | pilon          | 1.23        |
      | racon          | 1.4.20      |
      | flye           | 2.8.3       |
      | jellyfish      | 2.2.10      |
      | bedtools       | 2.30        |
      | ragtag         | 1.1.1       |
      | bioawk         | 1.0         |
      | seqkit         | 0.15.0      |
  - Blobtools environment
      | **Tools** | **Version** |
      | --------- | ----------- |
      | blobtools | 1.1.1       |
      | minimap2  | 2.17        |
  - Medaka environment
      | **Tools** | **Version** |
      | --------- | ----------- |
      | medaka    | 1.1.1       |
  - Pomoxis environment
      | **Tools** | **Version** |
      | --------- | ----------- |
      | pomoxis   | 0.2.5       |
  - Quast environment
      | **Tools** | **Version** |
      | --------- | ----------- |
      | quast     | 5.0.2       | 
  - RepeatMasking environment
      | **Tools**      | **Version** |
      | -------------- | ----------- |
      | repeatmasker   | 4.1.1       |
      | repeatmodeler  | 2.0.1       |
      | cd-hit         | 4.8.1       |
      | genometools    | 1.6.1       |   

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
### Script and Pipeline Preparation
- Open the `Whole_pipeline.sh` and enter the paths to the different files, folders indicated:
  - **SCPTS**: Enter the full path to where the package was downloaded and unzipped so that secondary scripts can be called e.g. `/path/to/PknowlesiDeNovo/package/Scripts`
  - **tConEnv**: Enter the name of main pipeline e.g. `"Main_pipeline_env"`
  - **bConEnv**: Enter the name of the blobtools conda environment e.g. `"blobtools_env"`
  - **qConEnv**: Enter the name of the quast conda environment
  - **pConEnv**: Enter the name of the pomoxis conda environment
  - **mConEnv**: Enter the name of the medaka conda environment
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
  - **busDock**: Enter full path to the docker image for the busco tool e.g. `"/path/to/busco_v5.0.0_cv1.sif"`
  - **LINEAGE**: Enter the lineage library you wish to use for busco e.g. `"plasmodium_odb10"` To see the full list available, go to: https://busco-data.ezlab.org/v5/data/lineages/ 
## Pipeline
### Data Preparation
This section mainly covers initial parsing, basecalling, demultiplexing and adapter removal from the raw Nanopore reads. 
- First confirm if the data you wish to process has been previously basecalled.
  - This allows you to re-try the entire analysis without undertaking the basecalling step
#### Basecalling
- This pipeline is tuned for Guppy GPU basecalling. While it is possible to carry out Guppy basecalling on a CPU, I would advise not to attempt this with a Plasmodium whole genome as it is quite simply too slow.
- If you update your Guppy version between analyses, I would suggest re-running this pipeline to ensure you obtain the best results.
- `${SCPTS}/Basecalling.sh $RPATH $DPATH $FLOWCELL $KIT $ONT`
  - Calls the Basecalling script in the provided path to where the package scripts are stored. 
  - Uses `RPATH` as the input of the raw reads; `DPATH` as the output of the basecalling
  - Outputs will be saved in the directory called `Basecalling` within the experiment name you have provided 
#### Demultiplexing
- While guppy also contains a demultiplexer, qcat was found to be better at correct binning of barcodes
- However, qcat is now deprecated and **Oxford Nanopore now recommends the use of the guppy demultiplexer rather than qcat**
- The general set up for the Guppy basecaller (*guppy_barcoder*) is present in the `Demultiplexing.sh`. It can be run by uncommenting the section and changing `KIT2` to `KIT` on the `${SCPTS}/Demultiplexing.sh $DPATH $DEMULP $THREADS $KIT2 $ONT 2>&1 | tee $LOGFF/Demultiplexing/${EXPMT}.txt` command
#### Adapter Removal

#### Checks
### De Novo Assembly
#### Alignment
#### De novo assembly
#### Genome size estimation
#### Decontamination
### Correction and Polishing
#### Racon Polishing
#### Medaka Correction
#### Pilon Correction
### [Optional *but recommended*]: Apicoplast and Mitochondrial Separation and Circlarisation
#### Remove API and MIT
### RepeatMasking
#### RepeatMasking
#### De-chrimerisation
### Quality Assessment
#### QUAST
#### BUSCO
#### Pomoxis
#### Assembly-stats
### Ready for annotation
## FAQ
### I haven’t used barcodes. What do I do?
## Acknowledgements
## References
