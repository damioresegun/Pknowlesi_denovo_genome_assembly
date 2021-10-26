# README
## Overview
This wiki is designed to work as a kind of manual which can hopefully take you through the major steps of the pipeline.
A pipeline to take input data from raw nanopore reads to being ready for manual or automated annotation. 
It is important to note that while this pipeline automates alot of commands, it cannot automate everything. 
User input is needed at various points to confirm which branch of the pipeline to follow as well as when to continue at certain 'break-points'
I'm hoping that this has been fairly prepped in a fairly straight-forward manner that is somewhat approachable
However, if you want to tweak the script, the vast majority is written in Bash with some secondary scripts in python.
Please be aware that the first portion of this script is tuned for barcoded ONT sequencing reads however, it is possible to use non-barcoded ONT sequenced read.
## Update on Pipeline
This repository mainly focuses on the generation of the de novo assemblies and their subsequent improvement and analyses. Further analyses can be found in: https://github.com/damioresegun/DNDS_analyses for duplication and subsitution ratios (Forked from Peter Thorpe). Additionally, analysis outputs can be found on the associated zenodo repository.
## Pipeline
![image](https://user-images.githubusercontent.com/33520829/112047928-cb046800-8b45-11eb-91a4-068118a73c4e.png)
## Table of Contents
  * [Installation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Installation)
    + [Processing](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Installation#processing)
    + [Conda](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Installation#conda)
    + [Other tools not on Conda](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Installation#other-tools-not-on-conda)
  * [Script and Pipeline Preparation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Script-and-Pipeline-Preparation)
  * [Data Preparation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Script-and-Pipeline-Preparation)
    + [Basecalling](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Data-Preparation#basecalling)
    + [Demultiplexing](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Data-Preparation#demultiplexing)
    + [Adapter Removal](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Data-Preparation#adapter-removal)
    + [Other Checks](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Data-Preparation#checks)
  * [De Novo Assembly](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/De-Novo-Assembly)
    + [Alignment](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/De-Novo-Assembly#alignment)
    + [De novo assembly](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/De-Novo-Assembly#de-novo-assembly)
  * [Genome size estimation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Genome-size-estimation)
    + [Known issue:](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Genome-size-estimation#known-issue)
  * [Decontamination](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Decontamination)
    + [Known Issues:](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Decontamination#known-issues)
  * [Correction and Polishing](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Correction-and-Polishing)
    + [Racon Polishing](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Correction-and-Polishing#racon-polishing)
    + [Medaka Correction](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Correction-and-Polishing#medaka-correction)
    + [Pilon Correction](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Correction-and-Polishing#pilon-correction)
      - [Known Issue](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Correction-and-Polishing#known-issue)
  * [Optional: Apicoplast and Mitochondrial Separation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/%5BOptional-but-recommended%5D:-Apicoplast-and-Mitochondrial-Separation)
    + [Remove API and MIT](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/%5BOptional-but-recommended%5D:-Apicoplast-and-Mitochondrial-Separation#remove-api-and-mit)
    + [Circularisation and Completion](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/%5BOptional-but-recommended%5D:-Apicoplast-and-Mitochondrial-Separation#circularisation-and-completion)
  * [RepeatMasking](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/RepeatMasking)
    + [Known issues](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/RepeatMasking#known-issues)
  * [De-chimerisation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/De-chimerisation)
  * [Quality Assessment](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Quality-Assessment)
    + [QUAST](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Quality-Assessment#quast)
    + [BUSCO](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Quality-Assessment#busco)
    + [Pomoxis](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Quality-Assessment#pomoxis)
    + [Assembly-stats](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Quality-Assessment#assembly-stats)
  * [Ready for annotation](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/Ready-for-annotation)
  * [Duplication and Substitution Analyses](https://github.com/damioresegun/DNDS_analyses)
  * [FAQ](https://github.com/damioresegun/Pknowlesi_denovo_genome_assembly/wiki/FAQ)
## Data Access
The completed genomes generated from this pipeline was uploaded to Zenodo and NCBI. 
- Zenodo repository: https://doi.org/10.5281/zenodo.5598264
- NCBI Accession Code:

## How to Cite
- Pre-publication: https://doi.org/10.5281/zenodo.5598264
- Post-publication:
## Acknowledgements
- Dr Peter Thorpe helped in various parts of the pipeline development; in particular the repeatmasking process; providing scripts to start off with and adapt
- Amir Szitenberg for his repeatmasking cookbook which is the backbone for the in-depth repeatmasking process utilised in this pipeline
- Dr Janet Cox-Singh helped in guiding the biological necessities required for this pipeline
## Contact
- Email: dro@st-andrews.ac.uk
- Alternative Email: damioresegun@gmail.com
- Twitter: @dami_oresegun
