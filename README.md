# NLR-Assembler

## Introduction



PLEASE NOTE THAT THIS SOFTWARE IS HIGHLY EXPERIMENTAL.

## Installation

Download [Python (3.10)](https://www.python.org/downloads/release/python-3109/) and the NLR-Assembler repository. All dependencies for NLR-Assembler are shown in the requirements.txt and can be installed using pip or anaconda. We recommend using an anaconda environment to generate the complete environment for NLR-Assembler.

The following commands can be used to setup an environment and install the requirements using anaconda:

    conda create -n NLR-Assembler python=3.10
    source activate NLR-Assembler
    pip3 install -r requirements.txt

## Running NLR-Assembler

Four commands are currently available using NLR-Assembler: index, group, contig-coverage and nlr-coverage. All can be running using the following generic command (once the conda environment has been activated):

    source activate NLR-Assembler
    python3 main.py <command> --input XXX --parameter XXX --output XXX

These commands are used at specific points in the NLR-Assembler pipeline (see below)

### Index

    python3 main.py index --fastq input.fasta --adapters whitelist.txt --cores 12 --split 1000000

#### Parameters
Parameter |Argument | Description|
|---|---|---|
--fastq | input.fasta | raw barcoded reads in fastq format
--adapters | whitelist.txt | A whitelist of all valid adapter sequences used by 10x Genomics available [here](https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt__;!!Nmw4Hv0!ycvz4SQfCxuE9_PRUpGVO--oCiODqHYPmg7Noau4s5gTnsYvrtpEN2IbEVUjFfmCfkU4dOyCC0VYXO0vyiahvQuV0XwT-YnZVKLBgA$) in txt format.










# The NLR-Assembler Pipeline

## Prerequisites

### Python v3.10

Available to download from https://www.python.org/downloads/release/python-3109/

### NLR-Assembler

Available to download from this repository (see documentation above).

### BWA

Available to download from http://bio-bwa.sourceforge.net/

### Samtools

Available to download from http://samtools.sourceforge.net/

### Blast+

Availavle to download from [NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

### Long Ranger Basic Pipeline

Available to download from [10x Genomics](https://support.10xgenomics.com/genome-exome/software/downloads/latest)

### Sed

A standard GNU package (available on most operating systems).

### A *De Novo* Assembler of Your Choice

A widerange of commercial and open-source *de novo* assemblers can be used to produce an assembly from raw data for the pipeline. We found [ABySS](https://github.com/bcgsc/abyss) (v2.3.5) to be the most effective tool for our data. 

### NLR-Annotator Version 2 (Validation only)

Availavble to download from https://github.com/steuernb/NLR-Annotator

## Preprocessing

### Read Indexing

### Genome Assembly

### RenSeq Bait Alignment

### Read mapping

## Improving RenSeq Assemblies

MutantHunter

## Validation

## Tips & Issues 

## Acknowledgements

## Contact


