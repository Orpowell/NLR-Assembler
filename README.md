# NLR-Assembler and Pipeline

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
    python3 main.py <command> ...

These commands are used at specific points in the NLR-Assembler pipeline (see below).

## Index

The index command generates an index of all barcoded reads where each read is assigned a colour (RGB value) according to its barcode. Reads with the same barcode are assigned the same colour in the index. Reads with no barcode matching the whitelist provided are assigned the colour black (RGB value : 0,0,0). Sequencing errors in barcodes are automatically corrected whilst generating the output. Note: index is designed to run using multiple cores (we used 12).

The final output is a csv file named 'barcode_index.csv'. 

Run the index command as follows:
 
    python3 main.py index --fastq input.fasta --adapters whitelist.txt --cores 12 --split 1000000

### parameters

parameter | argument | description|
|---|---|---|
--fastq | input.fasta | The raw barcoded reads in fastq format (Forward reads (R1) for 10x genomics data)
--adapters | whitelist.txt | A whitelist of all valid adapter sequences used by 10x Genomics available [here](https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt) in txt format.
--cores | 12 | The number of cores utilised as an integer (we recommend 12). If left blank all available cores will be used.
--split | 1000000 | The number of sequences per core (default is 1000000).

## Group

The group command generates an assembly by groupings contigs from a standard assembly that have been assembled from reads with similar sets of barcodes. A barcode profile is generated for each contig in the assembly from mapping of the reads to the assembly. The barcode profiles for each contig are weighted using Term-Frequency Inverse Document Frequency and used to generate a cosine similarity matrix comparing all contigs (1). Contigs with a high cosine similarity are grouped together and combined into a single contig, separated by a spacer of ambiguous nucleotides to produce the final assembly(2).

![Cosine Concept Figure](images/reduced_cosine_concept.pdf)

The final output is a fasta file named grouped_assemblies.fa

    python3 main.py group --samfile mapping.sam  --index barcode_index.csv --assembly assembly.fasta  --blast contig_bait.blastn

### Parameters

parameter | argument | description|
|---|---|---|
--samfile | mapping.sam | Raw reads mapped to the assembly with PCR duplicates removed in sam format.
--index | barcode_index.csv | The index file generated using the index command (see above).
--assembly | assembly.fasta | A "draft" assembly generated using de-barcoded reads in fasta format.
--blast | contig_bait.blastn | An alignment of RenSeq baits used to generate the raw reads to the generic assembly in BLAST6 format.

The specific details for generating each file are explained in the NLR-Assembler Pipeline section.

## Contig Coverage (Validation Only)



![NLR Coverage](images/concept_contig_coverage.pdf)

The final output is a csv file containing information on each contig group. In addition, overall statistics for grouped contigs are logged in the standard output.


    python3 main.py contig-coverage -assembly grouped_assemblies.fa -blast contig_coverage.blastn

### Parameters
parameter | argument | description|
|---|---|---|
--assembly | grouped_assemblies.fa | NLR-Assembler assembly in fasta format (output of the group command shown above)
--blast | contig_coverage.blastn | An alignment of the NLR-Assembler assembly to a reference genome in BLAST6 format.

The specific details for generating each file are explained in the NLR-Assembler Pipeline section.

## NLR Coverage (Validation Only)

    python3 main.py nlr-coverage --draft draft_nlr_coverage.blastn \
    --final final_nlr_coverage.blastn \
    --nlr nlr_sequences.fasta

Output info required

### Parameters

parameter | argument | description|
|---|---|---|
--draft | draft_nlr_coverage.blastn | An alignment of NLR sequences annotated from a reference genome to the draft assembly in BLAST6 format.
--final | final_nlr_coverage.blastn | An alignment of NLR sequences annotated from a reference genome to the final assembly in BLAST6 format.
--nlr | nlr_sequences.fasta | NLR sequences annotated from a reference genome using NLR-Annotator and converted to a single-line FASTA file.

The specific details for generating each file are explained in the NLR-Assembler Pipeline section.

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

### A *De Novo* Assembler of Your Choice

A widerange of commercial and open-source *de novo* assemblers can be used to produce an assembly from raw data for the pipeline. We found [ABySS](https://github.com/bcgsc/abyss) (v2.3.5) to be the most effective tool for our data.

### Seqkit

Available to download from https://bioinf.shenwei.me/seqkit/download/

### NLR-Annotator Version 2 (Validation only)

Availavble to download from https://github.com/steuernb/NLR-Annotator

### RenSeq Baits

The sequences of the RenSeq baits used for your specific experiments

### 10x Genomics Barcode Whitelist

A whitelist of all valid adapter sequences used by 10x Genomics for linked-read sequencing is available [here](https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt).

## Preprocessing

### Process Raw Data
Debarcoded and subsampled (~4% of data) as shown below:

    longranger basic --id=processed --fastqs=raw_data_directory

    gunzip processed/outs/barcoded.fastq.gz

    sed -n '1~200p;2~200p;3~200p;4~200p;5~200p;6~200p;7~200p;8~200p' processed/outs/barcoded.fastq > processed_reads.fastq

### *De Novo* Assembly of Processed reads
Generate a draft assembly using the processed reads and an assembler of your choice. For our data, we had the best success with [ABySS](https://github.com/bcgsc/abyss) and the following parameters. Then renumber the contigs in the draft assembly using [seqkit](https://bioinf.shenwei.me/seqkit/download/) and convert the file from a multi-line fasta file to a single line fasta file.

    abyss-pe name=assembly k=85 B=100G kc=3 H=3 v=-v pe='pea' pea=processed_reads.fastq lr='lra' lra=processed_reads.fastq

    seqkit replace -p '.+' -r '{nr}' assembly-contigs.fa -o renumbered_assembly-contigs.fa

    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < renumbered_assembly-contigs.fa > flat_renumbered_contigs.fa

    tail -n +2 flat_renumbered_contigs.fa > draft_assembly.fa

### Read Indexing

 Use the NLR-Assembler index command (see above) to generate an index the raw barcoded reads (Not the reads processed by Long Ranger Basic).

    python3 main.py index --fastq Raw_R1_reads.fasta --adapters whitelist.txt --cores 12 --split 1000000

### RenSeq Bait Alignment

Align RenSeq baits from your experiment to the draft assembly using [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata).

    blastn -query renseq_baits.fasta -subject draft_assembly.fasta -outfmt 6 -out bait_assembly_alignment.blastn

### Mapping

Map the (cleaned) raw reads used for *de novo* assembly against your draft assembly before removing PCR duplicates from the mappign. We recommend using [bwa](http://bio-bwa.sourceforge.net/) and [samtools](http://samtools.sourceforge.net/).

    bwa index draft_assembly.fasta

    bwa mem draft_assembly.fasta processed_reads.fastq > mapping.sam

    samtools view -Sh mapping.sam > filtered_mapping.sam

    samtools sort -n filtered_mapping.sam > mapping.sorted.bam

    samtools fixmate -m mapping.sorted.bam fixmate.bam

    samtools sort fixmate.bam > fixmate.sorted.bam

    samtools markdup -r -s fixmate.sorted.bam markdup.bam

    samtools view -Sh markdup.bam > markdup.sam


## Generating and Using the Final Assembly

Use the NLR-Assembler group command (see above) to generate the final assembly using the bait alignment, read_mapping, draft assembly and index generated during preprocessing.

    python3 main.py group --samfile mapping.sam \
    --index barcode_index.csv \
    --assembly draft_assembly.fasta \
    --blast bait_assembly_alignment.blastn

The final assembly can then be used in place of the standard wildtype assembly used in the (MutantHunter Pipeline)[https://github.com/steuernb/MutantHunter] to identify novel resistance gene candidates.

## Validation (requires a Reference Genome)

If a reference genome is available  NLR and contig coverage can be used to assess the quality of the final assembly using NLR-Assembler. Contig coverage determines the size of the region covered by a set of contigs grouped together by NLR-Assembler. NLR coverage determines the average percentage of annotated NLR sequences covered by the draft and final assemblies respectively.


### NLR Coverage
Lorem ipsum

    NLR-Annotator -i reference_genome.fa \
    -x mot.txt \
    -y store.txt \
    -o nlr_annotations.txt \
    -m nlr_annotations.motifs.bed \
    -f reference_genome.fa nlr_annotations.fa 0

    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <  nlr_annotations.fa > multiline_nlr_annotations.fa

    tail -n +2 multiline_nlr_annotations.fa > reformatted_nlr_annotations.fa 

    blastn -max_target_seqs 1 -query draft_assembly.fasta \
    -subject reformatted_nlr_annotations.fa \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out draft_assembly_nlr_coverage.blastn

    blastn -max_target_seqs 1 -query final_assembly.fasta \
    -subject reformatted_nlr_annotations.fa \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out final_assembly_nlr_coverage.blastn

    python3 main.py nlr-coverage -b draft_assembly_nlr_coverage.blastn \
    -c final_assembly_nlr_coverage.blastn \
    -n reformatted_nlr_annotations.fa

Results expalined


### Contig Coverage

Lorem Ipsum

    blastn -query $contigs -subject reference_genome.fasta \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
    -out query_coverage.blastn

    python3 main.py query-coverage -b query_coverage.blastn -c final_assembly.fasta

Results explained

## Tips & Issues 

## Acknowledgements

A massive thank you to Dr. Burkhard Steuernagel and Dr. Brande Wulff who co-supervised my dissertation project!

## Contact

If you wish to discuss NLR-Assembler, raise issues or collaborate with me please contact us [here](mailto:nlrassembler@oliverpowell.com)



