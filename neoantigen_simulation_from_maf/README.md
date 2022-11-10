# Introduction

This folder contains scripts to simulate neoantigens from a BED file. If an annotation GTF is supplied, then transcript-aware neoantigens are simulated. Works for SNV and indels; scripts for SV neoantigen simulation are being updated (i.e. using bedpe input).


# Requirements

## Software
* python 3.7.2 or greater
    * pysam
    * pybedtools 
* bedtools 2.27.1
* samtools 1.9
* netMHCpan 4.0

## Inputs

We require an input of variants in BED format. See below (`Setting up the input file`) for more details. 

This pipeline is meant to run on computing clusters using SLURM -- it submits individual jobs as separate scripts to conduct neoantigen simulation from the VCF file.

# Instructions

You can clone this repository 

## Setting up the config file

First, you will need to set up your config file to store various parameters for the run. You can find a sample config file in the repository

```
$ cat config.simulate_neoantigens.test.txt
splitlines      50
exongtf ./Homo_sapiens.GRCh37.75.exon.gtf
genegtf ./Homo_sapiens.GRCh37.75.gene.gtf
genomefa        ./genome.fa
sloplen 36
corescript      ./core_calculate_neoantigens_from_maf.sh
```

1. `splitlines` will break up the input file of VCFs (see below) into batches, each `splitlines` long
2. `exongtf` is the GTF of annotated exons
3. `genegtf` is the GTF of annotated genes
4. `genomefa` is the FASTA of genomes
5. `sloplen` gives the # of nucleotides in each direction of the variant to capture for translation
6. `corescript` this is the core script that will be run on eacf VCF. Default is the `core_calculate_neoantigens_from_maf.sh` script -- make sure that the path to this script is appropriate

## Setting up the input file

The input file should be a 2-column file formatted like so (note that `sample1_filename.bed` and `sample2_filename.bed` are just dummy files put here for the purposes of demonstration)

```
$ cat inputfile_template.txt
sample1 ./sample1_filename.bed
sample2 ./sample2_filename.bed
$
$ cat ./sample1_filename.bed
1   1000000 1000000 C   T   + Test_snv   SNP
2   1000000 1000000 -   TATA  + Test_insertion   SNP
```

the columns are 
1. chromosome
2. start
3. end (note that we expect 1-based coordinates)
4. ref allele
5. alt allele
6. Strand (just keep this as + for now)
7. Description of mutation (an arbitrary description can be used here)
8. Type of mutation (an arbitrary type can be used here, though we expect SNP, DEL, and INS here for SNVs or SNPs, deletions, and insertions respectively)




## Submitting the job

## Once the jobs are done, produce the wt-vs-mutant peptide comparison file

## test MHC-peptide affinity changes

# Forthcoming modifications

Forthcoming modifications include:
1. A snakemake workflow instead of the current batch-core setup
