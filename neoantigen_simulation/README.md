# Introduction

This folder contains scripts to simulate neoantigens from a VCF file. If an annotation GTF is supplied, then transcript-aware neoantigens are simulated. Works for SNV and indels; scripts for SV neoantigen simulation are being updated (i.e. using bedpe input).


# Requirements

* python 3.7.2 or greater
    * pysam
    * pybedtools 
* bedtools 2.27.1
* samtools 1.9
* netMHCpan 4.0

This pipeline is meant to run on computing clusters using SLURM -- it submits individual jobs as separate scripts to conduct neoantigen simulation from the VCF file.

# Instructions

## Setting up the config file

## Setting up the input file

## Submitting the job

# Forthcoming modifications

Forthcoming modifications include:
1. A snakemake workflow instead of the current batch-core setup
