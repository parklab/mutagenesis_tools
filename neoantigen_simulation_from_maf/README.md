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

## STEP 1: Setting up the config file

First, you will need to set up your config file to store various parameters for the run. You can find a sample config file in the repository

```
$ cat config_files/config.simulate_neoantigens.test.txt
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

## STEP 2: Setting up the input file

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

## STEP 3: Submitting the job

If you have a SLURM system, then the submission can be done like so
```
bash ./batch_calculate_neoantigens_from_maf.sh ./inputfile_template.txt ./config.simulate_neoantigens.test.txt
```

## Interlude

Note that Steps 4-5 are optional if you want to run `netMHCpan` using our scripts. 

## STEP 4: Once the jobs are done, produce the wt-vs-mutant peptide comparison file

Once your jobs are done, you will need to convert the results into a large table that summarizes the wild-type peptides and putative neoantigens on-site. 

You will need a file listing the mutant-vs-wildtype peptide comparisons. Again, these are just dummy readouts
```
$ cat completed.mutant_vs_wt_peptide_analyses.txt
sample1    ./sample1/sample1.mutant_vs_wt_peptide.txt
sample2    ./sample2/sample2.mutant_vs_wt_peptide.txt
```

We also expect outputs from OptiType (we used OptiType calls from TCGA). You will need the HLA assignments for a given sample like so
```
$ cat OptiTypeCallsHLA.tsv
A1,A2,B1,B2,C1,C2,Reads,Objective,aliquot_id
...
```

We provide a script to create this table. 

```
bash ./create_mutant_to_wt_with_hla_annot.sh <directory whwere you ran STEP 3> <HLA file from OptiType> <give the job a name>
```

If this script runs properly, then we will expect `$outname.all_samples_master_mutant_vs_wt.with_hla.txt` as an output file

## STEP 5: test MHC-peptide affinity changes

1. You will need to create a config for running pMHC. See the template provided in `config_files/config.pep_mhc_bind.test.txt`
    ```
    cat config_files/config.pep_mhc_bind.test.txt 
    corescript      ./pMHC_prediction/core_netmhcpan_mutant_vs_wt.revised.sh
    splitlines      10
    peplenmin       8
    neoantigen_mutant_wt_file       ./$outname.all_samples_master_mutant_vs_wt.with_hla.txt
    ```
2. Make sure you have `netMHCpan-4.0` installed in `./pMHC_prediction`
3. Run as below. This requires SLURM for parallelized jobs. 
```
bash ./pMHC_prediction/bash_netmhcpan_mutant_vs_wt.revised.sh $PWD/inputfile_template.txt <config> 
```

# Forthcoming modifications

Forthcoming modifications include:
1. A snakemake workflow instead of the current batch-core setup to simplify all tasks above
2. Support for simulating neoantigens directly from RNA-seq BAMs, splice junctions, etc. 
3. Support for other neoantigen prioritization tools
4. Suport for more diverse translation simulations (6-frame, alternative start codons, uORFs and unannotated ORFs, etc)
