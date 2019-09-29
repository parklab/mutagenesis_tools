#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o %j.generate_fusion_genes_modified.txt #stdout file
#SBATCH -e %j.generate_fusion_genes_modified.txt #stderr file


# inputs
infile=$1
config=$2
jobtask=$3
# bkptname=$5
scriptdir=$4




python /n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/tile_peptides_from_sequence_dev.py $jobtask.$gene.seq_around_lesion.$lesioname.fa