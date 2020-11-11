#!/bin/bash

#SBATCH -p short
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=6G
#SBATCH -o ./runlogs/out/%j.core_calculate_neoantigens.out.txt #stdout file
#SBATCH -e ./runlogs/err/%j.core_calculate_neoantigens.err.txt #stderr file

## Step 0: inputs
infile=$1 # this should contain the list of BARISTA outputs
config=$2
scriptdir=$3

# get the SRR accession number to be downloaded
taskid=$SLURM_ARRAY_TASK_ID
outname=$(sed -n $taskid'p' $infile | cut -f1)
intervals=$(sed -n $taskid'p' $infile | cut -f2)

# extract from Config:
exongtf=$(grep -v "#" $config | grep -P "^exon_gtf\t" | cut -f2)
genomefa=$(grep -v "#" $config | grep -P "^genome_fa\t" | cut -f2)
sloplen=$(grep -v "#" $config | grep -P "^sloplen\t" | cut -f2)

mkdir -p $outname
cd $outname

echo -e "sample\t$outname\n\
intervals\t$intervals\n\
exongtf\t$exongtf\n\
genomefa\t$genomefa\n\
sloplen\t$sloplen\n\
####" > run_params.txt

#
python $scriptdir/get_neoantigens_at_mutations.py