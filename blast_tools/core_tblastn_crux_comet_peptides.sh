#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-12:00
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o ./runlogs/out/%j.core_tblastn_crux.out.txt #stdout file
#SBATCH -e ./runlogs/err/%j.core_tblastn_crux.err.txt #stderr file

# Step 0: Module loads and variable loads
module load blast/2.6.0+

infilelist=$1
inconfig=$2

# Step 1: get the file we need and the BLAST database against we are searching the peptides
taskid=$SLURM_ARRAY_TASK_ID
jobid=$SLURM_ARRAY_JOB_ID
jobtask=$SLURM_JOB_ID

# what is the link?
infile=$(sed -n $taskid'p' $infilelist)
blastdb=$(grep -v "#" $inconfig | grep -s "blastdb" | cut -f2)
evalue=$(grep -v "#" $inconfig | grep -s "evalue" | cut -f2)

# Step 2: extract the directory
filedir=$(dirname $infile)
filename=$(basename $infile)
filename=${filename%.txt}
blastdbdir=$(dirname $blastdb)
blastdbname=$(basename $blastdb)

# Step 3: label each line
cut -f1,14,15,17,18 $infile | tr '\t' ',' | sed "s/^/>/g" | \
awk -v fnm="$filename" 'BEGIN {FS=OFS="\t"} {print $0,fnm}' - | tr '\t' ',' > $filedir/$filename.labels.temp
paste $infile $filedir/$filename.labels.temp > $filedir/$filename.labels.txt

# Step 4: prepare FASTA
paste $filedir/$filename.labels.temp <(cut -f14 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa

# Step 5: search with BLASTDB
BLASTDB=$BLASTDB:$blastdbdir

echo $blastdbname

tblastn -db $blastdbname \
-query $filedir/$filename.peptides.fa \
-out $filedir/$filename.tblastn_fast.txt \
-task tblastn-fast \
-max_intron_length 5000 \
-evalue $evalue \
-outfmt 6

# Step 6: recover the coordinates of the search and format as a BED file that we can then search against RNA-seq files
awk 'BEGIN {FS=OFS="\t"} {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12}' $filedir/$filename.tblastn_fast.txt > $filedir/$filename.tblastn_fast.bed

# Step 7: remove temp files
rm $filedir/$filename.labels.temp