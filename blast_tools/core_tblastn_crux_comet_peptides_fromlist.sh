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
filetype=$(grep -v "#" $inconfig | grep -s "type" | cut -f2)

# Step 2: extract the directory

# sname=$(dirname $infile)
# sname=$(dirname $sname)
# filename=$(basename $sname)
# filedir=$filename"_tblastn"
# blastdbdir=$(dirname $blastdb)
# blastdbname=$(basename $blastdb)

filename=$(basename $infile)
filename=${filename%.txt}
filedir=$filename"_tblastn"
blastdbdir=$(dirname $blastdb)
blastdbname=$(basename $blastdb)


#oldfiledir=$(dirname $infile)
#filedir=$(basename $oldfiledir)
#filedir=$filedir"_tblastn"
#filename=$(basename $infile)
#filename=${filename%.txt}

mkdir -p $filedir

# Step 3: Generate a FASTA
# this script assumes that you are submitting a two-column list for peptides:
paste <(for i in $(cut -f1 $infile); do echo -e ">$i"; done) <(cut -f2 $infile) | \
sed "s/\t/\n/g" > $filedir/$filename.peptides.fa

# paste $infile $filedir/$filename.labels.temp > $filedir/$filename.labels.txt

# Step 5: search with BLASTDB
BLASTDB=$BLASTDB:$blastdbdir

echo $blastdbname

tblastn -db $blastdbname \
-query $filedir/$filename.peptides.fa \
-out $filedir/$filename.tblastn_fast.txt \
-max_intron_length 5000 \
-evalue $evalue \
-outfmt 7


# -task tblastn-fast

# parsing the search
# Step 1: find all of the 0-hit peptides, locate their sequences, and recover their spectra
grep -B 4 -s "# 0 hits" $filedir/$filename.tblastn_fast.txt | grep -s "Query" | sed "s/: /\t/g" | cut -f2 | sort -k1,1 | uniq > $filedir/$filename.zero_hit.temp
# for i in $(cut )
# Step 2: find all of the peptides with no 100-pct search 

# Step 6: recover the coordinates of the search and format as a BED file that we can then search against RNA-seq files
# awk 'BEGIN {FS=OFS="\t"} {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12}' $filedir/$filename.tblastn_fast.txt > $filedir/$filename.tblastn_fast.bed

# Step 7: remove temp files
# rm $filedir/$filename.labels.temp