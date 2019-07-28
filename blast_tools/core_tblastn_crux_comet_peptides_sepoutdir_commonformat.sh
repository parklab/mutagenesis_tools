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

sname=$(dirname $infile)
sname=$(dirname $sname)
filename=$(basename $sname)
filedir=$filename"_tblastn"
blastdbdir=$(dirname $blastdb)
blastdbname=$(basename $blastdb)


#oldfiledir=$(dirname $infile)
#filedir=$(basename $oldfiledir)
#filedir=$filedir"_tblastn"
#filename=$(basename $infile)
#filename=${filename%.txt}

mkdir -p $filedir

# Step 3: label each line based on the type of file being submitted
# takes format with
# scan number\tpeptide\tsequence\tsource(s)\tunrestricted number of score columns
# expects

#sname=$(dirname $infile)
#sname=$(dirname $sname)
#sname=$(basename $sname)
#

if [ $filetype == "barista" ]; then
    # format barista
    # awk 'BEGIN {FS=OFS="\t"} {print $5,$18,$20,$1,$2,$3,$4}' $infile
    # label
    awk -v fname=$filename 'BEGIN {FS="\t";OFS=";"} {print ">scan_"$5,"mz_"$7,"mas_"$8,fname}' $infile | paste - <(cut -f18 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
    awk -v fname=$filename 'BEGIN {FS=OFS="\t"} {print $5,$7,$8,fname,$20,$21}' $infile | paste <(grep -s ">" $filedir/$filename.peptides.fa | sed "s/>//g") - > $filedir/$filename.labels.txt
    # paste <(tr '\t' ';' $filedir/$filename.labels.temp) <(cut -f18 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
elif [ $filetype == "percolator" ]; then
    # format percolator
    awk -v fname=$filename 'BEGIN {FS="\t";OFS=";"} {print ">scan_"$2,"mz_"$4,"mas_"$5,fname}' $infile | paste - <(cut -f11 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
    awk -v fname=$filename 'BEGIN {FS=OFS="\t"} {print $2,$4,$5,fname,$11,$12}' $infile | paste <(grep -s ">" $filedir/$filename.peptides.fa | sed "s/>//g") - > $filedir/$filename.labels.txt
else
    # default: assume that we are formatting a comet file
    awk -v fname=$filename 'BEGIN {FS="\t";OFS=";"} {print ">scan_"$1,"mz_"$3,"mas_"$4,fname}' $infile | paste - <(cut -f14 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
    awk -v fname=$filename 'BEGIN {FS=OFS="\t"} {print $1,$3,$4,fname,$14,$15}' $infile | paste <(grep -s ">" $filedir/$filename.peptides.fa | sed "s/>//g") - > $filedir/$filename.labels.txt
    #awk -v fname=$filename 'BEGIN {FS=OFS="\t"} {print "scan_"$1,"mz_"$3,"mas_"$4,fname}' $infile | tr '\t' ';' > $filedir/$filename.labels.temp
    #paste $filedir/$filename.labels.temp <(cut -f14 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
    # old style
    #cut -f1,14,15,17,18 $infile | tr '\t' ',' | sed "s/^/>/g" | \
    #awk -v fnm="$filename" 'BEGIN {FS=OFS="\t"} {print $0,fnm}' - | tr '\t' ',' > $filedir/$filename.labels.temp
    ## Step 4: prepare FASTA
    #paste $filedir/$filename.labels.temp <(cut -f14 $infile) | sed "s/\t/\n/g" > $filedir/$filename.peptides.fa
fi

# paste $infile $filedir/$filename.labels.temp > $filedir/$filename.labels.txt

# Step 5: search with BLASTDB
BLASTDB=$BLASTDB:$blastdbdir

echo $blastdbname

tblastn -db $blastdbname \
-query $filedir/$filename.peptides.fa \
-out $filedir/$filename.tblastn_fast.txt \
-task tblastn-fast \
-max_intron_length 5000 \
-evalue $evalue \
-outfmt 7

# Step 6: recover the coordinates of the search and format as a BED file that we can then search against RNA-seq files
awk 'BEGIN {FS=OFS="\t"} {print $2,$9,$10,$1,$3,$4,$5,$6,$7,$8,$11,$12}' $filedir/$filename.tblastn_fast.txt > $filedir/$filename.tblastn_fast.bed

# Step 7: remove temp files
rm $filedir/$filename.labels.temp