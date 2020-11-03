#!/bin/bash

# load the conda environment with alignment tools
source /homes/vvv776/anaconda3/etc/profile.d/conda.sh
conda activate /homes/vvv776/anaconda3/envs/alignment
# loads bwa, fermi, stringtie, minimap2, and samtools

# inputs
inmaffiles=$1
config=$2
scriptdir=$3

taskid=$LSB_JOBINDEX
filename=$(sed -n $taskid'p' $inmaffiles | cut -f2) # the starting BAM file
jobtask=$(sed -n $taskid'p' $inmaffiles | cut -f1) # name of the outputs

echo -e "$inputfile" > $outname.input_params.txt
echo -e "$config" >> $outname.input_params.txt
echo -e "$scriptdir" >> $outname.input_params.txt

# GTFs from the config file
genegtf=$(grep -v "#" $config | grep -s "genegtf" | cut -f2)
exongtf=$(grep -v "#" $config | grep -s "exongtf" | cut -f2)
# get the coordinate system
coordsys=$(grep -v "#" $config | grep -s "coord_sys" | cut -f2)
# do we need to shuffle genes to look for chimeras?
shuffle=$(grep -v "#" $config | grep -s "shuffle" | cut -f2)

# output directory?
mkdir -p $jobtask
cd $jobtask
mkdir -p temp_files # tempfiles

cut -f5-7,9-10,35-36 $filename | sed "s/\t/_/g" > $jobtask.bkpt.txt
nbkpts=$(grep -c "^" $jobtask.bkpt.txt)
echo -e "$nbkpts variants found"
filenamebase=$(basename $filename)

if [ $coordsys == "1" ]; then
    # awk 'BEGIN {FS=OFS="\t"} {print $5,$6-2,$6-1,$5,$7,$7+1}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$6,$5,$7,$7+1}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
    # $8,$9,$10,$11,$12,$13,$14}' $filename
    # awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$7,$8,$9,$10,$11,$12,$13,$14}' $filename | tail -n +2 > bkpt.internal.$jobtask.temp
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$7}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.internal.$jobtask.temp
    # https://stackoverflow.com/questions/1602035/how-to-print-third-column-to-last-column for faster iteration
    # paste <(cut -f5-7 $filename) <(cut -f35 $filename) <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
elif [ $coordsys == "0" ]; then
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$6,$5,$7,$7+1}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
    #awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$6,$5,$7,$7+1,$8,$9,$10,$11,$12,$13,$14}' $filename | paste - $jobtask.bkpt.txt | tail -n +2 > bkpt.$jobtask.temp
#     awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' $filename | tail -n +2 > bkpt.internal.$jobtask.temp
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$7}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.internal.$jobtask.temp
    # paste <(cut -f5-7 $filename) <(cut -f35 $filename) <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
else
    echo -e "INVALID COORDINATE SYSTEM SPECIFIED. Please specify \"0\" for zero-based coordinates\n
    (eg: X   0   1 in zero-based coordinates is the first position of chromosome X)\n
    or  \"1\" for one-based coordinates\n
    (eg: X   1   1 in one-based coordinates the first position of chromosome X)\n"
fi

echo -e "identified locations of variants"


###PART 1: GENERATE FUSION TRANSCRIPT GTF

## convert PART 1 into a script! 

# remove any existing version...
[ -e bkpt.$jobtask.fusion_transcripts.gtf ] && rm bkpt.$jobtask.fusion_transcripts.gtf

# run collect_exons_around_break.sh to execute items 1-5
if [ $shuffle == "FALSE" ]; then
    bash $scriptdir/collect_exons_around_snvindel.ebi_version.sh \
    ./bkpt.internal.$jobtask.temp \
    $config \
    $jobtask \
    $scriptdir
else
    bash $scriptdir/collect_exons_around_break_shuffle.ebi_version.sh \
    ./bkpt.$jobtask.temp \
    ./bkpt.internal.$jobtask.temp \
    $config \
    $jobtask \
    $scriptdir    
fi

# collect the neoantigens
for i in $(find ./ -name "*unique_neoas_mutant.txt" | sort -k1,1); do
    neoafile=$i
    for j in $(cut -f1 $neoafile); do
        echo -e "$neoafile\t$j\t$filenamebase" | sed "s/_/\t/g"
    done
done > $jobtask.mutant_neoantigens_per_gene.txt

# collect the WT peptides
for i in $(find ./ -name "*lost_neoas_mutant.txt" | sort -k1,1); do
    neoafile=$i
    for j in $(cut -f1 $neoafile); do
        echo -e "$neoafile\t$j\t$filenamebase" | sed "s/_/\t/g"
    done
done > $jobtask.unique_wt_neoantigens_per_gene.txt

# match the two
for i in $(find ./ -name "*unique_neoas_mutant.txt" | sort -k1,1); do
    a=$(basename $i)
    a=${a%.unique*}
    wtfile=$(find ./ -name $a"*lost_neoas_mutant.txt")
    echo -e "$i\t$wtfile"
    paste $i $wtfile | awk -v name=$a 'BEGIN {FS=OFS="\t"} {print name,$0}'
done > $jobtask.unique_mutant_vs_wt_neoantigens_per_gene.txt

# end