#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o ./out.%j.matchmutantwt.txt #stdout file
#SBATCH -e ./err.%j.matchmutantwt.txt #stderr file

#"""
#introduction:
#improvements for the core version
#1. submit properly-formatted files. Too much time spent on formatting within core
#2. each core job will run on one sample.
#3. Should add a config option to ask whether .temp files will be kept or not
#
#Steps
#1. bedtools intersect the mutant lesion and wt lesions
#2. using awk statements (NO WHILE LOOP NEEDED), filter for lines where (a) peptide length matches (see https://stackoverflow.com/questions/4648851/check-field-length-using-awk) and (b) hla allele name matches
#3. on the lines filtered in #2, compute the differences in the peptide strings
#4. paste desired columns to get desired output
#Move all temp files into storage. 
#"""


# package loads
module load \
bedtools/2.27.1 \
samtools/1.9
# set path of netMHCpan
PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/pMHC_prediction/netMHCpan-4.0/

# get files
infile=$1
config=$2
scriptdir=$3 # directory where the script is running

# get the sample name
taskid=$SLURM_ARRAY_TASK_ID
jobid=$SLURV_ARRAY_JOB_ID
jobtask=$SLURM_JOB_ID
insample=$(sed -n $taskid'p' $infile | cut -f1)

# get the names of the wt and mutant lesion files
mutantfile=$(grep -v "#" $config | grep -P "^mutantfile\t" | cut -f2)
wtfile=$(grep -v "#" $config | grep -P "^wtfile\t" | cut -f2)

# isolate the intervals that match the sample. store as a BED file
grep -P "$insample" $mutantfile | cut -f1-5,12,15-16 | sort -k1,1 -k2,2n -k3,3n > $insample.mutant_sig.bed.temp
grep -P "$insample" $wtfile | cut -f1-5,12,15-16 | sort -k1,1 -k2,2n -k3,3n > $insample.wt.bed.temp

# map lesions and store only those lines where the peptides are the same length
bedtools intersect -wb -sorted -a $insample.mutant_sig.bed.temp -b $insample.wt.bed.temp | \
awk 'BEGIN {FS=OFS="\t"} {if (($5 == $13) && length($6) == length($14)) print $0}' > $insample.mutant_sig_mapped_wt.bed.temp

# now, compare the string differences
cut -f6,14 $insample.mutant_sig_mapped_wt.bed.temp > $insample.mutant_sig_mapped_wt.peptides.temp

while read line; do
    # get peptides
    mutpep=$(echo $line | tr ' ' '\t' | cut -f1)
    wtpep=$(echo $line | tr ' ' '\t' | cut -f2)
    # compute differences
    diffresid=$(diff <(fold -w1 <<< "$mutpep") <(fold -w1 <<< "$wtpep") | awk '/[<>]/{printf $2}')
    echo $diffresid
done < $insample.mutant_sig_mapped_wt.peptides.temp > $insample.peptide_diff.temp

# header
echo -e "#chr\tstart\tend\tsample\thla\tmutpep\twtpep\tmutaff\twtaff\tmutpct\twtpct\tmutstat\twtstat\tdiff" > $insample.header.temp
# paste results
paste $insample.mutant_sig_mapped_wt.bed.temp $insample.peptide_diff.temp | \
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$14,$7,$15,$8,$16,$NF}' | \
cat $insample.header.temp - > $insample.mutant_sig_mapped_wt_pep_diff.bed

# mark strong binders
awk 'BEGIN {FS=OFS="\t"} {if ($10 <= 0.5) print $0}' $insample.mutant_sig_mapped_wt_pep_diff.bed > $insample.strong_binder_mutant_sig_mapped_wt_pep_diff.bed
awk 'BEGIN {FS=OFS="\t"} {if ($10 <= 2.0 && $10 > 0.5) print $0}' $insample.mutant_sig_mapped_wt_pep_diff.bed > $insample.weak_binder_mutant_sig_mapped_wt_pep_diff.bed

# move temp files
mkdir -p tempdir
mv $insample.*.temp tempdir


#while read line; do
#	coord=$(echo $line | tr ' ' '\t' | cut -f1-3)
#	mutantpep=$(echo $line | tr ' ' '\t' | cut -f6)
#	wtpep=$(echo $line | tr ' ' '\t' | cut -f13)
#	muthla=$(echo $line | tr ' ' '\t' | cut -f5)
#	wthla=$(echo $line | tr ' ' '\t' | cut -f12)
#	mutsample=$(echo $line | tr ' ' '\t' | cut -f4)
#	wtsample=$(echo $line | tr ' ' '\t' | cut -f11)
#	muraff=$(echo $line | tr ' ' '\t' | cut -f7)
#	wtaff=$(echo $line | tr ' ' '\t' | cut -f14)
#	mutpeplen=${#mutantpep}
#	wtpeplen=${#wtpep}
#	# get string differences
#	diffresid=$(diff <(fold -w1 <<< "$mutantpep") <(fold -w1 <<< "$wtpep") | awk '/[<>]/{printf $2}')
#	#
#	# print out entries if samples match, HLAs match, and peptide lengths are similar
#	if [[ $mutsample =~ $wtsample ]] || [[ $wtsample =~ $mutsample ]]; then
#		if [ $muthla == $wthla ] && [ $mutpeplen -eq $wtpeplen ]; then 
#			echo -e "$coord\t$mutantpep\t$wtpep\t$muthla\t$mutsample\t$wtsample\t$muraff\t$wtaff\t$diffresid"
#		fi
#	fi
#	counter=$((counter + 1))
#	a=$(expr $counter % 1000)
#	if [[ $a -eq 0 ]]; then
#		echo -e "$counter lines complete" 1>&2 # print to stderr
#	fi
#done < $infile >> $outfile
#
## infile=sbpeptide_mutant_matching_wt_lesion.txt
#
#
## 
