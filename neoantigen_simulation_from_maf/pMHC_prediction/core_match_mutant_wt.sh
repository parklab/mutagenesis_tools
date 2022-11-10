#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-24:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o ./out.%j.matchmutantwt.txt #stdout file
#SBATCH -e ./err.%j.matchmutantwt.txt #stderr file


infile=$1
outfile=$2

echo -e "#chr\tstart\tend\tmutpep\twtpep\thla\tmutsample\twtsample\taffinity_mutant\taffinity_wt\tbasediff" > $outfile

counter=0

"""
improvements for the core version
1. join statements to combine the two sets of samples 
"""

while read line; do
	coord=$(echo $line | tr ' ' '\t' | cut -f1-3)
	mutantpep=$(echo $line | tr ' ' '\t' | cut -f6)
	wtpep=$(echo $line | tr ' ' '\t' | cut -f13)
	muthla=$(echo $line | tr ' ' '\t' | cut -f5)
	wthla=$(echo $line | tr ' ' '\t' | cut -f12)
	mutsample=$(echo $line | tr ' ' '\t' | cut -f4)
	wtsample=$(echo $line | tr ' ' '\t' | cut -f11)
	muraff=$(echo $line | tr ' ' '\t' | cut -f7)
	wtaff=$(echo $line | tr ' ' '\t' | cut -f14)
	mutpeplen=${#mutantpep}
	wtpeplen=${#wtpep}
	# get string differences
	diffresid=$(diff <(fold -w1 <<< "$mutantpep") <(fold -w1 <<< "$wtpep") | awk '/[<>]/{printf $2}')
	#
	# print out entries if samples match, HLAs match, and peptide lengths are similar
	if [[ $mutsample =~ $wtsample ]] || [[ $wtsample =~ $mutsample ]]; then
		if [ $muthla == $wthla ] && [ $mutpeplen -eq $wtpeplen ]; then 
			echo -e "$coord\t$mutantpep\t$wtpep\t$muthla\t$mutsample\t$wtsample\t$muraff\t$wtaff\t$diffresid"
		fi
	fi
	counter=$((counter + 1))
	a=$(expr $counter % 1000)
	if [[ $a -eq 0 ]]; then
		echo -e "$counter lines complete" 1>&2 # print to stderr
	fi
done < $infile >> $outfile

# infile=sbpeptide_mutant_matching_wt_lesion.txt


# 
