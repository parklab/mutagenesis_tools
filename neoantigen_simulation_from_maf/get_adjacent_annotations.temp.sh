#!/bin/bash

# module load bedtools/2.27.1

intvl=$1
annot=$2
outname=$3

echo -e "$intvl"
echo -e "$annot"
echo -e "$outname"

echo -e "$PWD"

# intersect to get the gtf
# bedtools intersect -wb -a intvl -b $annot > .temp.$outname.bed

# ls -l

# find the preceding and succeeding exons
# a=0
# while read line; do 
#	iline=$(echo -e "$i" | sed 's/\t/\\\t/g' | sed 's/ /\\ /g')
#	grep -B 1 -A 1 -s $i $exongtf | cut -f1,4- > .temp.$outname.$a.bed
#	cat .temp.$outname.$a.bed
#	a=$(( a + 1 ))
# done < <(cut -f4- .temp.$outname.bed)
