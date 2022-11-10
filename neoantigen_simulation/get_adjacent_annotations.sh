#!/bin/bash

module load bedtools/2.27.1

intvl=$1
annot=$2
outname=$3

echo -e "$intvl"
echo -e "$annot"
echo -e "$outname"


# intersect to get the gtf
bedtools intersect -wb -a $intvl -b $annot > .temp.$outname.intersect.bed

# find the preceding and succeeding exons
a=0
while read i; do 
	iline=$(echo -e "$i" | sed 's/\t/\\\t/g' | sed 's/ /\\ /g')
	echo -e "$iline"
	grep -B 1 -A 1 -s "$iline" $annot | cut -f1,4- > .temp.$outname.$a.bed
	cat .temp.$outname.$a.bed
	a=$(( a + 1 ))
done < <(cut -f5- .temp.$outname.intersect.bed)
