#!/bin/bash

module load bedtools/2.27.1

intvl=$1
annot=$2
outname=$3
getalladjacent=$4


echo -e "$intvl"
echo -e "$annot"
echo -e "$outname"
echo -e "$getalladjacent\tgetalladjacent"
if [[ $getalladjacent == "" ]]; then
	echo -e "no value for getalladjacent specified"
	getalladjacent=0
fi

if [[ $getalladjacent -gt 0 ]]; then
	getalladjacent=1
	echo -e "getting all exons"
else
	getalladjacent=0
	echo -e "getting only the immediately adjacent exons"
fi


# intersect to get the gtf
# bedtools intersect -wb -a $intvl -b $annot > .temp.$outname.intersect.bed

# find the preceding and succeeding exons
a=0
while read i; do
	iline=$(echo -e "$i" | sed 's/\t/\\\t/g' | sed 's/ /\\ /g')
	transcriptid=$(echo -e "$i" | cut -f9 | cut -d' ' -f4 | sort -k1,1 | uniq)
	echo -e "$transcriptid"
	echo -e "$iline"
	if [[ $getalladjacent -eq 0 ]]; then
		grep -B 1 -A 1 -s "$iline" $annot | cut -f1,4- | grep -s "$transcriptid" > .temp.$outname.$a.bed
	else
		grep -B 10000 -A 10000 -s "$iline" $annot | cut -f1,4- | grep -s "$transcriptid" > .temp.$outname.$a.bed
	fi
	cat .temp.$outname.$a.bed
	a=$(( a + 1 ))
done < <(cut -f5- $intvl)
# done < <(cut -f5- .temp.$outname.intersect.bed)
