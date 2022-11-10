#!/bin/bash

infile=$1
inhla=$2
outname=$3

if [[ $outname == '' ]]; then
    outname="output_file_name"
fi

echo -e "combining all mutant_vs_wt file names together into a single master file..."
echo -e "creating master file"
for i in $(cut -f1 $infile); do
        tempfile=$(grep -s $i $infile | cut -f2)
        paste <(grep -v "#" $tempfile | cut -f1-5 | sed "s/\t/~/g") \
        <(grep -v "#" $tempfile | rev | cut -f2 | rev) \
        <(grep -v "#" $tempfile | cut -f6-8) | \
        awk -v sample=$i 'BEGIN {FS=OFS="\t"} {print sample,$0}'
done > $outname.all_samples_master_mutant_vs_wt.txt

echo -e "setting up hla file"
for i in $(cut -f1 $infile); do
        grep -s $i $inhla | \
        head -1 | \
        sed "s/,/\t/g" | \
        cut -f1-6,9 | \
        awk -v sample=$i 'BEGIN {FS=OFS="\t"} {print sample,$0}'
done > $outname.samples_to_hla_type.txt


# for additional
for i in $(cut -f1 $indir/../20201115_all_completed.mutant_vs_wt_peptide_analyses.uniq.txt); do
        grep -s $i 20220506_samples_to_hla_table.combined.txt | \
        awk -v sample=$i 'BEGIN {FS=OFS="\t"} {print sample,$0}'
done > $outname.samples_to_hla_type.txt

echo -e "joining"
join -t$'\t' \
<(sort -k1b,1 outname.samples_to_hla_type.txt) \
<(sort -k1b,1 $outname.all_samples_master_mutant_vs_wt.txt) > $outname.all_samples_master_mutant_vs_wt.with_hla.txt