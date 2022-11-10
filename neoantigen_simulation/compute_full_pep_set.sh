#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o %j.out.txt
#SBATCH -e %j.err.txt

infile=$1
headerfile=$2

newheader="#"
for i in $(cut -f1 $headerfile); do newheader=$newheader""$i"\t"; done

for i in $(cut -f1 $infile); do
	outfile=$i.matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.txt
	simpfile=$i.matched_mutant_wt_pep_pair.simple_affinity_comparison.txt
	netmhcpanpair=$(grep -P $i"\t" $infile | cut -f3)
	affpair=$(grep -P $i"\t" $infile | cut -f2)
	join -t$'\t' <(sort -k1b,1 $netmhcpanpair) <(sort -k1b,1 $affpair) | \
	awk '$4 == $26 && $6 == $27 && $10 == $29' | \
	uniq | \
	cat <(echo -e "$newheader") - | \
	cut -f1-22 > $outfile 
	awk 'BEGIN {FS=OFS="\t"} {print $21,$20,$4,$2,$3,$6,$7,$10,$11,$8,$12}' $outfile > $simpfile
done

cat *.matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.txt | sort -k1,1n | \
uniq > all_files_matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.txt

cat *.matched_mutant_wt_pep_pair.simple_affinity_comparison.txt | sort -k1,1n | \
uniq > all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.txt


	
cat header_for.all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.txt all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.txt > all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.reheader.txt
cat header_for.all_files_matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.txt all_files_matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.txt > all_files_matched_mutant_wt_pep_pair.netmhcpan_result_affinity_comparison.reheader.txt
