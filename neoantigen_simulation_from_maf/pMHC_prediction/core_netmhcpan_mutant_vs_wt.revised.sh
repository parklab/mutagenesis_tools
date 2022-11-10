#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-12:00
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --mem=24G
#SBATCH -o ./runlogs/%j.core_netmhcpan.out.txt #stdout file
#SBATCH -e ./runlogs/%j.core_netmhcpan.err.txt #stderr file



# #SBATCH -p park
# #SBATCH -A park_contrib



# path to netmhcpan
PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/pMHC_prediction/netMHCpan-4.0/

# inputs
infile=$1
config=$2
scriptdir=$3

# get the files to be run
taskid=$SLURM_ARRAY_TASK_ID
i=$(sed -n $taskid'p' $infile | cut -f1)

# get the master wt-vs-mutant file
neoafile=$(grep -P "neoantigen_mutant_wt_file\t" $config | cut -f2)
# length filter
peptidelenmin=$(grep -P "peplenmin\t" $config | cut -f2)
# note that the min has to be
if [[ $peptidelenmin -lt 8 ]]; then
    echp -e "minimum peptide length has to be 8 AA. Setting to 8"
    peptidelenmin=8
fi

# first, get the relevant entries in the neoafile for the sample
grep -P "^$i\t" $neoafile > $i.mutant_wt_neoa_pairs.txt

# extract hla types
hlainput=$(cut -f2-7 $i.mutant_wt_neoa_pairs.txt | uniq | sed "s/\t/,HLA-/g" | sed "s/\*//g" | sed "s/^/HLA-/g")

# extract the mutant and wt neoantigens -- only extract the pairs with comparable affinities
# awk -v len=$peptidelenmin 'BEGIN {FS=OFS="\t"} {if (length($11) >= len || length($12) >= len) print "pep_"NR,$0}' $i.mutant_wt_neoa_pairs.txt > $i.neoas_pairs.txt
awk -v len=$peptidelenmin 'BEGIN {FS=OFS="\t"} {if (length($11) >= len || length($12) >= len) print "pep_"NR,$0}' $i.mutant_wt_neoa_pairs.txt > $i.neoas_to_analyze.txt
awk 'BEGIN {FS=OFS="\t"} {print $1,$12,$13}' $i.neoas_to_analyze.txt | sed "s/\t\t/\t-\t/g" | sed "s/\t$/\t-/g" > $i.mutant_wt_pairs.txt.temp
cut -f2 $i.mutant_wt_pairs.txt.temp | grep -v "X" | grep -v -P "^$" | awk -v len=$peptidelenmin 'BEGIN {FS=OFS="\t"} {if (length($1) >= len) print $0}' > $i.mutant.txt.temp
cut -f3 $i.mutant_wt_pairs.txt.temp | grep -v "X" | grep -v -P "^$" | awk -v len=$peptidelenmin 'BEGIN {FS=OFS="\t"} {if (length($1) >= len) print $0}' > $i.wt.txt.temp

# run netMHCpan
netMHCpan -a $hlainput -f $i.mutant.txt.temp -s 1 -inptype 1 -BA 1 -v > $i.netmhcpan.mutant_neoas.out.nonX.txt &
netMHCpan -a $hlainput -f $i.wt.txt.temp -s 1 -inptype 1 -BA 1 -v > $i.netmhcpan.wt_neoas.out.nonX.txt &
wait

# now to condense the input into a table

grep -P "  PEPLIST" $i.netmhcpan.mutant_neoas.out.nonX.txt | \
sed "s/<=//g" | sed 's/  */\t/g' | cut -f2- | \
cat <(echo -e "Pos\tHLA\tPeptide\tCore\tOf\tGp\tGl\tIp\tIl\tIcore\tIdentity\tScore\tAff(nM)\t%Rank\tBindLevel") - > $i.netmhcpan.mutant_neoas.out.nonX.table.txt.temp &
grep -P "  PEPLIST" $i.netmhcpan.wt_neoas.out.nonX.txt | \
sed "s/<=//g" | sed 's/  */\t/g' | cut -f2- | \
cat <(echo -e "Pos\tHLA\tPeptide\tCore\tOf\tGp\tGl\tIp\tIl\tIcore\tIdentity\tScore\tAff(nM)\t%Rank\tBindLevel") - > $i.netmhcpan.wt_neoas.out.nonX.table.txt.temp &
wait

cut -f2,3,12- $i.netmhcpan.mutant_neoas.out.nonX.table.txt.temp | awk 'BEGIN {FS=OFS="\t"} {if (NF < 6) print $0,"NA"; else print $0}' > $i.mutant.pred.table.temp
cut -f2,3,12- $i.netmhcpan.wt_neoas.out.nonX.table.txt.temp | awk 'BEGIN {FS=OFS="\t"} {if (NF < 6) print $0,"NA"; else print $0}'> $i.wt.pred.table.temp

# match the pairs to one another
join -a1 -t$'\t' -1 2 -2 2 <(sort -k2b,2 $i.mutant_wt_pairs.txt.temp) <(sort -k2b,2 $i.mutant.pred.table.temp) | \
sort -k3b,3 | \
join -a1 -t$'\t' -1 3 -2 2 - <(sort -k2b,2 $i.wt.pred.table.temp) | uniq | \
awk 'BEGIN {FS=OFS="\t"} {if ($4 == $9) print $3,$2,$1,$4,$5,$6,$7,$8,$10,$11,$12,$13}' | \
cat <(echo -e "pairid\tmutant_neoa\twt_neoa\thla\tmutant_score\tmutant_affinity\tmutant_pct_rank\tmutant_binder\twt_score\twt_affinity\twt_pct_rank\twt_binder")  - > $i.mutant_wt_netmhcpan_pairs.txt

# now merge the above with the original table to yield the relevant output
cut -f1-4,6,8,10,12 $i.mutant_wt_netmhcpan_pairs.txt | sort -k1b,1 | \
join -a1 -t$'\t' <(sort -k1b,1 $i.neoas_to_analyze.txt) - | cut -f1-13,16- > $i.mutant_wt_pairs_affinity_measurements.txt.temp

# python script to calculate number of string differences
python $scriptdir/n_aa_diff.py $i.mutant_wt_pairs_affinity_measurements.txt.temp 12 13 > $i.mutant_wt_pairs_affinity_measurements.txt


mkdir -p $i.tempfiles
mv $i.*.temp $i.tempfiles
