#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=6G
#SBATCH -o ./runlogs/out/out.%j.generate_fusion_genes_modified.txt #stdout file
#SBATCH -e ./runlogs/err/err.%j.generate_fusion_genes_modified.txt #stderr file


# package loads
module load \
bedtools/2.27.1 \
samtools/1.9
# set path of netMHCpan
PATH=$PATH:/n/data1/hms/dbmi/park/vinay/pipelines/pMHC_prediction/netMHCpan-4.0/

# inputs
infile=$1 # contains the file with sample names, their HLA alleles, and their corresponding .maf files
config=$2 # config file
scriptdir=$3 # directory where the script is running

# extract the neoantigen BED file from the config file
neoantigenbed=$(grep -v "#" $config | grep -s "neoantigens" | cut -f2)
colnumneo=$(grep -v "#" $config | grep -s "colnumber_neoantigens" | cut -f2)
truncpep=$(grep -v "#" $config | grep -s "truncate_peptide" | cut -f2)
strongbindlim=$(grep -v "#" $config | grep -s "strong_bind_pct_lim" | cut -f2)
weakbindlim=$(grep -v "#" $config | grep -s "weak_bind_pct_lim" | cut -f2)
colnmhcp=2
echo -e "obtained parameters
neoantigen: $neoantigenbed
column number of neoantigens in bed: $colnumneo
whether to discard truncated peptides: $truncpep
percentage limit for strong binders: $strongbindlim
percentage limit for weak binders: $weakbindlim
peptide column for MHC: $colnmhcp" 1>&2

# input sample
taskid=$SLURM_ARRAY_TASK_ID
jobid=$SLURV_ARRAY_JOB_ID
jobtask=$SLURM_JOB_ID
insample=$(sed -n $taskid'p' $infile | cut -f1)
inhla=$(sed -n $taskid'p' $infile | cut -f2 | cut -d',' -f1-6 | sed "s/,/,HLA-/g" | sed "s/^/HLA-/g" | sed "s/\*//g")
# print out the inputs
echo -e "$jobtask\t$insample\t$inhla" >> jobtask_insample_inhla.txt

# output directory?
mkdir -p $insample"_out"
cd $insample"_out"

# find the neoantigens of interest. If we need truncated neoantigens, then produce them
grep -P "$insample" $neoantigenbed > $insample.neoantigens.bed.temp

if [[ $truncpep -eq 1 ]]; then
    echo -e "obtaining non-truncated neoantigens"
    cat <(grep -v "\*" $insample.neoantigens.bed.temp) \
    <(grep -P "\*\t" $insample.neoantigens.bed.temp | grep -v -P "\*[A-Z]+\*") | \
    sort -k1,1 -k2,2n | uniq | sort -k$colnumneo"b",$colnumneo > $insample.neoantigens.bed
else
    sort -k$colnumneo"b",$colnumneo $insample.neoantigens.bed.temp > $insample.neoantigens.bed
fi
echo -e "obtained peptides" 1>&2

# check if we need to complete a run.
nlines=$(grep -c "^" $insample.neoantigens.bed)
# if we have no lines, then we exit the script
if [[ $nlines -eq 0 ]]; then
    echo -e "$nlines desired neoantigens found. DONE."
    echo -e "Exiting script" 1>&2
    exit 1
fi

# if we do have lines, then we continue
echo -e "$nlines neoantigens found for $insample. Predicting affinities to $inhla"

# create input
cut -f$colnumneo $insample.neoantigens.bed > $insample.peptides.input.txt

# run netMHCpan
netMHCpan \
-a $inhla \
-f $insample.peptides.input.txt \
-rth $strongbindlim \
-rlt $weakbindlim \
-s 1 \
-inptype 1 \
-BA 1 > $insample.netmhcpan_out.txt # we provide binding affinities
echo -e "netMHCpan completed" 1>&2

# count strong and weak binders
# grep -s "<= SB" $insample.netmhcpan_out.txt | sed "s/<= //g" | sed 's/  */\t/g' | sed "s/^\t//g" | sort -k14,14n | uniq > $insample.strong_binders.txt
# grep -s "<= WB" $insample.netmhcpan_out.txt | sed "s/<= //g" | sed 's/  */\t/g' | sed "s/^\t//g" | sort -k14,14n | uniq > $insample.weak_binders.txt
grep -s "PEPLIST" $insample.netmhcpan_out.txt | sed "s/<= //g" | sed 's/  */\t/g' | sed "s/^\t//g" | sort -k14,14n | uniq > $insample.all_peptides.temp
grep -s "Pos" $insample.netmhcpan_out.txt | sed 's/  */\t/g' | sed "s/^\t//g" | sort -k14,14n | uniq > $insample.mhcout_header.txt
# cat $insample.mhcout_header.txt $insample.strong_binders.txt $insample.weak_binders.txt | cut -f2- | sort -k$colnmhcp"b",$colnmhcp > $insample.sig_binders.txt
cat $insample.mhcout_header.txt $insample.all_peptides.temp | cut -f2-  > $insample.all_peptides.txt
cat $insample.mhcout_header.txt <(grep -P "[SW]B" $insample.all_peptides.temp) | cut -f2- > $insample.sig_binders.txt

echo -e "obtained strong and weak binders" 1>&2

# count entries
# nstrongbind=$(grep -c "^" $insample.strong_binders.txt)
nstrongbind=$(grep -c -P "\tSB" $insample.sig_binders.txt)
# nweakbind=$(grep -c "^" $insample.weak_binders.txt)
nstrongbind=$(grep -c -P "\tWB" $insample.sig_binders.txt)
echo -e "\tpredictions completed."
echo -e "\t\t$nstrongbind strong binders found"
echo -e "\t\t$nweakbind weak binders found"
echo -e "\tmapping strong and weak binders to lesion positions"

# match peptides to positions
# join -t$'\t' -1 $colnumneo -2 $colnmhcp $insample.neoantigens.bed $insample.sig_binders.txt | cut -f2- > $insample.sig_peptide_coord_matched.txt

# match all peptides to positions
join -t$'\t' -1 $colnumneo -2 $colnmhcp $insample.neoantigens.bed <(sort -k$colnmhcp"b",$colnmhcp $insample.all_peptides.txt) | cut -f2- > $insample.all_peptide_coord_matched.txt

# match significant peptides to positions
grep -P "\t[WS]B" $insample.all_peptide_coord_matched.txt > $insample.sig_peptide_coord_matched.txt
echo -e "joined peptides to positions" 1>&2

# print out numbers of coordinates
ncoordstrong=$(grep -P "\tSB" $insample.sig_peptide_coord_matched.txt | cut -f1-3 | sort -k1,1 | uniq)
ncoordweak=$(grep -P "\tWB" $insample.sig_peptide_coord_matched.txt | cut -f1-3 | sort -k1,1 | uniq)
echo -e "\t\t$ncoordstrong lesions mapped to strong binders"
echo -e "\t\t$ncoordweak lesions mapped to strong binders"

echo -e "\t$insample DONE" 1>&2


#
## 
#
#########
#    #
#    # run netMHCpan if entries exist
#    nlines=$(grep -c "^" $i.neoantigens.txt)
#    if [[ $nlines -gt 0 ]]; then
#        echo -e "$nlines non-truncated neoantigens found for $i. predicting affinities to $hlainput ..."
#        cut -f5 $i.neoantigens.txt > $i.neoantigens.netmhcpan_input.temp
#        netMHCpan \
#        -a $hlainput \
#        -f $i.neoantigens.netmhcpan_input.temp \
#        -s 1 \
#        -inptype 1 \
#        -BA 1 > $i.netmhcpan_out.txt
#        echo -e "\tpredictions completed. Summarizing strong and weak binders..."
#        #
#        # summarize outputs
#        sb=$(grep -s "<= SB" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4 | sort -k1,1 | uniq | grep -c "^")
#        wb=$(grep -s "<= WB" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4 | sort -k1,1 | uniq | grep -c "^")
#        
#        # sb=$(grep -s "<= SB" $i.netmhcpan_out.txt | uniq | grep -c "^")
#        # wb=$(grep -s "<= WB" $i.netmhcpan_out.txt | uniq | grep -c "^")
#        
#        
#        
#        echo -e "\t$sb strong binders and $wb weak binders found. mapping back to lesion position"
#        # grep -P "<= [SW]B" $i.netmhcpan_out.txt > $i.netmhcpan_out.strong_weak_binders.txt
#        # print strong binders
#        #echo -e "#" > $i.strong_binders.txt
#        #for j in $(grep -P "<= SB" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4); do
#        #    paste <(grep -P "\t$j" $i.neoantigens.txt) <(grep -P "$j" $i.netmhcpan_out.txt | grep -s "<= SB") >> $i.strong_binders.txt
#        #done
#        
#        #
#        echo -e "#" > $i.all_sig_binders.txt
#        while read peptidehit; do
#            # find the coordinate
#            # echo -e $peptidehit | tr ' ' '\t' | cut -f3
#            echo -e $peptidehit | tr ' ' '\t' > .peptidehit.mhc.temp
#            peptidehitseq=$(cut -f3 .peptidehit.mhc.temp)
#            # echo -e "searching for $peptidehitseq"
#            grep -P "$peptidehitseq" $i.neoantigens.txt > .peptidehit.coord.temp
#            while read coord; do
#                echo -e $coord | tr ' ' '\t' > .coord.peptide.temp
#                paste .peptidehit.mhc.temp .coord.peptide.temp >> $i.all_sig_binders.txt
#            done < .peptidehit.coord.temp
#        done < <(grep -s "PEPLIST" $i.netmhcpan_out.txt | grep -s "<= " | sed "s/<= //g" | sed 's/  */\t/g')
#        
#        # <(grep -s "PEPLIST" $i.netmhcpan_out.txt | grep -s "<= " | sed "s/<= //g" | sed 's/  */\t/g')
#        # <(grep -s "PEPLIST" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4)
#        
#            
#            #grep -P "\t$peptidehit\t" $neoantigengene | > .peptidehit.coord.temp
#            #grep -P "$peptidehit" $i.netmhcpan_out.txt > .peptidehit.mhcp.temp
#            #while read mhcp; do
#            #    echo -e $mhcp | tr ' ' '\t' > .mhcp.line
#            #    paste .peptidehit.coord.temp .mhcp.line >> $i.mhcp_out.snvcoord.txt
#            #done < .peptidehit.mhcp.temp
#        # done < <(grep -s "PEPLIST" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4)
#        # grep -s "<=\tSB" $i.all_binders.txt > $i.mhcp_out.snvcoord.strong_binders.txt
#        # grep -s "<=\tWB" $i.all_binders.txt > $i.mhcp_out.snvcoord.strong_binders.txt
#        
#        # print weak binders
#        #echo -e "#" > $i.weak_binders.txt
#        #for j in $(grep -P "<= WB" $i.netmhcpan_out.txt | sed 's/  */\t/g' | cut -f4); do
#        #    paste <(grep -P "\t$j" $i.neoantigens.txt) <(grep -P "$j" $i.netmhcpan_out.txt | grep -s "<= WB") >> $i.weak_binders.txt
#        #    # grep -P "\t$i" $i.neoantigens.txt
#        #done
#        #
#        
#        echo -e "\tstoring output..."
#        #
#        outfolder=$i"_peptide_mhc_results"
#        mkdir -p $outfolder
#        mv $i.hla_types.txt $i.neoantigens.txt $i.netmhcpan_out.txt $i.all_sig_binders.txt $outfolder
#        #
#        echo -e "sample: $i\tneoantigens: $nlines\tstrong_binders: $sb\tweak_binders: $wb"
#        echo -e "$i\t$nlines\t$sb\t$wb" >> batch1plus_neoantigen_counts.txt
#        #
#        rm $i.neoantigens.netmhcpan_input.temp
#        echo "\tdone. next"
#    else
#        echo -e "$nlines non-truncated neoantigens found for $i. skipping ..."
#        echo -e "$i\t0\t$sb\t$wb" >> batch1plus_neoantigen_counts.txt
#    fi
#
