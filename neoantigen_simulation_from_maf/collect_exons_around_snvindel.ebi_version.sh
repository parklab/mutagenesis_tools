#!/bin/bash

# load the conda environment with alignment tools # commented b/c we activate it above...
# source /homes/vvv776/anaconda3/etc/profile.d/conda.sh
# conda activate /homes/vvv776/anaconda3/envs/alignment



# OCTOBER 13, 2020 -- THE ERROR PREVENTING SOME MUTATIONS FROM BEING TRANSLATED IS OCCURRING HERE...

# inputs
infile=$1
config=$2
jobtask=$3
scriptdir=$4

genegtf=$(grep -v "#" $config | grep -s "genegtf" | cut -f2)
exongtf=$(grep -v "#" $config | grep -s "exongtf" | cut -f2)

mkdir -p temp_files # temp files -- this may have already been generated
mkdir -p peptide_sequences # temp files -- this may have already been generated


# 1. intersect the breakpoint with the gene GTFs and isolate the list of genes that we need to consider
bedtools intersect -wb -a $infile -b $genegtf > bkpt.$jobtask.gene.gtf
# frame as unique transcripts?
bedtools intersect -wb -a bkpt.$jobtask.gene.gtf -b $exongtf > bkpt.$jobtask.exon.gtf
bedtools intersect -wb -a <(cut -f1-3 bkpt.$jobtask.gene.gtf) -b $exongtf > bkpt.$jobtask.exon.gtf
# 2. get unique genes
cut -f20 bkpt.$jobtask.gene.gtf | sed "s/\"; /\t/g" | cut -f1 | sed "s/gene_id \"//g" | sort -k1,1 | uniq > $jobtask.unique_genes.txt
cut -f12 bkpt.$jobtask.exon.gtf | sed "s/\"; /\t/g" | cut -f2 | sed "s/transcript_id \"//g" | sort -k1,1 | uniq > $jobtask.unique_transcripts.txt

# 3. For each gene,
# a. get the lesions that intersect and generate the regions of the genes that lie outside of the lesions
# b. obtain the exons that lie within those extra-lesion regions and their FASTA sequences
# c. obtain the lesion sequence by applying the rule suggested by the type of lesion
# d. concatenate the extra-lesion and lesion sequences appropriately
echo -e "jobtask\tgene\tfasta\tseq_bed\tlesion_bed\tsurrounding_lesion_bed\tsurrounding_lesion_fasta\tproposed_neoa\tcontrol_neoa" > jobtask_gene_fasta.txt

# while read gene; do
for gene in $(cut -f1 $jobtask.unique_genes.txt); do
    # need to fix strand information
    # isolate the gene. NOTE that multiple sites may be lesioned
    echo -e "Now generating sequences for gene $gene"
    # grep -P "\"$gene\"" bkpt.$jobtask.gene.gtf > .temp.gene.gtf
    grep -P "\"$gene\"" bkpt.$jobtask.gene.gtf | sort -k1,1 | uniq > .temp.gene.gtf
    
    # generate the intervals outside the lesion
    echo -e "Now generating intervals that bracket lesion"
    cat <(cut -f12,15,18 .temp.gene.gtf | uniq) <(cut -f1,3,18 .temp.gene.gtf) > .temp.outside_lesion.starts.bed
    cat <(cut -f1,2,18 .temp.gene.gtf) <(cut -f12,16,18 .temp.gene.gtf | uniq) > .temp.outside_lesion.ends.bed
    paste <(cut -f1-2 .temp.outside_lesion.starts.bed) <(cut -f2-3 .temp.outside_lesion.ends.bed) | \
    awk -v gene=$gene 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,gene,$4,"WT","-","-","-"}' > .temp.outside_lesion.bed
    
    # cat .temp.outside_lesion.bed
    
    # add a line here to swap the start-and-end sites
    # awk 'BEGIN {FS=OFS="\t"} {if ($2 > $3) print $1,$3,$2,$0; else print $1,$2,$3,$0}' .temp.outside_lesion.old.bed | cut -f1-3,7- > .temp.outside_lesion.bed
    
    # get exons from the "outside_lesion" portion
    echo -e "Now generating exonic sequences that bracket the lesion"
    
    bedtools intersect -a $exongtf -b .temp.outside_lesion.bed > test_exo_lesion_temp.$gene.bed
    awk -v gene=$gene 'BEGIN {FS=OFS="\t"} {print $1,$4,$5,gene,$7,"WT","-","-","-"}' test_exo_lesion_temp.$gene.bed > test_exo_lesion_temp.$gene.2.bed
    sort -k1,1 -k2,2n test_exo_lesion_temp.$gene.2.bed > test_exo_lesion.$gene.bed
    
    cat test_exo_lesion.$gene.bed
    
    bedtools intersect -a $exongtf -b .temp.outside_lesion.bed | \
    awk -v gene=$gene 'BEGIN {FS=OFS="\t"} {print $1,$4,$5,gene,$7,"WT","-","-","-"}' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - -c 4,5,6,7,8,9 -o distinct | \
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3"_"$4,$5,$6,$7,$8,$9}' > .temp.outside_lesion.exons.bed
    
    # now, get the sequence within the lesion
    echo -e "Now getting sequence within lesion"
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$18,$7,$8,$9,$10}' .temp.gene.gtf > .temp.inside_lesion.bed
    cat .temp.outside_lesion.exons.bed .temp.inside_lesion.bed | sort -k1,1 -k2,2n > .temp.lesion_gene.bed
    # send the lesioned gene to the generate_lesion_sequence_from_exons
    bash $scriptdir/generate_sequence_around_lesion_from_exons.ebi_version.sh .temp.lesion_gene.bed $config $jobtask $gene $scriptdir
done
# done < $jobtask.unique_genes.txt


