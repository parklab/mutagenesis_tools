#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o %j.generate_fusion_genes_modified.txt #stdout file
#SBATCH -e %j.generate_fusion_genes_modified.txt #stderr file


# inputs
infile=$1
config=$2
jobtask=$3
# bkptname=$5
scriptdir=$4

echo -e $scriptdir


genegtf=$(grep -v "#" $config | grep -s "genegtf" | cut -f2)
exongtf=$(grep -v "#" $config | grep -s "exongtf" | cut -f2)


# replace genegtf and exongtf with config and interval sequence

mkdir -p temp_files
mkdir -p peptide_sequences


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

while read gene; do
    # need to fix strand information
    # isolate the gene. NOTE that multiple sites may be lesioned
    echo -e "Now generating sequences for gene $gene"
    grep -P "\"$gene\"" bkpt.$jobtask.gene.gtf > .temp.gene.gtf
    
    # generate the intervals outside the lesion
    cat <(cut -f12,15,18 .temp.gene.gtf | uniq) <(cut -f1,3,18 .temp.gene.gtf) > .temp.outside_lesion.starts.bed
    cat <(cut -f1,2,18 .temp.gene.gtf) <(cut -f12,16,18 .temp.gene.gtf | uniq) > .temp.outside_lesion.ends.bed
    paste <(cut -f1-2 .temp.outside_lesion.starts.bed) <(cut -f2-3 .temp.outside_lesion.ends.bed) | \
    awk -v gene=$gene 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,gene,$4,"WT","-","-","-"}' > .temp.outside_lesion.bed
    
    # get exons from the "outside_lesion" portion
    bedtools intersect -a $exongtf -b .temp.outside_lesion.bed | \
    awk -v gene=$gene 'BEGIN {FS=OFS="\t"} {print $1,$4,$5,gene,$7,"WT","-","-","-"} ' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - -c 4,5,6,7,8,9 -o distinct | \
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$1"_"$2"_"$3"_"$4,$5,$6,$7,$8,$9}' > .temp.outside_lesion.exons.bed
    
    # now, get the sequence within the lesion
    awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$18,$7,$8,$9,$10}' .temp.gene.gtf > .temp.inside_lesion.bed
    cat .temp.outside_lesion.exons.bed .temp.inside_lesion.bed | sort -k1,1 -k2,2n > .temp.lesion_gene.bed
    # send the lesioned gene to the generate_lesion_sequence_from_exons
    bash $scriptdir/generate_sequence_around_lesion_from_exons.sh .temp.lesion_gene.bed $config $jobtask $gene
done < $jobtask.unique_genes.txt




#
## since this is breakpoint 1, we need to get the left-most coordinate of the gtf in order to generate the fusion gene culminating in b1 
#awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$7,$8,$3,$10,$11,$12,$NF" breakpoint \""$4"\";"}' bkpt.$jobtask.gene.b1.gtf.temp > bkpt.$jobtask.gene.b1_fusion.gtf
#nfus1=$(grep -c "^" bkpt.$jobtask.gene.b1_fusion.gtf)
#echo -e "$nfus1 genes found for end 1 of the breakpoint"
#
## 2. find the stretch from the gene gtf within the affected genes from #1
#bedtools intersect -wb -a <(cut -f4-6,7 $infile) -b $genegtf > bkpt.$jobtask.gene.b2.gtf.temp
## since this is breakpoint 1, we need to get the left-most coordinate of the gtf in order to generate the fusion gene culminating in b1 
#awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$7,$2,$9,$10,$11,$12,$NF" breakpoint \""$4"\";"}' bkpt.$jobtask.gene.b2.gtf.temp > bkpt.$jobtask.gene.b2_fusion.gtf
#nfus2=$(grep -c "^" bkpt.$jobtask.gene.b2_fusion.gtf)
#echo -e "$nfus2 genes found for end 2 of the breakpoint"
#
## since we are not shuffling genes, we will just identify each unique gene in our dataset and extract all of its associated exons
#cat bkpt.$jobtask.gene.b1_fusion.gtf bkpt.$jobtask.gene.b2_fusion.gtf | cit =f9 | sed "s/; /\t/g" | cut -f1 | sort -k1,1 | uniq > $jobtask.unique_genes.txt
## go gene by gene and extract the exons. Then, collect the exons for each breakpoint end.
#echo -e "jobtask\tgene1\tb1\tgene2\tb2" > $jobtask.breakpoint_ends.txt
#while read gene; do
#    # get exons for end 1
#    bash $scriptdir/get_exons_within_gene_gtf.sh $exongtf bkpt.$jobtask.gene.b1_fusion.gtf $jobtask $gene 1
#    
#    #bedtools intersect -wb -a $exongtf -b <(grep -P "$gene" bkpt.$jobtask.gene.b1_fusion.gtf) | sort -k1,1 -k4,4n > bkpt.$jobtask.$gene.orig_exon.b1_fusion.gtf.temp
#    #paste <(cut -f1-9 bkpt.$jobtask.$gene.orig_exon.b1_fusion.gtf.temp) <(cut -f18 bkpt.$jobtask.$gene.orig_exon.b1_fusion.gtf.temp | sed "s/; /;\t/g" | rev | cut -f1 | rev) | sed "s/;\t/; /g" > bkpt.$jobtask.$gene.orig_exon.b1_fusion.gtf
#    #nfus1=$(grep -c "^" bkpt.$jobtask.$gene.orig_exon.b1_fusion.gtf)
#    #echo -e "$nfus1 exons found for end 1 of the breakpoint at gene $gene"
#    
#    # get exons for end 2
#    bash $scriptdir/get_exons_within_gene_gtf.sh $exongtf bkpt.$jobtask.gene.b2_fusion.gtf $jobtask $gene 2
#    
#    # extract the breakpoint names
#    
#    
#    # print out the names of the GTFs, the genes they correspond to, and the job task all on a table
#    echo -e "$jobtask\t$gene\t$PWD/bkpt.$jobtask.gene.b1_fusion.gtf\t$gene\t$PWD/bkpt.$jobtask.gene.b2_fusion.gtf" >> $jobtask.breakpoint_ends.txt
#    
#    ## grep -P "$gene" bkpt.$jobtask.gene.b1_fusion.gtf > bkpt.$jobtask.$gene.b1.gtf
#    #bedtools intersect -wb -a $exongtf -b <(grep -P "$gene" bkpt.$jobtask.gene.b2_fusion.gtf) | sort -k1,1 -k4,4n > bkpt.$jobtask.$gene.orig_exon.b2_fusion.gtf.temp
#    #paste <(cut -f1-9 bkpt.$jobtask.$gene.orig_exon.b2_fusion.gtf.temp) <(cut -f18 bkpt.$jobtask.$gene.orig_exon.b2_fusion.gtf.temp | sed "s/; /;\t/g" | rev | cut -f1 | rev) | sed "s/;\t/; /g" > bkpt.$jobtask.$gene.orig_exon.b2_fusion.gtf
#    #nfus1=$(grep -c "^" bkpt.$jobtask.$gene.orig_exon.b2_fusion.gtf)
#    #echo -e "$nfus1 exons found for end 1 of the breakpoint at gene $gene"
#    
#    ## get FASTA sequences for each end # need to define fasta
#    #bash $scriptdir/generate_lesion_sequence.sh bkpt.internal.$jobtask.temp $mutaction $config $jobtask
#    #
#    #bedtools getfasta -s -fi $fasta -bed <(grep -s "$i" bkpt.$jobtask.gene.b1_fusion.exons.gtf) > bkpt.$jobtask.gene.b1_fusion.exons.fa
#    #bedtools getfasta -s -fi $fasta -bed <(grep -s "$i" bkpt.$jobtask.gene.b2_fusion.exons.gtf) > bkpt.$jobtask.gene.b1_fusion.exons.fa
#    
#    
#done < $jobtask.unique_genes.txt
#
