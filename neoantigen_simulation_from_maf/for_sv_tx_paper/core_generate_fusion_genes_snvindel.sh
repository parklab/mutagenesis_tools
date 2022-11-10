#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o ./runlogs/out/out.%j.generate_fusion_genes_modified.txt #stdout file
#SBATCH -e ./runlogs/err/err.%j.generate_fusion_genes_modified.txt #stderr file


# submit support request to be added to #SBATCH -A park_contrib


# package loads
module load \
bedtools/2.27.1 \
samtools/1.9

# inputs
inbkpt=$1 # contains the MAF file of SNVs and indels. 
config=$2
scriptdir=$3

# GTFs from the config file
genegtf=$(grep -v "#" $config | grep -s "genegtf" | cut -f2)
exongtf=$(grep -v "#" $config | grep -s "exongtf" | cut -f2)
# get the coordinate system
coordsys=$(grep -v "#" $config | grep -s "coord_sys" | cut -f2)
# do we need to shuffle genes to look for chimeras?
shuffle=$(grep -v "#" $config | grep -s "shuffle" | cut -f2)

# input breakpoint
taskid=$SLURM_ARRAY_TASK_ID
jobid=$SLURV_ARRAY_JOB_ID
jobtask=$SLURM_JOB_ID
filename=$(sed -n $taskid'p' $inbkpt)
# sed -n $taskid'p' $inbkpt > bkpt.$jobtask.temp
# generate a TEMP file to convert the MAF into a BEDPE
# we will be working with the 0-based system by default
echo -e "$jobtask\t$filename" >> jobtask_file.txt

# output directory?
mkdir -p $jobtask
cd $jobtask
mkdir -p temp_files

cut -f5-7,9-10,35-36 $filename | sed "s/\t/_/g" > $jobtask.bkpt.txt
nbkpts=$(grep -c "^" $jobtask.bkpt.txt)
echo -e "$nbkpts variants found"
filenamebase=$(basename $filename)

if [ $coordsys == "1" ]; then
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6-2,$6-1,$5,$7,$7+1}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
    # $8,$9,$10,$11,$12,$13,$14}' $filename
    # awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$7,$8,$9,$10,$11,$12,$13,$14}' $filename | tail -n +2 > bkpt.internal.$jobtask.temp
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$7}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.internal.$jobtask.temp
    # https://stackoverflow.com/questions/1602035/how-to-print-third-column-to-last-column for faster iteration
    # paste <(cut -f5-7 $filename) <(cut -f35 $filename) <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
elif [ $coordsys == "0" ]; then
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$6,$5,$7,$7+1}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
    #awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6-1,$6,$5,$7,$7+1,$8,$9,$10,$11,$12,$13,$14}' $filename | paste - $jobtask.bkpt.txt | tail -n +2 > bkpt.$jobtask.temp
#     awk -v bkpt=$bkptname 'BEGIN {FS=OFS="\t"} {print $5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' $filename | tail -n +2 > bkpt.internal.$jobtask.temp
    awk 'BEGIN {FS=OFS="\t"} {print $5,$6,$7}' $filename | paste - $jobtask.bkpt.txt <(cut -f8-14 $filename) | tail -n +2 > bkpt.internal.$jobtask.temp
    # paste <(cut -f5-7 $filename) <(cut -f35 $filename) <(cut -f8-14 $filename) | tail -n +2 > bkpt.$jobtask.temp
else
    echo -e "INVALID COORDINATE SYSTEM SPECIFIED. Please specify \"0\" for zero-based coordinates\n
    (eg: X   0   1 in zero-based coordinates is the first position of chromosome X)\n
    or  \"1\" for one-based coordinates\n
    (eg: X   1   1 in one-based coordinates the first position of chromosome X)\n"
fi

echo -e "identified locations of variants"

# get the name of the breakpoint
 #bkptname=$(cut -f34 bkpt.$jobtask.temp)
 
# bkptname=$(cut -f34 $filename)

## start of comment

####==================================================================

###PART 1: GENERATE FUSION TRANSCRIPT GTF

## convert PART 1 into a script! 

# remove any existing version...
[ -e bkpt.$jobtask.fusion_transcripts.gtf ] && rm bkpt.$jobtask.fusion_transcripts.gtf

# run collect_exons_around_break.sh to execute items 1-5
if [ $shuffle == "FALSE" ]; then
    bash $scriptdir/collect_exons_around_snvindel.sh \
    ./bkpt.internal.$jobtask.temp \
    $config \
    $jobtask \
    $scriptdir
else
    bash $scriptdir/collect_exons_around_break_shuffle.sh \
    ./bkpt.$jobtask.temp \
    ./bkpt.internal.$jobtask.temp \
    $config \
    $jobtask \
    $scriptdir    
fi

# collect the neoantigens
for i in $(find ./ -name "*unique_neoas_mutant.txt"); do
    neoafile=$i
    for j in $(cut -f1 $neoafile); do
        echo -e "$neoafile\t$j\t$filenamebase" | sed "s/_/\t/g"
    done
done > $jobtask.mutant_neoantigens_per_gene.txt

# collect the WT peptides
for i in $(find ./ -name "*lost_neoas_mutant.txt"); do
    neoafile=$i
    for j in $(cut -f1 $neoafile); do
        echo -e "$neoafile\t$j\t$filenamebase" | sed "s/_/\t/g"
    done
done > $jobtask.unique_wt_neoantigens_per_gene.txt


## if the shuffle code works, then we will just go line-by-line, simulate FASTA sequences, and obtain the lesion. 
## then, we can forget about simulating transcripts. We can save that for later
## now, simulate protein sequences and pull reading frames
#for i in $(cut -f6 jobtask_gene_fasta.txt); do    
#    python /n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/propose_nres_peptides.py $i
#done


#
## items 8-10 will be their own scripts so that the lesion can be appropriately simulated as a sequence
## for mutations occurring in noncoding regions, provide some desired genomic length around the break.  
## 8. create every possible combination of transcripts from b1 and b2 to make the full fusion transcript,
## generate what the sequence between b1 and b2 should look like based on the type of lesion,
## and output the results
##
## the truly algorithmic way to reconstruct the intervening sequence would be to go back to the BAM file and infer the kind of sequence that mapped to the site.
## right now, we'll just specify the mutation type. If the input file is a MAF, then we'll extract the alternate allele that we need. if we have multiple, then we'll use both
## for instructions: we can produce the GTF for the intervening sequence between the breakpoints. However, our main output will be the FASTA sequence.
## We'll generate FASTAs for b1 and b2, but we'll generate the FASTA for the intervening sequence, decide how we should handle it (reverse it? produce another copy of it? ditch it?)
## then, the final output will be the FASTA, with different entries corresponding to the different parts of the sequence
## we'll produce sequences for both senses.
#
## note that producing the genes within the lesion is essentially the "inverse" of producing the genes outside the lesion! 
#
## maybe we can forget about assembling fusion transcripts and just assemble the exons. Transcript assembly takes way too long and could be done more efficiently and algorithmically based on RNA counts
## just collect the genes and assemble the different fusion genes.
#
#
#
#
#ntx1=$(grep -c "^" transcripts.b1.$jobtask.txt)
#ntx2=$(grep -c "^" transcripts.b2.$jobtask.txt)
#echo -e "$ntx1 transcripts terminate at end 1"
#echo -e "$ntx2 transcripts terminate at end 2"
#
## if both breakpoints are at transcripts
#echo -e "assembling fusion transcripts..."
#if [[ $ntx1 > 0 && $ntx2 > 0 ]]; then
#    while read tx1; do
#        grep -P "$tx1\"" bkpt.$jobtask.exon.b1_fusion.gtf > bkpt.$jobtask.b1_fusion_transcript.gtf.temp
#        while read tx2; do
#            grep -P "$tx2\"" bkpt.$jobtask.exon.b2_fusion.gtf > bkpt.$jobtask.b2_fusion_transcript.gtf.temp
#            fusionname=$tx1"_"$tx2"_"$bkptname
#            echo -e $fusionname >> fusion_tx.$jobtask.txt
#            cat bkpt.$jobtask.b1_fusion_transcript.gtf.temp bkpt.$jobtask.b2_fusion_transcript.gtf.temp | \
#            awk -v fxname=$fusionname 'BEGIN {FS=OFS="\t"} {print $0" fusion_transcript_id \""fxname"\""}' >> bkpt.$jobtask.fusion_transcripts.gtf
#        done < transcripts.b2.$jobtask.txt
#    done < transcripts.b1.$jobtask.txt
#fi
#
## if one breakpoint is at a transcript
#if [[ $ntx1 == 0 && $ntx2 > 0 ]]; then
#    while read tx2; do
#        grep -P "$tx2\"" bkpt.$jobtask.exon.b2_fusion.gtf > bkpt.$jobtask.b2_fusion_transcript.gtf.temp
#        fusionname="INTERGENIC_"$tx2"_"$bkptname
#        echo -e $fusionname >> fusion_tx.$jobtask.txt
#        cat bkpt.$jobtask.b2_fusion_transcript.gtf.temp | \
#        awk -v fxname=$fusionname 'BEGIN {FS=OFS="\t"} {print $0" fusion_transcript_id \""fxname"\""}' >> bkpt.$jobtask.fusion_transcripts.gtf
#    done < transcripts.b2.$jobtask.txt
#fi
#
#if [[ $ntx2 == 0 && $ntx1 > 0 ]]; then
#    while read tx1; do
#        grep -P "$tx1\"" bkpt.$jobtask.exon.b1_fusion.gtf > bkpt.$jobtask.b1_fusion_transcript.gtf.temp
#        fusionname=$tx1"_INTERGENIC_"$bkptname
#        echo -e $fusionname >> fusion_tx.$jobtask.txt
#        cat bkpt.$jobtask.b1_fusion_transcript.gtf.temp | \
#        awk -v fxname=$fusionname 'BEGIN {FS=OFS="\t"} {print $0" fusion_transcript_id \""fxname"\""}' >> bkpt.$jobtask.fusion_transcripts.gtf
#    done < transcripts.b1.$jobtask.txt
#fi
#
#if [[ $ntx2 == 0 && $ntx1 == 0 ]]; then
#    echo -e "" >> bkpt.$jobtask.fusion_transcripts.gtf
#    echo -e "no fusion transcripts found; intergenic transcript suspected"
#    # here, you should instead propose a default 1-MB window 
#    exit 1
#fi
#
## 9. Based on the expanse of the gene which was altered, return a GTF with the sequence to be altered and provide the proposed sequence
## for the action, add in the config file a key about what to do with certain mutations.
## three possible actions: remove, copy, or swap
#mutype=$(cut -f10 ./bkpt.$jobtask.temp)
#mutaction=$(grep -v "#" $config | grep -s $mutype | cut -f2)
## bash generate_lesion_sequence.sh $lesionpos $genomeref $mutaction $config # need to write this code up.
## for now, since we care about SNV/indels,
#bash $scriptdir/generate_lesion_sequence.sh bkpt.internal.$jobtask.temp $mutaction $config $jobtask
#
##### next start of comment
##
### these transcripts should then be piped into the NONCODING TRANSCRIPTION IDENTIFICATION PIPELINE
### add a line to incorporate the change
##
### 10. Simulate the new transcript sequence
### get the FASTA sequences. Make sure you adjust based on strand
##echo -e "getting FASTA sequences of exons of fusion transcripts..."
##for i in $(cut -f1 fusion_tx.$jobtask.txt | sort -k1,1 | uniq); do
##    # bedtools getfasta -s -fi $fasta -bed <(grep -s "$i" bkpt.$jobtask.fusion_transcripts.gtf) > bkpt.$i.fusion_transcripts.fa #GET SEPARTE SEQUENCES FOR EACH END
##    bedtools getfasta -s -fi $fasta -bed <(grep -s "$i" bkpt.$jobtask.b1_fusion_transcript.gtf.temp) > bkpt.$i.b1_fusion_transcripts.fa #GET SEPARTE SEQUENCES FOR EACH END
##    bedtools getfasta -s -fi $fasta -bed <(grep -s "$i" bkpt.$jobtask.b2_fusion_transcript.gtf.temp) > bkpt.$i.b2_fusion_transcripts.fa #GET SEPARTE SEQUENCES FOR EACH END
##    # AT THIS STEP, APPEND THE LESION SEQUENCE SUCH THAT 1-LESION-2 IS PRESENT # adjust the > header so that the sequence is one unit? # no need.
##    cat bkpt.$i.b1_fusion_transcripts.fa $jobtask"_allele.fa" bkpt.$i.b2_fusion_transcripts.fa > bkpt.$i.fusion_transcripts.fa
##    python /n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/propose_nres_peptides.py bkpt.$i.fusion_transcripts.fa # PROBABLY ASK THE 
##    mv bkpt.$i.fusion_transcripts.fa $jobtask"_fastas" # *$jobtask.fusion_transcripts.fa fasta_files
##    # for this analysis, append the exon FASTA sequences together and do the tiling analysis. Then, propose peptides
##    mv bkpt.$i.*fusion_transcripts*_peptides_*aa.txt $jobtask"_peptides"
##    # submit the exons to the /Users/vinayakvsv/Desktop/mnt.transfer/pipelines/mutagenesis/propose_nres_peptides.py script
##    # this python script contains packages; a shell script to use these packages should be built to accept FASTA input
##    # store the sequences, the peptide residues, and the subset involving the fusion
##done
##
##
##
### some of the fields fail to generate "closest" exons. The intergenics?" NO. something else...
##
### add to the global
##echo -e "adding to global..."
##cat bkpt.$jobtask.fusion_transcripts.gtf >> bkpt.$jobid.fusion_transcripts.gtf
##
#### end of comment
#
