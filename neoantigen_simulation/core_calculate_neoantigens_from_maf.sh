#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o ./runlogs/out/%j.core_calculate_neoantigens.out.txt #stdout file
#SBATCH -e ./runlogs/err/%j.core_calculate_neoantigens.err.txt #stderr file

## Step 0: inputs
infile=$1 # this should contain the list of BARISTA outputs
config=$2
scriptdir=$3

# source
source activate cell_line_analysis

# get the SRR accession number to be downloaded
taskid=$SLURM_ARRAY_TASK_ID
outname=$(sed -n $taskid'p' $infile | cut -f1)
intervals=$(sed -n $taskid'p' $infile | cut -f2)

# extract from Config:
exongtf=$(grep -v "#" $config | grep -P "^exongtf\t" | cut -f2)
genomefa=$(grep -v "#" $config | grep -P "^genomefa\t" | cut -f2)
sloplen=$(grep -v "#" $config | grep -P "^sloplen\t" | cut -f2)

mkdir -p $outname
cd $outname

echo -e "sample\t$outname\n\
intervals\t$intervals\n\
exongtf\t$exongtf\n\
genomefa\t$genomefa\n\
sloplen\t$sloplen\n\
####" > run_params.txt

#
# python $scriptdir/get_neoantigens_at_mutations.py
if [[ $exongtf == '' ]]; then
    echo -e "getting genome context" >> run_params.txt
    python $scriptdir/get_neoantigens_at_mutations.py --input $intervals --fasta $genomefa --name $outname --distance $sloplen
else
    echo -e "getting annotation-sequence context" >> run_params.txt
    python $scriptdir/get_neoantigens_at_mutations.py --input $intervals --fasta $genomefa --name $outname --distance $sloplen --annotation $exongtf
fi

# move tempfiles
echo -e "moving tempfiles away..."
mkdir -p tempfiles
mv .temp.*.bed tempfiles

# to debug
for i in $(awk 'BEGIN {FS=OFS="\t"} {if (length($8) < 5) print $2}' $outname.mutation_sequence_contexts.txt); do
    grep -P "\t$i\t" $intervals
done > $outname.to_debug.bed
grep -c "^" $outname.to_debug.bed | paste - <(echo -e "variants don't have an annotation sequence context generated. Please recheck") >> run_params.txt