#!/bin/bash

# the samplefile will contain a list of samples that we would like to run
samplefile=$1
# the config file will contain
# 1. a path to the mutant peptides
# 2. a path to the WT peptides
# 3. a genome file
config=$2

# get the number of samples that we need
nsamples=$(grep -c "^" $samplefile)
# how many should we run in one job array
nsamplequeue=$(grep -s "nsamplequeue" $config | cut -f2)
# number of such job arrays that we will submit
nbatches=$(expr $nsamples / $nsamplequeue)

echo -e "samplefile\t$samplefile\nconfig\t$config" > input_params.txt

# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

scriptdir=$(dirname "${BASH_SOURCE[0]}")

# each job in the array will split into many other smaller jobs
echo -e "#jobtask\tsample\tfhla" > jobtask_insample_inhla.txt

mkdir -p intables_batch_mutant_wt_peptides
intablename=$(basename $samplefile)
intablename=${intablename%.*}
echo $intablename
split --lines=$nsamplequeue --additional-suffix=".txt" $samplefile ./intables_batch_mutant_wt_peptides/$intablename"_"


for i in intables_batch_mutant_wt_peptides/$intablename"_"*; do
    subbatchnfiles=$(grep -c "^" $i)
    echo -e "Submitting $subbatchnfiles samples from $i"
    sbatch --array=1-$subbatchnfiles%$subbatchnfiles $scriptdir/core_match_mutant_wt.updated.sh \
               $i \
               $config \
               $scriptdir
done
