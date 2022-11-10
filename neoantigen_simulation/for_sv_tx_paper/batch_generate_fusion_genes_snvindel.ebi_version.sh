#!/bin/bash

# inputs
breakpoints=$1
config=$2

# how many individual MAFs do we have to work with
nbkpts=$(grep -c "^" $breakpoints)
nsamplequeue=$(grep -s "nsamplequeue" $config | cut -f2)
corescript=$(grep -P "^corescript\t" $config | cut -f2)

# write up the inputs to a separate file
echo -e "breakpoints\t$breakpoints\n
config\t$config" > input_params.txt

# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

# what is the directory of this script...
scriptdir=$(dirname "${BASH_SOURCE[0]}")

# associate the job id with the file used
echo -e "jobtask\tfile" >> jobtask_file.txt

# how many batches
nbatches=$(expr $nbkpts / $nsamplequeue)

# create of jobs to submit (so that we can submit queues)
mkdir -p intables_batch_get_coverage_at_intervals
intablename=$(basename $breakpoints)
intablename=${intablename%.*}
echo $intablename
split --lines=$nsamplequeue --additional-suffix=".txt" $breakpoints intables_batch_get_coverage_at_intervals/$intablename"_"

for isubfile in intables_batch_get_coverage_at_intervals/$intablename"_"*; do
    subbatchnfiles=$(grep -c "^" $isubfile)
    subname=$(basename $isubfile)
    subname=${subname%.txt}
    echo $isubfile
    echo -e "Submitting intervals from $nfiles files from $isubfile"
    # for j in $(cut -f1 $i); do
    #     cat $j
    # done > $i.temp.bed
    scriptdir=$(dirname "${BASH_SOURCE[0]}")
    
    ### CONVERT THIS TO AN LSF-COMPATIBLE SUBMISSION
    bsub -J generate_snv_indel_neoa[1-$subbatchnfiles] -o ./runlogs/out/out.$subname.txt -e ./runlogs/err/err.$subname.txt -M 48G \
    bash $corescript $isubfile $config $scriptdir >> runlogs/runlogs.$subname.txt
    echo -e 'waiting for 10 seconds...'
    sleep 10
    
    #sbatch --array=1-$subbatchnfiles%$subbatchnfiles $scriptdir/core_generate_fusion_genes_snvindel.ebi_version.sh \
    #           $isubfile \
    #           $config \
    #           $scriptdir
    ## no wonder only a fraction of my jobs submitted! I inputted $breakpoints instead of $i! # this was fixed
    
done
