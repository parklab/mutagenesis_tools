#!/bin/bash

breakpoints=$1
config=$2

nbkpts=$(grep -c "^" $breakpoints)
nsamplequeue=$(grep -s "nsamplequeue" $config | cut -f2)
# splitlines=$(grep -v "#" $config | grep -s "splitlines" | cut -f2)

echo -e "breakpoints\t$breakpoints\n
config\t$config" > input_params.txt

# get the number of files that we are running
# nfiles=$(grep -c "^" $intable)

# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

scriptdir=$(dirname "${BASH_SOURCE[0]}")

# each job in the array will split into many other smaller jobs
#sbatch --array=1-$nbkpts%$nsamplequeue
echo -e "jobtask\tfile" >> jobtask_file.txt

#sbatch --array=1-$nbkpts%$nsamplequeue $scriptdir/core_generate_fusion_genes_snvindel.sh \
#               $breakpoints \
#               $config \
#               $scriptdir
#                   # $i $config $scriptdir $whichRNA


nbatches=$(expr $nbkpts / $nsamplequeue)

mkdir -p intables_batch_get_coverage_at_intervals
intablename=$(basename $breakpoints)
intablename=${intablename%.*}
echo $intablename
split --lines=$nsamplequeue --additional-suffix=".txt" $breakpoints intables_batch_get_coverage_at_intervals/$intablename"_"

for isubfile in intables_batch_get_coverage_at_intervals/$intablename"_"*; do
    subbatchnfiles=$(grep -c "^" $isubfile)
    echo $isubfile
    echo -e "Submitting intervals from $nfiles files from $isubfile"
    # for j in $(cut -f1 $i); do
    #     cat $j
    # done > $i.temp.bed
    scriptdir=$(dirname "${BASH_SOURCE[0]}")
    
    sbatch --array=1-$subbatchnfiles%$subbatchnfiles $scriptdir/core_generate_fusion_genes_snvindel.sh \
               $isubfile \
               $config \
               $scriptdir
    # no wonder only a fraction of my jobs submitted! I inputted $breakpoints instead of $i! # this was fixed
    
done
