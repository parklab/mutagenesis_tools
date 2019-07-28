#!/bin/bash

# input:
# a file containing the directory of RNA-seq BAMs and the breakpoint
# a configuration file with necessary inputs

infiles=$1
config=$2
nsamplequeue=$(grep -v "#" $config | grep -s "nsamplequeue" | cut -f2)

# echo -e "intable\t$intable\nconfig\t$config"> input_params.txt

# get the number of files that we are running
nfiles=$(grep -c "^" $infiles)

# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

# each job in the array will split into many other smaller jobs
sbatch --array=1-$nfiles%$nsamplequeue /n/data1/hms/dbmi/park/vinay/pipelines/blast_tools/core_tblastn_crux_comet_peptides.sh \
               $infiles \
               $config
                   # $i $config $scriptdir $whichRNA
