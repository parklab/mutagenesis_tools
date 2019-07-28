#!/bin/bash

# input:
# a file containing the directory of RNA-seq BAMs and the breakpoint
# a configuration file with necessary inputs

infiles=$1
config=$2
nsamplequeue=$(grep -v "#" $config | grep -s "nsamplequeue" | cut -f2)
splitlines=$(grep -v "#" $config | grep -s "splitlines" | cut -f2)

## echo -e "intable\t$intable\nconfig\t$config"> input_params.txt
#
## get the number of files that we are running
#nfiles=$(grep -c "^" $infiles)
#
# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out
#
#echo "hello"$nfiles
#echo "hello"$nsamplequeue
#
## each job in the array will split into many other smaller jobs
#sbatch --array=1-$nfiles%$nsamplequeue /n/data1/hms/dbmi/park/vinay/pipelines/blast_tools/core_tblastn_crux_comet_peptides_sepoutdir.sh \
#               $infiles \
#               $config
#                   # $i $config $scriptdir $whichRNA


nfiles=$(grep -c "^" $infiles)
echo $nfiles
nbatches=$(expr $nfiles / $splitlines)

mkdir -p intables
intablename=$(basename $infiles)
intablename=${intablename%.*}
echo $intablename
split --lines=$splitlines --additional-suffix=".txt" $infiles intables/$intablename"_"

for i in intables/$intablename"_"*; do
    subbatchnfiles=$(grep -c "^" $i)
    echo $i
    # nfiles=$(grep -c "^" $i)
    echo -e "Submitting $nfiles files from $i"
    scriptdir=$(dirname "${BASH_SOURCE[0]}")
    sbatch --array=1-$subbatchnfiles%$nsamplequeue /n/data1/hms/dbmi/park/vinay/pipelines/blast_tools/core_tblastn_crux_comet_peptides_sepoutdir.sh \
                   $i \
                   $config
done
