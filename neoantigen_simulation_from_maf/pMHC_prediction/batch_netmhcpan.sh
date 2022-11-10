#!/bin/bash

samplefile=$1
config=$2

nsamples=$(grep -c "^" $samplefile)
nsamplequeue=$(grep -s "nsamplequeue" $config | cut -f2)
nbatches=$(expr $nsamples / $nsamplequeue)

echo -e "samplefile\t$samplefile\nconfig\t$config" > input_params.txt

# create a runlogs directory
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

scriptdir=$(dirname "${BASH_SOURCE[0]}")

# each job in the array will split into many other smaller jobs
echo -e "#jobtask\tsample\tfhla" > jobtask_insample_inhla.txt

mkdir -p intables_batch_netmhcpan
intablename=$(basename $samplefile)
intablename=${intablename%.*}
echo $intablename
split --lines=$nsamplequeue --additional-suffix=".txt" $samplefile ./intables_batch_netmhcpan/$intablename"_"

for i in intables_batch_netmhcpan/$intablename"_"*; do
    subbatchnfiles=$(grep -c "^" $i)
    echo -e "Submitting $subbatchnfiles samples from $i"
    sbatch --array=1-$subbatchnfiles%$subbatchnfiles $scriptdir/core_netmhcpan.sh \
               $i \
               $config \
               $scriptdir
done
