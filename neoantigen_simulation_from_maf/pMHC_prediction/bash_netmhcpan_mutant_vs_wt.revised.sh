#!/bin/bash
# this script can serve as a bash template for batch-core operations on O2.

# load up necessary modules. Add more as necessary, but the basic ones are here
module load gcc/6.2.0

# the basic inputs are (1) an input file containing the tab-separted columns of input and (2) the configuration file
infile=$1 # each row of the input file contains the tab-separted list of input files that we need. 
config=$2

splitlines=$(grep -v "#" $config | grep -P "^splitlines\t" | cut -f2)
corescript=$(grep -v "#" $config | grep -P "corescript\t" | cut -f2)

# for very large files, sometimes you may need to split the job up into smaller batches that can then each be submitted as job arrays
nfiles=$(grep -c "^" $infile) # number of jobs to be run is the number of 
echo $nfiles
nbatches=$(expr $nfiles / $splitlines)

mkdir -p intables # to store the batches of files we are going to run
intablename=$(basename $infile)
intablename=${intablename%.*}
echo $intablename
split --lines=$splitlines --additional-suffix=".txt" $infile intables/$intablename"_" # 

# submit the job arrays one batch at a time
mkdir -p runlogs
mkdir -p runlogs/err
mkdir -p runlogs/out

dt=$(date '+%d_%m_%Y_%H_%M_%S')
echo -e "job\trun" > ./runlogs/runlogs.$dt.txt

jobidall=""
for i in intables/$intablename"_"*; do
	nfiles=$(grep -c "^" $i)
	echo -e "Submitting $nfiles files from $i ... Current sequence of job ids: $jobidall"
	scriptdir=$(dirname "${BASH_SOURCE[0]}") # the script directory, useful for downstream applications
	#sbatch --array=1-$nfiles%$nfiles $corescript $i $scriptdir # $config $scriptdir 
 	sbatch --array=1-$nfiles%$nfiles $corescript $i $config $scriptdir # THIS MODIFICATION REQUIRES THAT PER-JOB CONFIG IS SPECIFIED IN INFILES. IF NECESSARY, THE SAME CONFIG AS ABOVE CAN BE PLACED IN
	jobid=${currentjob#*job }
	echo -e "job\t$currentjob " >> ./runlogs/runlogs.$dt.txt
	jobidall=$jobidall":"$jobid
done
