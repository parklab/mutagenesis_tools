#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-12:00
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o %j.core_makeblastdb.out.txt #stdout file
#SBATCH -e %j.core_makeblastdb.err.txt #stderr file

module load blast/2.6.0+

infasta=$1
name=$2
outdir=$3

makeblastdb -in $infasta -input_type fasta -dbtype nucl -title $name -parse_seqids -out $name
BLASTDB=$BLASTDB:$outdir

# makeblastdb -in /n/data1/hms/dbmi/park/vinay/referenceGenomes/hg19_ucsc.fa -input_type fasta -dbtype nucl -title hg19_ucsc_blastdb -parse_seqids -out hg19_ucsc_blastdb
# makeblastdb -in /n/data1/hms/dbmi/park/vinay/referenceGenomes/hg19_ucsc.fa -input_type fasta -dbtype nucl -title hg19_ucsc_blastdb -parse_seqids -out hg19_ucsc_blastdb

# makeblastdb -in ./hg19_ucsc.fa -dbtype nucl -parse_seqids -out hg19_ucsc_blastdb.db -logfile makeblastdb_hg19_ucsc_db.txt
