#!/bin/bash
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 0-12:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=12G
#SBATCH -o %j.out.txt
#SBATCH -e %j.err.txt

paste <(cut -f1 all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.reheader.txt | tail -n +2 | sed "s/~/\t/g" | cat <(echo -e "chr\tstart\tend\tref\talt") - ) <(cut -f2- all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.reheader.txt) > all_files_matched_mutant_wt_pep_pair.simple_affinity_comparison.reheader.bed

