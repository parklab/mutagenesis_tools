#!/usr/bin/env python

import os,sys,re
import importlib
import argparse
import time

# start timing
start_time = time.time()

# import the main function
import mutant_neoantigen_simulation as neoa

# set up argparse
# place an argparse interface here
parser = argparse.ArgumentParser(description='Simulate neoantigens around a sequence')
parser.add_argument('-i','--input',
                    help='BED file input containing variants. We suggest the "chr,start,end,ref,alt,strand" format. NOTE: Your coordinates should be specified as one-based (i.e. a position is marked by a number, not by the span between two numbers). For insertions, the coordinate should be marked with two positions, the start right before the insertion, and the end right after the insertion. REQUIRED',type=str)
parser.add_argument('-a','--annotation',
                    help='Annotation file. Optional, but suggested if you want to conduct meaningful sequence translations. REQUIRED',type=str,default=None)
parser.add_argument('-f','--fasta',
                    help='Genome fasta. REQUIRED',type=str)
parser.add_argument('-n','--name',
                    help='Job name. DEFAULT: job',type=str,default="job")
parser.add_argument('-d','--distance',
                    help='Distance on either side of the mutation to extract sequence context. DEFAULT: 36 bases',type=int,default=36)
parser.add_argument('-j','--jobs',
                    help='Number of jobs to run in parallel',type=int,default=10)


# add an option to handle a zero-based coordinates


args = parser.parse_args()

print("Arguments:")
print(args)

intervals = args.input
genomefa = args.fasta
annot = args.annotation
jobname = args.name
dist = args.distance

print(intervals,genomefa,annot,jobname,dist)

# first, create the set of mutants
mutant_file = neoa.mutation_set(intervals,
                                ref_genome = genomefa,
                                gtf = annot)

# then, for each mutant, derive the annotation-defined sequence contexts
# if annotation is specified, then get the annotation context
if annot is not None:
    mutant_file.get_seq_context_at_annot(equalslop=dist,clip_at_stop_codons=True,store_different_only=True,remove_intersect_file_at_end=False)
    # translate the sequence contexts
    # mutant_file.get_translations_at_context()
else:
    mutant_file.get_seq_context(equalslop=dist)



# finally, write out two separate tables -- mutants with sequence contexts, and mutants with pairs of peptides
# sequence contexts
with open(jobname+".mutation_sequence_contexts.txt",'w') as f, open(jobname+".mutant_vs_wt_peptide.txt",'w') as a:
    outheader = ['#chr','start','stop','ref','alt','strand','wt_seq_context','mut_seq_context','meta']
    outheader2 = ['#chr','start','stop','ref','alt','mut_peptide','wt_peptide','meta']
    f.write('\t'.join(outheader)+"\n")
    a.write('\t'.join(outheader2)+"\n")
    for i in mutant_file:
        f.write(i.__str__()+"\t"+','.join(i.annot_strand)+"\n")
        # f.write(i.__str__()+"\n")
        if len(i.mutant_vs_wt_pairs) > 0:
            # only activate this option if we have reported mutant-vs-wt pairs
            for j in i.mutant_vs_wt_pairs:
                print(j)
                a.write("\t".join(i.coordinate) + "\t" +  "\t".join([j[0],j[1],str(j[2])]) + "\t" + "\t".join([str(m) for m in i.metadata]) + "\n")

# end time
end_time = time.time()
elapsed_time = end_time - start_time
print("elapsed time in seconds for",jobname,elapsed_time)
