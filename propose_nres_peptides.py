import os,sys,re
import numpy,math,random
import itertools
import numpy as np
from collections import defaultdict
import argparse

# sys.path.insert(")
# import make_all_nmer_substitutions as nmersub
# from make_all_nmer_substitutions import reverse_seq,translation,codon_table,tile_sequence,complement_seq
from make_all_nmer_substitutions_copy_dev import reverse_seq,translation,codon_table,tile_sequence,complement_seq

def import_sequence(infa,concat_seq=True):
    inseqname = []
    inseq = []
    # infile = open(infa)
    # for i in infile:
    with open(infa) as infafile:
        for i in infafile:
            if i[0] != ">":
                inseq.append(i.strip("\n"))
            else:
                inseqname.append(i.strip("\n")[1:])
    if concat_seq:
        seq_full = ''.join(inseq)
        seqname_full = ';'.join(inseqname)
    else:
        seq_full = inseq
        seqname_full = inseqname
    return(seq_full,seqname_full)

def reverse_complement(inseq):
    return(reverse_seq(complement_seq(inseq)))

def translate_sequence(inseq):
    return(translation(inseq))

def start_stop_codons():
    codon_table_in = codon_table()
    start_codon = [i for i in codon_table_in if codon_table_in[i] == "M"]
    stop_codon = [i for i in codon_table_in if codon_table_in[i] == "*"]
    return(start_codon,stop_codon)

def define_frames(inseq,search_start_codons=True,search_stop_codons=True,start_seq_site_lim=3):
    # we will save the start and stop sites
    start_seq_sites = []
    end_seq_sites = []
    frames = []
    start_codon,stop_codons = start_stop_codons()
    # print(start_codon,stop_codons)
    ## codons = nmersub.tile_sequence(inseq = inseq,k=3,overlap=False)
    ## search for start codon sequences
    ## re_startcodon = re.compile(start_codon)
    inseq = inseq.lower()
    if search_start_codons:
        start_codons_regex = '|'.join(['(?=('+str(i)+'))' for i in start_codon])
        start_seq_sites = [i.start() for i in re.finditer(start_codons_regex,inseq)]
        # if groups of start seq sites are within `start_seq_site_lim` bases of each other, pick the most upstream start site

    else:
        start_seq_sites = [0]
    # search for stop-codon sequences
    if search_stop_codons:
        # generate stop codons regex
        stop_codons_regex = '|'.join(['(?=('+str(i)+'))' for i in stop_codons])
        # print(stop_codons_regex)
        ## search for non-overlapping instances of start codons--but what if these start codons are gateways to new frames?
        stop_seq_sites = [i.start()+3 for i in re.finditer(stop_codons_regex,inseq)] # add a 3 to provide the index after the stop codon ends
    else:
        stop_seq_sites = [len(inseq)]

    # note that if we aren't looking for frames but would like to find an in-frame translation, adjust the start and stop sites accordingly so that we can get appropriate frames
    if not search_start_codons and not search_stop_codons:
        other_start = ((stop_seq_sites[0] - start_seq_sites[1]) % 3) - 1
        other_end = stop_seq_sites[0] + ((stop_seq_sites[0] - start_seq_sites[1]) % 3) - 1
        start_seq_sites.append(other_start)
        stop_seq_sites.append(other_end)
    # print(start_seq_sites,stop_seq_sites)


    # search every combination of start and stop codons; keep frames that are divisible by 3.
    frame_positions = [(i,j) for i in start_seq_sites for j in stop_seq_sites if (j - i) % 3 == 0 and j > i]
    # print(frame_positions)
    frames = [inseq[i[0]:i[1]] for i in frame_positions]

    return(frames)

def write_peptides(peptide_array,name):
    # outfile = open(os.getcwd()+"/"+name+".txt",'w')
    outfile = open(name+".txt",'w')
    for i in peptide_array:
        for j in i:
            outfile.write(j+"\n")
    outfile.close()

def length_fiter(peptide,minlen=8):
    return(len(peptide) >= minlen)

def write_peptide_frames_old(peptide_array,name,filename,lenfilter=10):
    # write as a FASTA file
    outfile = open(filename+".fa",'a')
    count = 0
    for i in peptide_array:
        # outfile.write(">"+name+"_"+str(count)+"\n")
        for j in i:
            if len(j) >= lenfilter:
                outfile.write(">"+name+"_"+str(count)+"\n")
                outfile.write(j+"\n")
            count += 1
    outfile.close()

def write_peptide_frames(peptide_array,name,lenfilter=10):
    # write as a FASTA file
    outfile = open(name+".fa",'w')
    count = 0
    for i in peptide_array:
        outfile.write(">"+name+"_"+str(count)+"\n")
        # if len(j) >= lenfilter:
        #     outfile.write(j+"\n")
        if len(i) >= lenfilter:
            outfile.write(i+"\n")
        count += 1
    outfile.close()


def define_peptides(inseq,name,peptide_lengths=[8,9,10,11],search_start_codons=True,search_stop_codons=True,search_alt_strand=False,digest_peptides=True):
    # if search_start_codons=True, then get all possible start sites--any position with an M. Else, start from the beginning of the sequence
    # if search_stop_codons=True, then get all possible end sites--any codon with a stop codon. Else, end at the beginning of the sequence
    # if search_alt_strand=True, then get the reverse complement of the sequence
    if search_alt_strand:
        inseq_use = reverse_complement(inseq)
    else:
        inseq_use = inseq
    inseq_frames = define_frames(inseq=inseq_use,search_start_codons=search_start_codons,search_stop_codons=search_stop_codons)
    # inseq_translated_frames = [translatstopgae_sequence(i) for i in inseq_frames]


    inseq_translated_frames = [translate_sequence(i) for i in inseq_frames]


    # omit polypeptides with more than two stop codons
    inseq_translated_frames = [i.strip("*") for i in inseq_translated_frames if len(re.findall(r'\*',i)) == 1]
    print("\t%d frames found"%len(inseq_translated_frames))
    # print(inseq_translated_frames)
    # inseq_translated_frames = [] # some way to save this as a nested array?
    print("Emitting full frames")
    # write_peptide_frames_old([inseq_translated_frames],filename="internal_frame_translations",name=name+"_peptide_frames",lenfilter=0)
    write_peptide_frames_old([inseq_translated_frames],filename="internal_frame_translations",name=name+"_peptide_frames",lenfilter=8)
    # write_peptide_frames([inseq_translated_frames],name=name+"_peptide_frames",lenfilter=0)
    # write_peptide_frames(inseq_translated_frames,name=name+"_peptide_frames",lenfilter=0)
    if digest_peptides:
        inseq_translated_frames_peptides = []
        print("\tobtaining peptides ..."),
        for i in peptide_lengths:
            # I should indicate the coordinates for each frame in the future...
            a = [[''.join(k) for k in tile_sequence(inseq=j,k=i,overlap=False)] for j in inseq_translated_frames] # this line is throwing an error
            # inseq_translated_frames_peptduleraides.append(a)
            inseq_translated_frames_peptides.append(a)
            # need to allow the fill peptide if len(peptide) < tile
            # print(inseq_translated_frames_peptides)
            write_peptides(a,name=name+"_peptides_"+str(i)+"aa")
            print("\t %d %d-mers ..."%(len(a),i)),
        # return(inseq_translated_frames_peptides)
    # else:
    return(inseq_translated_frames)

"""
# 3/26/2019
Next additions to streamline peptide outputs
1. Several of the frames are overlapping; combine them into a single consensus frame
2. Remove specific frames that are too small (below the peptide length limit)
"""

def main():
    infa = sys.argv[1] # infasta
    print("Now reading %s ..."%infa)
    inseq,inseqname = import_sequence(infa,concat_seq=False)
    infilename = os.path.basename(infa).split(".txt")[0] # infa.split(".fa")[0]
    print("\t ... defining peptides ... "),

    # remove file if it exists
    filename="internal_frame_translations"
    if os.path.exists(filename+".fa"):
        os.remove(filename+".fa")


    ##
    # if we have multiple, then we don't concat
    if type(inseq) is list:
        for i,j in zip(inseq,inseqname):
            j_rename = "_".join(j.split(":"))
            print(j_rename)
            inpeptides = define_peptides(inseq=i,name=j_rename,digest_peptides=False)
    else:
        inpeptides = define_peptides(inseq=inseq,name=infilename,digest_peptides=False)
    ##

    # inpeptides = define_peptides(inseq=inseq,name=infilename,digest_peptides=False)
    print("peptides obtained!")

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
else:
    print("Imported \'propose_nres_peptides\'")

# incorporte a function to just produce a translation and not to
