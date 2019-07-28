#!/usr/bin/env python

# import biopython

import os,sys,re
# import regex

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

def search_instances_re(inre,inseq,max_n_substitution):
    # uses regex
    # from https://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string
    regex_re = "(" + inre + "){e<=" + str(max_n_substitution) + "}"
    matches = regex.finditer(regex_re,inseq)
    return(matches)

def check_substitutions(mismatch_site,inseq):
    # LOOP START
    # step 1: make recoding
    # step 2: check translation
    # step 3: determine edit distance from original
    # LOOP END
    # step 4: find the recodings with the smallest edit distance
    # step 4: return outputstring
    return()

