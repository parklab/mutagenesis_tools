#!/usr/bin/env python

import os,sys,re
import subprocess
import argparse
import pysam
import glob
import itertools
from scipy.stats import chisquare
import pybedtools as pbt
from multiprocessing import Pool

sys.path.insert(0,'/n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/')
import tile_peptides_from_sequence_dev as tilepep


# test object
class mutation:
    # set attributes
    def __init__(self, contig, pos, ref=None, alt=None,**kwargs):
        # check kwargs
        # determine if we should use zero_based
        if 'one_based' in kwargs:
            is_one_based = kwargs['one_based']
        else:
            is_one_based = False
        if type(is_one_based) is not bool:
            is_one_based = False
            
        # determine strand
        if 'strand' in kwargs:
            self.strand = kwargs['strand']
        else:
            self.strand = "+"
            
        # determine metadata
        if 'metadata' in kwargs:
            self.metadata = kwargs['metadata']
        else:
            self.metadata = ['']
        
        # set the coordinates
        self.contig = str(contig) # chromosome/contig where the mutation occurs
        if type(pos) is str or type(pos) is int:
            if is_one_based:
                self.pos = (int(pos),int(pos))
            else:
                pos1 = min(0,int(pos) - 1)
                pos2 = int(pos)
                self.pos = (min(pos1,pos2),max(pos1,pos2))
        elif type(pos) is tuple or type(pos) is list:
            self.pos = pos
        else:
            print("WARNING: pos is not str, int, tuple, or list. Specifying None for the position")
            self.pos = None
        
        # reference and alt variants
        self.ref = ref
        self.alt = alt
        
        # instantiate options for local WT/mutant sequence context
        self.wt_sequence_context = None
        self.mutant_sequence_context = None
        
        
#         # type of variant
#         if self.pos[1] - self.pos[0] <= 1:
#             self.type_of_variant = "snv/snp"
#         elif self.pos[1] - self.pos[0] > 1 and self.pos[1] - self.pos[0] <= 500:
#             self.type_of_variant = "small_indel"
#         elif self.pos[1] - self.pos[0] > 500:
#             self.type_of_variant = "sv"
#         else:
#             self.type_of_variant = None
            
    def __repr__(self):
        # string representation
        mut_info = [self.contig,
                   self.pos[0],
                   self.pos[1],
                   self.ref,
                   self.alt,
                   self.strand,
                   self.wt_sequence_context,
                   self.mutant_sequence_context] + self.metadata
            
        return('\t'.join([str(i) for i in mut_info]))
    
    # set reference genome and genome file
    def set_reference_genome(self,ref_genome,**kwargs):
        self.ref_genome = ref_genome
        
        if 'genome_file' in kwargs:
            self.genome_file = kwargs['genome_file']
        else:
            self.genome_file = None
    
    # set gtf
    def set_gtf(self,gtf):
        self.gtf = gtf
    
    ### operations
    
class mutation_file:
    # this option will store a list of mutation objects.
    # test if we can paralllelize operations involving mutation objects (instantiation of mutation objects, mutagenesis, etc)
    pass