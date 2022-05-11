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
import make_all_nmer_substitutions_copy_dev as nmersub

# test object
# tips on kwargs/args for parallelization: https://stackoverflow.com/questions/54029931/how-to-pass-args-and-kwargs-to-multiprocessing-pool
# tips for multiprocessing with objects: https://stackoverflow.com/questions/43002766/apply-a-method-to-a-list-of-objects-in-parallel-using-multi-processing

# tips on kwargs/args for parallelization: https://stackoverflow.com/questions/54029931/how-to-pass-args-and-kwargs-to-multiprocessing-pool
# tips for multiprocessing with objects: https://stackoverflow.com/questions/43002766/apply-a-method-to-a-list-of-objects-in-parallel-using-multi-processing

# test object
class mutation:
    # set attributes
    def __init__(self, contig, pos, ref=None, alt=None,**kwargs):
        # check kwargs
        # determine if we should use zero_based
        if 'zero_based' in kwargs:
            is_zero_based = kwargs['zero_based']
        else:
            is_zero_based = False
        if type(is_zero_based) is not bool:
            is_zero_based = False

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

        # determine if gtf, reference genome, or genome file is set
        if 'gtf' in kwargs:
            self.gtf = kwargs['gtf']
        else:
            self.gtf = None
        if 'ref_genome' in kwargs:
            self.ref_genome = kwargs['ref_genome']
        else:
            self.ref_genome = None
        if 'genome_file' in kwargs:
            self.genome_file = kwargs['genome_file']
        else:
            self.genome_file = None

        ###

        # set the coordinates
        self.contig = str(contig) # chromosome/contig where the mutation occurs
        if type(pos) is str or type(pos) is int:
            if not is_zero_based:
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

        # shorthand for all coordinates
        self.coordinate = [self.contig] + list(self.pos) + [self.ref,self.alt]
        self.coordinate = tuple([str(i) for i in self.coordinate])

        # instantiate options for local WT/mutant sequence context
        self.wt_sequence_context = self.ref
        self.mutant_sequence_context = self.alt

        # instantiate options for amino acid sequence context
        self.wt_residue_context = None
        self.mutant_residue_context = None
        self.mutant_vs_wt_pairs = []

        # instantiate options for the strand of the annotation
        self.annot_strand = [] # if we have both...

        # set the name
        self.name = '~'.join(self.coordinate)

    def __repr__(self):
        # string representation
        mut_info = [self.contig,
                   self.pos[0],
                   self.pos[1],
                   self.ref,
                   self.alt,
                   self.strand,
                   self.name,
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

    def as_cosmic_context(self):
        # set reverse-complement table
        base_to_rc = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
        if self.ref in ['G','A']:
            print("converting to COSMIC six-channel C/T>N mutation...")
            start_nt_context = ''.join([base_to_rc[j.upper()] for j in self.wt_sequence_context])
            start_nt_context = start_nt_context[::-1]
            end_nt_context = ''.join([base_to_rc[j.upper()] for j in self.mutant_sequence_context])
            end_nt_context = end_nt_context[::-1]
            self.wt_sequence_context = start_nt_context
            self.mutant_sequence_context = end_nt_context

    ### operations
    def mutate_sequence(self,left_slop,right_slop,**kwargs):
        # check if a given sequence is given
        if 'seq' in kwargs:
            seq = kwargs['seq']
            # if seq is provided, then over-ride the self.wt_sequence context
        else:
            seq = None

        if self.wt_sequence_context is None and seq is None:
            return
        elif self.wt_sequence_context is not None and seq is None:
            seq = self.wt_sequence_context

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = False
        if 'strand_aware' in kwargs:
            strand_aware = kwargs['strand_aware']
        else:
            strand_aware = False
        # print('strand awareness:',strand_aware)
        # print('verbose',verbose)

#         right_window = seq[::-1]
        if verbose:
            print(seq,left_slop,right_slop,"seq to mutate")

        # we need a different way to specify the left/right windows
        left_window = seq[:left_slop]
        if self.ref != '-':
            # SNVs, SNPs, or deletions. NOT insertions
            # left_window = seq[:(left_slop - 1)]
            # right_window = seq[(left_slop - 1):] # a correction from before...
            left_window = seq[:left_slop]
            right_window = seq[left_slop:] # a correction from before...
#            right_window = seq[(left_slop - 1 + len(self.alt)):]
        else:
            # the case of insertions
#             left_window = seq[:(left_slop+1)]
            # left_window = seq[:left_slop]
#             right_window = seq[(left_slop+1):]
#             right_window = seq[(left_slop-2):] # to handle how insertions are coded
            left_window = seq[:(left_slop+1)]
            right_window = seq[left_slop:] # to handle how insertions are coded
        ####

        if verbose:
            print("left_window",left_window,seq)
            print("right_window",right_window,seq)

        if self.alt != '-':
            if verbose:
                print(left_window,"join with",self.alt,"and",right_window[len(self.ref):])
#             new_mutant_sequence = left_window + self.alt + right_window[1:] # should cover SNVs or truncations
            new_mutant_sequence = left_window + self.alt + right_window[len(self.ref):] # should cover SNVs or truncations
        else:
            # for deletions
            if verbose:
                print(left_window,"join with",right_window[len(self.ref):])
#             new_mutant_sequence = left_window + right_window[len(self.alt):] # should cover deletions
            new_mutant_sequence = left_window + right_window[len(self.ref):] # should cover deletions

        if verbose:
            print("New mutant sequence",new_mutant_sequence)

        wtseq = seq
        mutseq = new_mutant_sequence
        return(wtseq,mutseq)



    def get_sequence_context(self,**kwargs):
        ### kwargs to get arguments
        # since we require a reference genome, we will check for it. Else, we exit
        if self.ref_genome is None:
            print("No reference genome specified with 'ref_genome', so cannot retrieve sequence")
            return

        # left slop
        if 'leftslop' in kwargs:
            leftslop = int(kwargs['leftslop'])
        else:
            leftslop = 0

        if 'rightslop' in kwargs:
            rightslop = int(kwargs['rightslop'])
        else:
            rightslop = 0

        if 'use_cosmic' in kwargs:
            use_cosmic = kwargs['use_cosmic']
        else:
            use_cosmic = False

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = False

        # strand-aware?
        if 'strand_aware' in kwargs:
            strand_aware = kwargs['strand_aware']
        else:
            strand_aware = False
        # print('strand awareness:',strand_aware)

        # only trigger "equalslop" if it is specified; otherwise, default to left/right slop
        if 'equalslop' in kwargs:
            leftslop = rightslop = int(kwargs['equalslop'])
        if verbose:
            print(leftslop,"leftslop")
            print(rightslop,"rightslop")

        # add an option to return the sequence instead of altering the object attribute
        if 'return_sequence_only' in kwargs:
            return_sequence_only = kwargs['return_sequence_only']
        else:
            return_sequence_only = False
        if type(return_sequence_only) is not bool:
            return_sequence_only = False

        # update positions

        ## start calculation
        # this operation will
        # 1. fish out the 5-3 sequence around the bases while minding the strand and using leftslop/rightslop to adjust the length accordingly
        # 2. If ref/alt bases are specified, then two different sequences will be generated -- a WT sequence, and a mutant sequence
        # 3. This function won't return anything -- it will just update the mutant_sequence_context and wt_sequence_context attributes

        # first, set as bedtool
        # if coordinate is specified, then override the use of the object's own coordinate
        if 'coord' in kwargs:
            coord = kwargs['coord']
        else:
            # coord = [self.contig,
            #          self.pos[0] - leftslop,
            #          self.pos[1] + rightslop + 1,
            #          self.strand]
            coord = [self.contig,
                     self.pos[0] - leftslop - 1,
                     self.pos[1] + rightslop,
                     self.name,
                     '0',
                     self.strand] # because the command requires zero-based...
        if verbose:
            print(coord,"coord")
        coord = ' '.join([str(i) for i in coord])
        # print("coord",coord)
        coord_bt = pbt.BedTool(coord,from_string=True)

        # second, set fasta
        fasta = pbt.example_filename(self.ref_genome)

        # now, read sequence
        a = coord_bt.sequence(fi=fasta,s=False)
        a_fasta_out = open(a.seqfn).read()
        wt_seq = a_fasta_out.split("\n")[1]
        self.wt_sequence_context = wt_seq
        # print("New sequence:",self.wt_sequence_context)

        # if we have specified an alt variant
        # print(self.wt_sequence_context,self.mutant_sequence_context)
        if self.alt is not None:
            # based on the alt, we will either remove the reference base, swap it for the alt base, or insert new bases
            wt,mut = self.mutate_sequence(leftslop,rightslop)
            if return_sequence_only:
                return(wt,mut)
            else:
                self.wt_sequence_context = wt
                self.mutant_sequence_context = mut

        # if strand-aware, then...
        if strand_aware and self.strand == '-':
            # print("printing out negative strand...")
            wt_temp = nmersub.complement_seq(nmersub.reverse_seq(self.wt_sequence_context))
            mut_temp = nmersub.complement_seq(nmersub.reverse_seq(self.mutant_sequence_context))
            # print(self.wt_sequence_context,wt_temp)
            self.wt_sequence_context = wt_temp
            self.mutant_sequence_context = mut_temp

        if use_cosmic:
            self.as_cosmic_context()

    # at adjacent annotations
#     def get_adjacent_annotations(self,leftslop=0,rightslop=0,remove_intersect_file_at_end=True):
    def get_adjacent_annotations(self,leftslop=0,rightslop=0,**kwargs):

        # get requisite scripts
        scriptdir = '/n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/neoantigen_simulation_from_maf/' # replace with a function to get the current script's directory
        adjacent_annot = scriptdir + 'get_adjacent_exon_annotations.sh'

        # remove temp files?
        if 'remove_intersect_file_at_end' in kwargs:
            remove_intersect_file_at_end = bool(kwargs['remove_intersect_file_at_end'])
        else:
            remove_intersect_file_at_end = False # True

        # first, format the mutation position into a bedfile and get the exon(s) which intersect with the mutation
        coord = [self.contig,
                 self.pos[0] - leftslop - 1,
                 self.pos[1] + rightslop,
                 self.strand] # because the command requires zero-based... # this can occasionally get truncated exons not immediately affected by the mutation
        coord = [self.contig,
                 self.pos[0] - 1,
                 self.pos[1],
                 self.strand] # because the command requires zero-based... # this only focuses on which exons touch the mutation

        # the below code will get the exons that intersect with this mutation
        coordentry = ' '.join([str(i) for i in coord])
        coordname = '_'.join([str(i) for i in coord])
        tempcoordbed = '.temp.'+coordname+'.bed'
        tempcoordbedintersect = '.temp.'+coordname+'.intersect_save.bed'
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        coord_bt.intersect(self.gtf,wb=True).saveas(tempcoordbedintersect)

        # second, call the sets of adjacent annotations and read them in
        bashCommand = ['bash',adjacent_annot,tempcoordbedintersect,self.gtf,coordname,'1'] # this line generates all of the files # this will get all annotations
        try:
            output = subprocess.call(bashCommand)
            annot_files = glob.glob('.temp.'+coordname+'.[0-9].bed') # this line recovers all of the beds
        except:
            print("bedtools intersection failed for coordinate",coordentry)

        if remove_intersect_file_at_end:
            print("removing %s"%tempcoordbedintersect)
            rmcommand = ['rm',tempcoordbedintersect]
            subprocess.call(rmcommand)
        return(annot_files)

#    def extract_sequence_at_annotations(self,infile,leftslop,rightslop,remove_infile_at_end=True):
    def extract_sequence_at_annotations(self,infile,leftslop,rightslop,remove_infile_at_end=False):

        # infile is the bedfile of the intervals containing the mutation and the adjacent event(s)
        # steps
        # 1. separate the interval that intersects the mutation from the interval that does not (these should be adjacent)
        # coord
        coord = [self.contig,
                 self.pos[0] - 1,
                 self.pos[1],
                 self.strand] # make into zero-based

        print(coord,infile,"file to use")
        coordentry = ' '.join([str(i) for i in coord])
        print(coordentry,'coord going into the intersect-for-annotation')
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        print(coord_bt)
        print('post bedfile coord going into the intersect-for-annotation')

        # slopped-coord
        slopped_coord = [self.contig,
                 self.pos[0] - leftslop - 1,
                 self.pos[1] + rightslop,
                 self.strand] # because the command requires zero-based...

        # annotation
        annot_bt = pbt.BedTool(infile)
        print(annot_bt)
        print(infile,'annot going into the intersect-for-annotation')
        annot_intersect = annot_bt.intersect(coord_bt,wa=True) # the exon that overlaps with the interval
        annot_adjacent = annot_bt.intersect(coord_bt,wa=True,v=True)

        # let's check why these coordinates failed...


        if len(annot_intersect) == 0 or annot_intersect == ' ' or annot_intersect == '':
            print(annot_intersect,"Error: no intersection with annotation found for %s. Defaulting to genome context..."%self.__repr__())
            print(len(annot_intersect))
            # add the option to get the wt-context from genome here and do the mutation...
            print("===")
            return(None,None)
        else:
            print(infile,'passes the intersect-for-annotation')

        if annot_intersect[0][0] != self.contig:
            print("Error: contigs mismatch")
            return(None,None)
        else:
            print("contigs match %s and %s"%(annot_intersect[0][0],self.contig)) # this works...
        ### MAY 7, 2022: THE ABOVE WORKS; IT IS NOT THE CAUSE FOR WHY SVS GET 'LOST' FROM THE SEQUENCES

        # 2. within the intersecting interval, determine whether the slopped sequence is entirely within the intersecting interval
        get_left_adjacent = 0
        get_right_adjacent = 0
        annot_intersect_strand = annot_intersect[0][4]
        print(annot_intersect_strand,"annot_intersect_strand")

        if slopped_coord[1] < int(annot_intersect[0][1]): # in essence -- is the slopped window smaller than the current exon? if so, then
            print("slopped interval is left-bounded by the current ") # in essence, the slopped-coordinate is broader than the left-most exon, so we need to trim the slopped coordinate on the left-end
            get_left_adjacent = abs(slopped_coord[1] - int(annot_intersect[0][1]))
            # actually, set the strand to the current mutation...
            seq = [slopped_coord[0],str(int(annot_intersect[0][1]) - 1),slopped_coord[2],annot_intersect_strand] # turning this into a zero-based coordinate
            # specifically, this stores the sequence coordinates
        else:
            seq = slopped_coord [:-1] + [annot_intersect_strand] # essentially, keep the strand
        print(seq,'to_use_1')

        if slopped_coord[2] > int(annot_intersect[0][2]):
            print("slopped interval is right-bounded")
            get_right_adjacent = abs(slopped_coord[2] - int(annot_intersect[0][2]))
            seq = [seq[0],seq[1],annot_intersect[0][2],seq[3]] # modify the right-bound accordingly
            #         print(seq,"after controlling for right bounding")
        # if the right-window doesn't exceed the annotation window size, then keep the window as is
        # else:
        #     seq = slopped_coord [:-1] + [annot_intersect_strand] # essentially, keep the strand
        # else:
        #     seq = slopped_coord [:-1] + [We ] # essentially, keep the strand

        print(seq,'to_use_2')

        seq = [str(j) for j in seq]
        print(seq,'to_use_3')

        # 3. For any "overlaps", look for non-intersecting intervals that can come before or after and take those sequences
        seq1 = ['']
        seq2 = ['']

        print("NOW ASSESSING WHICH ADJACENT ANNOTATIONS TO ASSESS")
        print(annot_adjacent)
        for i in annot_adjacent:
            print("********=======********=======********=======")
            print(i)
            print('interval to assess adjacent against',slopped_coord)
            print("********=======********=======********=======********=======********=======********=======")

            # BUG IS SOMEWHERE HERE -- NOT SPECIFYING OVERLAP LENGTH CORRECTLY
            # i is the adjacent exon
#             print(i,slopped_coord,"compare annot with slopped coordinate")
            print(int(i[1]),slopped_coord[1],slopped_coord,get_left_adjacent,'for evaluating left-slop')
            if int(i[1]) < slopped_coord[1] and int(i[2]) < slopped_coord[1] and get_left_adjacent > 0:
                # if the current "adjacent exon" has a start-site before the current slop-coordinate site
                # and if the current "adjacent exon" has an end site before the current slop-coordinate site...

                print("Exon to left-extend into:",i) # this checks out

                # this looks for the upstream adjacent exon if we need that overlap...
#                 seq1 = [i[0],int(i[2]) - get_left_adjacent - 1,i[2]]
                # seq1 = [i[0],int(i[2]) - get_left_adjacent,i[2],i[3]]
                seq1 = [i[0],                                       # contig
                        int(i[2]) - get_left_adjacent + 1,          # include all annotations that are within  'get_left_adjacent' of the current window
                        i[2],
                        i[3]] # does this fix the one-off error?
                # this ends up encompassing all adjacent exons within that particular range

                print(seq1,"seq1")
                seq1 = [str(j) for j in seq1]

            # BUG IS HERE... may 8, 2022:
            # print(i)
            print(i.__repr__(),slopped_coord,'check right-slop') # HERE IS THE BUG -- WHYWe  ARE WE NOT KEEPING THE RIGHT-HAND ANNOTATIONS?
            print(int(i[2]),slopped_coord[2],slopped_coord,get_right_adjacent,'for evaluating right-slop')
            if int(i[1]) > slopped_coord[2] and int(i[2]) > slopped_coord[1] and get_right_adjacent > 0: # this bug!
            # if int(i[2]) > slopped_coord[2] and int(i[2]) > slopped_coord[1] and get_right_adjacent > 0: # this bug!
                # this looks for the upstream adjacent exon if we need that overlap...
                print("Exon to right-extend into:",i)

#                 seq2 = [i[0],int(i[1]) - 1,int(i[1]) + get_right_adjacent]
                seq2 = [i[0],int(i[1]) - 1,int(i[1]) + get_right_adjacent,i[3]]

                print(seq2,"seq2")

                seq2 = [str(j) for j in seq2]

        # 4. Format the collection of sequences into a single sequence and set as the wt_sequence_context and mutant sequences
        seq_total = [seq1,seq,seq2] # this one works okay, though why do none of the
        print(seq_total,"seq_total");
        # get the strand of the annotations...
        annot_strand_to_use = list(set([i[-1] for i in seq_total if i[-1] != '.'] or i[-1] != ''))
        # print(annot_strand_to_use,"annot_strand before filtering...")
        annot_strand_to_use = [i for i in self.annot_strand if i != '']
        print(annot_strand_to_use,"annot_strand")
        self.annot_strand += annot_strand_to_use # list(set([i[-1] for i in seq_total if i[-1] != '.'] or i[-1] != ''))
        # annot_strand_to_use = annot_strand_to_use[0]
        # self.annot_strand = [i for i in self.annot_strand if i != '']

        seq_total = [' '.join(i) for i in seq_total if i != ['']]
        print(seq_total,'sequ_total_to_infer')
        seq_total_fa = ''
        for i in seq_total:
            # coord bedpt
            coord_bt = pbt.BedTool(i,from_string=True)
            print("\t\t\t\t\tprocessing sequence...")
            print(coord_bt,'sequence_to_obtain')
            # set fasta
            fasta = pbt.example_filename(self.ref_genome)
            # now, read sequence
            try:
                a = coord_bt.sequence(fi=fasta) # bug is somewhere here
            except:
                print("cannot extract a sequence for %s. Skipping..."%i)
                continue
            a_fasta_out = open(a.seqfn).read()
            print(a_fasta_out)
            print(a_fasta_out[:10],'sequence obtained by bedtools for for',i)
            wt_seq = a_fasta_out.split("\n")[1]
            # mutant sequence
            seq_total_fa += wt_seq

        # if the strand is negative, get the reverse-complement
        # if annot_strand_to_use[0] == '-':
        # if annot_intersect_strand == '-':
        #     print("getting reverse-complement of sequence...")
        #     seq_total_fa_final_temp = nmersub.reverse_seq(seq_total_fa)
        #     seq_total_fa_final = nmersub.complement_seq(seq_total_fa_final_temp)
        #     print("wt sequence %s is now stored as %s"%(seq_total_fa,seq_total_fa_final))
        # else:
        #     seq_total_fa_final = seq_total_fa

        # Finally, we need to clean up the temporary bed files to prevent over-cluttering if specified
        if remove_infile_at_end:
            print("removing %s"%infile)
            rmcommand = ['rm',infile]
            subprocess.call(rmcommand)
        print(seq_total_fa[1:100],'seq_to_return')
        print(annot_intersect_strand,'annot_intersect_strand')
        return(seq_total_fa,annot_intersect_strand)
        # the parent function 'get_sequence_context_at_annotation' will then run the mutate_sequence operation to simulate mutant


    def extract_sequence_at_annotations_revised(self,infile,leftslop,rightslop,remove_infile_at_end=False):

        # infile is the bedfile of the intervals containing the mutation and the adjacent event(s)
        # steps
        # 1. separate the interval that intersects the mutation from the interval that does not (these should be adjacent)
        # coord
        if self.ref == '-':
            coord = [self.contig,
                     self.pos[0],
                     self.pos[1],
                     '.'] # keep as is for insertions
            # slopped_coord = [self.contig,
            #          self.pos[0] - leftslop,
            #          self.pos[1] + rightslop,
            #          self.strand] # because the command requires zero-based...
            slopped_coord = [self.contig,
                     self.pos[0] - leftslop,
                     self.pos[1] + rightslop,
                     '.'] # because the command requires zero-based...
        else:
            # coord = [self.contig,
            #          self.pos[0] - 1,
            #          self.pos[1],
            #          self.strand] # make into zero-based
            coord = [self.contig,
                     self.pos[0] - 1,
                     self.pos[1],
                     '.'] # make into zero-based
            # slopped_coord = [self.contig,
            #          self.pos[0] - leftslop - 1,
            #          self.pos[1] + rightslop,
            #          self.strand] # because the command requires zero-based...
            slopped_coord = [self.contig,
                     self.pos[0] - leftslop - 1,
                     self.pos[1] + rightslop,
                     '.'] # because the command requires zero-based...



        print(coord,infile,"file to use")
        coordentry = ' '.join([str(i) for i in coord])
        print(coordentry,'coord going into the intersect-for-annotation')
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        print(coord_bt)
        print('post bedfile coord going into the intersect-for-annotation')

        # slopped-coord


        # annotation
        annot_bt = pbt.BedTool(infile)
        print(annot_bt)
        print(infile,'annot going into the intersect-for-annotation')
        annot_intersect = annot_bt.intersect(coord_bt,wa=True) # the exon that overlaps with the interval
        annot_adjacent = annot_bt.intersect(coord_bt,wa=True,v=True)

        # make sure we haven't encountered any annotation failures...
        if len(annot_intersect) == 0 or annot_intersect == ' ' or annot_intersect == '':
            print(annot_intersect,"Error: no intersection with annotation found for %s. Defaulting to genome context..."%self.__repr__())
            print(len(annot_intersect))
            # add the option to get the wt-context from genome here and do the mutation...
            print("===")
            return(None,None)
        else:
            print(infile,'passes the intersect-for-annotation')

        if annot_intersect[0][0] != self.contig:
            print("Error: contigs mismatch")
            return(None,None)
        else:
            print("contigs match %s and %s"%(annot_intersect[0][0],self.contig)) # this works...
        #
        # 2. for the current intersected interval -- if the slopped coordinate is still smaller than the current exon, curtail the sequence
        # note that we have just one annotation file here...
        # set the left-window coordinates
        window_coordinates = [coord] # this is being defined okay
        left_stop = False
        if int(slopped_coord[1]) > int(annot_intersect[0][1]):
            print('the left hand side of the current intersected exon is sufficient')
            window_coordinates = [[slopped_coord[0],slopped_coord[1],coord[1],annot_intersect[0][4]]] + window_coordinates
            # note that slopped_coord is 0-based
            left_stop = True
        else:
            # note that annot_intersect is 1-based
            window_coordinates = [[slopped_coord[0],int(annot_intersect[0][1]) - 1,coord[1],annot_intersect[0][4]]] + window_coordinates # add the rest of the left-hand side of the window
            # BUG HERE -- I have to go one base above...
            budget = leftslop - (int(coord[1]) - int(annot_intersect[0][1]))
#            budget = leftslop - (int(coord[1]) - (int(annot_intersect[0][1]) - 1))
            # what coordinates do we need
            left_intervals_to_store = []
            left_intervals_to_use = [i for i in annot_adjacent if int(i[2]) < int(annot_intersect[0][1])]
            # now, loop...
            counter = 0
            while budget > 0 and len(left_intervals_to_use) > 0: # while we still have a budget...
                print('current conditions',counter,budget,len(left_intervals_to_use),infile)
            # while budget > 0 or len(left_intervals_to_use) > 0: # while we still have a budget...
                # print('current set of left intervals:',left_intervals_to_use,'budget so far',budget)
                current_intvl_len = int(left_intervals_to_use[-1][2]) - int(left_intervals_to_use[-1][1])
                print('current set of left intervals:',left_intervals_to_use,'budget so far',budget,'current interval len',current_intvl_len)
                if current_intvl_len <= budget:
                    coord_to_store = [left_intervals_to_use[-1][0],
                                      str(int(left_intervals_to_use[-1][1]) - 1), # we need to one-base the coordinates
                                      str(left_intervals_to_use[-1][2]),
                                      left_intervals_to_use[-1][4]]
                    budget -= current_intvl_len
                else:
                    coord_to_store = [left_intervals_to_use[-1][0],
                                      str(int(left_intervals_to_use[-1][2]) - budget - 1),
                                      str(left_intervals_to_use[-1][2]),
                                      left_intervals_to_use[-1][4]]
                    print("this current left-interval will spend our budget. Coord to store:",coord_to_store)
                    # budget -= budget
                    budget = 0
                left_intervals_to_store.append(coord_to_store)
                left_intervals_to_use = left_intervals_to_use[:-1]
                counter += 1
            left_stop = True
            left_intervals_to_store = left_intervals_to_store[::-1]
            window_coordinates = left_intervals_to_store + window_coordinates

        # set the right-window coordinates
        right_stop = False
        print('right-hand side of slopped coordinate:',slopped_coord[2])
        print('right-hand side of annotation:',annot_intersect[0][2])
        if int(slopped_coord[2]) < int(annot_intersect[0][2]):
            print('the right hand side of the current intersected exon is sufficient')
            if self.ref == '-':
                window_coordinates = window_coordinates + [[slopped_coord[0],coord[2]-1,slopped_coord[2],annot_intersect[0][4]]]
            else:
                window_coordinates = window_coordinates + [[slopped_coord[0],coord[2],slopped_coord[2],annot_intersect[0][4]]]
            right_stop = True
        else:
            if self.ref == '-':
                window_coordinates = window_coordinates + [[slopped_coord[0],coord[2]-1,annot_intersect[0][2],annot_intersect[0][4]]] # add the rest of the left-hand side of the window
            else:
                window_coordinates = window_coordinates + [[slopped_coord[0],coord[2],annot_intersect[0][2],annot_intersect[0][4]]] # add the rest of the left-hand side of the window
            budget = rightslop - (int(annot_intersect[0][2]) - int(coord[2]))
            # what coordinates do we need
            right_intervals_to_store = []
            right_intervals_to_use = [i for i in annot_adjacent if int(i[1]) > int(annot_intersect[0][2])]
            # now, loop...
            while budget > 0 and len(right_intervals_to_use) > 0: # while we still have a budget...
            # while budget > 0 or len(right_intervals_to_use) > 0: # while we still have a budget...
                current_intvl_len = int(right_intervals_to_use[-1][2]) - int(right_intervals_to_use[-1][1])
                print('current set of right intervals:',right_intervals_to_use,'budget so far',budget,'current interval len',current_intvl_len)
                if current_intvl_len <= budget:
                    coord_to_store = [right_intervals_to_use[-1][0],
                                      str(right_intervals_to_use[-1][1]),
                                      str(right_intervals_to_use[-1][2]),
                                      right_intervals_to_use[-1][4]]
                    budget -= current_intvl_len
                else:
                    coord_to_store = [right_intervals_to_use[-1][0],
                                      str(int(right_intervals_to_use[-1][1])),
                                      str(int(right_intervals_to_use[-1][1]) + budget),
                                      right_intervals_to_use[-1][4]]
                    print("this current right-interval will spend our budget. Coord to store:",coord_to_store)
                    # budget -= budget
                    budget = 0
                right_intervals_to_store.append(coord_to_store)
                right_intervals_to_use = right_intervals_to_use[:-1]
            right_stop = True
            window_coordinates = window_coordinates + right_intervals_to_store

        if right_stop and left_stop:
            seq_total = window_coordinates
        print(seq_total,"seq_total");

        # conert everything to str
        seq_total = [[str(j) for j in i] for i in seq_total]
        # convert the coord to str so that we can filter out appropriately...

        # get the strand of the annotations...
        annot_strand_to_use = list(set([i[-1] for i in seq_total if i[-1] != '.'] or i[-1] != ''))
        print(annot_strand_to_use,"annot_strand before filtering...")
        # annot_strand_to_use = [i for i in self.annot_strand if i != '']
        annot_intersect_strand = annot_strand_to_use[0]
        print(annot_intersect_strand,"annot_strand")
        self.annot_strand += annot_intersect_strand # list(set([i[-1] for i in seq_total if i[-1] != '.'] or i[-1] != ''))
        # annot_strand_to_use = annot_strand_to_use[0]
        # self.annot_strand = [i for i in self.annot_strand if i != '']

        seq_total = [' '.join(i) for i in seq_total if i != ['']]
        coord_mutation_str = ' '.join([str(j) for j in coord]) # store the coordinate of the mutation here...

        print(seq_total,'sequ_total_to_infer')
        seq_total_fa = ''
        for i in seq_total:
            # coord bedpt
            if i != coord_mutation_str:
                coord_bt = pbt.BedTool(i,from_string=True)
                print("\t\t\t\t\tprocessing sequence...")
                print(coord_bt,'sequence_to_obtain')
                # set fasta
                fasta = pbt.example_filename(self.ref_genome)
                # now, read sequence
                try:
                    a = coord_bt.sequence(fi=fasta) # bug is somewhere here
                except:
                    print("cannot extract a sequence for %s. Skipping..."%i)
                    continue
                a_fasta_out = open(a.seqfn).read()
                print(a_fasta_out)
                print(a_fasta_out[:10],'sequence obtained by bedtools for for',i)
                wt_seq = a_fasta_out.split("\n")[1]
                # mutant sequence
                seq_total_fa += wt_seq
            else:
                seq_total_fa += self.ref
        # generate the mutant sequence here...
        mut_seq_total_fa = ''
        for i in seq_total:
            # coord bedpt
            if i != coord_mutation_str:
                coord_bt = pbt.BedTool(i,from_string=True)
                print("\t\t\t\t\tprocessing sequence...")
                print(coord_bt,'sequence_to_obtain')
                # set fasta
                fasta = pbt.example_filename(self.ref_genome)
                # now, read sequence
                try:
                    a = coord_bt.sequence(fi=fasta) # bug is somewhere here
                except:
                    print("cannot extract a sequence for %s. Skipping..."%i)
                    continue
                a_fasta_out = open(a.seqfn).read()
                print(a_fasta_out)
                print(a_fasta_out[:10],'sequence obtained by bedtools for for',i)
                wt_seq = a_fasta_out.split("\n")[1]
                # mutant sequence
                mut_seq_total_fa += wt_seq
            else:
                mut_seq_total_fa += self.alt
        # remove any dashes or dots
        seq_total_fa = ''.join([i for i in seq_total_fa if i != '-'])
        mut_seq_total_fa = ''.join([i for i in mut_seq_total_fa if i != '-'])


        # Finally, we need to clean up the temporary bed files to prevent over-cluttering if specified
        if remove_infile_at_end:
            print("removing %s"%infile)
            rmcommand = ['rm',infile]
            subprocess.call(rmcommand)
        print(seq_total_fa[1:100],'seq_to_return')
        print(mut_seq_total_fa[1:100],'mut_seq_to_return')
        print(annot_intersect_strand,'annot_intersect_strand')
        return(seq_total_fa,mut_seq_total_fa,annot_intersect_strand)





    def get_sequence_context_at_annotation(self,**kwargs):
        # NOTE that this only works for the first adjacent feature on either side
        # since we require a GTF, we will check for it. Else, we exit
        if self.gtf is None:
            print("No GTF specified, so cannot retrieve sequence at annotation. Perhaps try get_sequence_context to get the sequence context without annotations")
            return
        if self.ref_genome is None:
            print("No reference genome specified, so cannot retrieve sequence")
            return

        # should we translate?
        if 'get_translation' in kwargs:
            get_translation = bool(kwargs['get_translation'])
        else:
            get_translation = True

        # this operation does the same as 'get_sequence_context', except it will use exonic sequences to decide
        if 'leftslop' in kwargs:
            leftslop = kwargs['leftslop']
        else:
            leftslop = 0

        if 'rightslop' in kwargs:
            rightslop = kwargs['rightslop']
        else:
            rightslop = 0

        # only trigger "equalslop" if it is specified; otherwise, default to left/right slop
        if 'equalslop' in kwargs:
            leftslop = kwargs['equalslop']
            rightslop = kwargs['equalslop']

        # do we return sequences or set as attributes
        if 'return_sequence_only' in kwargs:
            return_sequence_only = kwargs['return_sequence_only']
        else:
            return_sequence_only = False
        if type(return_sequence_only) is not bool:
            return_sequence_only = False

        # this operation will
        # 1. fish out the 5-3 sequence around the bases while minding the strand and using leftslop/rightslop to adjust the length accordingly
        # 2. Since annotation is required, if the slop position exceeds the current annotated unit, the adjacent annotated units overlapping the window will also be searched
        # 3. If ref/alt bases are specified, then two different sequences will be generated
        # 4. This function won't return anything -- it will just update the mutant_sequence_context and wt_sequence_context attributes

        # 1. get adjacent annotations
        adjacent_annot = self.get_adjacent_annotations(leftslop,rightslop,**kwargs)
        adjacent_annot.sort() # add a sort term in case...
        print(adjacent_annot,"adjacent_annot") # THIS WORKS

        # exit()

        # 1 and 2. get the adjacent sequences...
        wt_running = ''
        mut_running = ''
        seq_contexts = []
        for i in adjacent_annot:
            print(i)
            print(i,"file to run")
            with open(i,'r') as f:
                for j in f:
                    print(j,i,'line from file')
            print("==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888==========888888888888")
            # inseq,strand_seq = self.extract_sequence_at_annotations(i,leftslop=leftslop,rightslop=rightslop) # CHECK THIS...
            inseq,mutseq,strand_seq = self.extract_sequence_at_annotations_revised(i,leftslop=leftslop,rightslop=rightslop) # CHECK THIS revision
            print(self.__repr__(),inseq,'obtained sequence -- is it none') # BUG HAS OCCURRED BY HERE -- NO ANNOTATION SEQUENCE FOUND
            print(self.strand,'obtained sequence -- is it none') # BUG HAS OCCURRED BY HERE -- NO ANNOTATION SEQUENCE FOUND
            print(strand_seq,'obtained sequence -- is it none') # BUG HAS OCCURRED BY HERE -- NO ANNOTATION SEQUENCE FOUND
            print("==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777==========777777777777")
            # if inseq is None and self.wt_sequence_context is not None:
            if inseq is None:
                print(self,"no annotation sequence found...")
            else:
                # # 3. fish out the wt sequence and the mutated wt
                # print('found',inseq,'sequence for mutation',self.__repr__())
                # i_wt,i_mut = self.mutate_sequence(seq=inseq,
                #                                   left_slop=leftslop,
                #                                   right_slop=rightslop)
                # print('created',i_mut,'sequence for mutation',self.__repr__()) #
                if '-' in strand_seq:
                    print("getting reverse-complement of sequences...",inseq,mutseq)
                    i_wt_temp = nmersub.reverse_seq(inseq)
                    inseq = nmersub.complement_seq(i_wt_temp)
                    i_mut_temp = nmersub.reverse_seq(mutseq)
                    mutseq = nmersub.complement_seq(i_mut_temp)

                    print("wt and mut sequences are now stored as %s and %s"%(inseq,mutseq))
                # seq_contexts.append((i_wt,i_mut))
                # may 10, 2022: already adding the sequences...
                seq_contexts.append((inseq,mutseq))
        # get unique sequences
        print(seq_contexts)
        seq_contexts = set(seq_contexts)
        print(seq_contexts)
        for i,j in seq_contexts:
            wt_running += i + ","
            mut_running += j + ","
        # add to wt contexts...
        wt_running = wt_running.strip(",")
        mut_running = mut_running.strip(",")

        if return_sequence_only:
            return(wt_running,mut_running)
        else:
            self.wt_sequence_context = wt_running
            self.mutant_sequence_context = mut_running

        # reduce the annot strands to uniques...
        self.annot_strand = list(set(self.annot_strand))

        # if translate...
        if get_translation:
            self.generate_translation(**kwargs)


        # generate translation
    def generate_translation(self,shear_translations=True,compare_mutant_wt=True,clip_at_stop_codons=False,reset_mutant_vs_wt_pairs=True,store_different_only=True,min_pep_len=7,**kwargs):
        wt_peptides_full = []
        mut_peptides_full = []


        # here, generate the strand-aware sequences
        for i in itertools.zip_longest(self.wt_sequence_context.split(','),
                                       self.mutant_sequence_context.split(','), fillvalue=''):

            print("sequences to mutate:",i,"for",self.__str__())

            # wt sequence
            try:
                wt_peptides,wt_frame_pos,wt_frames = tilepep.return_peptides(inseq=str(i[0]))
            except:
                if len(i[1]) >= 3:
                    print("encountering an error with",i[0])
                else:
                    print("too short",i[0])
                wt_peptides,wt_frame_pos,wt_frames = [],[],[]

            # mutant sequence
            try:
                mut_peptides,mut_frame_pos,mut_frames = tilepep.return_peptides(inseq=str(i[1]))
            except:
                if len(i[1]) >= 3:
                    print("encountering an error with",i[1])
                else:
                    print("too short",i[1])
                mut_peptides,mut_frame_pos,mut_frames = [],[],[]


            # if len(i[0]) >= 3:
            #     print("Now translating",i[0],"the wild type sequence")
            #
            # else:
            #     wt_peptides,wt_frame_pos,wt_frames = [],[],[]
            #
            # if len(i[1]) >= 3:
            #     print("Now translating",i[0],"the mutant sequence")
            #     mut_peptides,mut_frame_pos,mut_frames = tilepep.return_peptides(inseq=i[1])
            # else:
            #     mut_peptides,mut_frame_pos,mut_frames = [],[],[]

            # do we store the full frames or the sheared translations
            if shear_translations:
                wt_peptides_full += wt_peptides
                mut_peptides_full += mut_peptides
            else:
                wt_peptides_full += wt_frames
                mut_peptides_full += mut_frames

        self.mutant_residue_context = mut_peptides_full
        self.wt_residue_context = wt_peptides_full

        if compare_mutant_wt:
            if reset_mutant_vs_wt_pairs:
                self.mutant_vs_wt_pairs = []
            # assemble pairs
            # for i,j in itertools.zip_longest(mut_peptides_full,wt_peptides_full, fillvalue=''):
            for i,j in zip(mut_peptides_full,wt_peptides_full):
                peplen_list = list(itertools.zip_longest(i,j, fillvalue=''))
                for a,b in peplen_list:
                    # if we have to clip at stopcodons, then...
                    print(a,b)
                    if clip_at_stop_codons:
                        a1 = a.split('*')[0]
                        b1 = b.split('*')[0]
                    else:
                        a1 = a
                        b1 = b
                    n_differences = sum([1 for k in list(itertools.zip_longest(a1,b1, fillvalue='')) if k[0] != k[1]])
                    print((a1,b1,a1==b1,n_differences))

                    if n_differences == 0 and store_different_only:
                        continue

                    if len(a1) < min_pep_len and len(b1) < min_pep_len:
                        continue

                    # store unique pairings
                    if (a1,b1) not in self.mutant_vs_wt_pairs:
                        self.mutant_vs_wt_pairs.append((a1,b1,n_differences))

                    # obj.mutant_vs_wt_pairs.append((a1,b1,a1==b1,n_differences))


####################


class mutation_set():
    def __init__ (self,bedfile,**kwargs):

        # bedfile assumes that we have ref/alt positions as columns 4/5 (respectively)
        # we assume strand at position 6
        # we assume other metadata at other columns


        self.mutations = []
        bedobj = pbt.BedTool(bedfile)

        for i in bedobj:
            icontig = i[0]
            ipos = (int(i[1]),int(i[2]))
            iref = i[3]
            ialt = i[4]
            istrand = i[5]
            imeta = i[6:]
            i_new = mutation(contig=icontig,
                             pos=ipos,
                             ref=iref,
                             alt=ialt,
                             strand = istrand,
                             metadata = imeta,
                             **kwargs)
            self.mutations.append(i_new)

    def __iter__(self):
        self.n = 0
        return self

    def __next__(self):
        if self.n < len(self.mutations):
            outmut = self.mutations[self.n]
            self.n += 1
        else:
            self.n = 0
            raise StopIteration
        return(outmut)

    # retrieve a particular index of a mutation
    def get_mutation_at_index(self,ind):
        return(self.mutations[min(ind,len(self.mutations) - 1)])

    # retrieve the mutant/WT context of each position
    def worker_seq_context(self,obj):
        mutobj = obj[0]
        kwargs = obj[1]
        mutobj.get_sequence_context(**kwargs)
        # note: it is not possible to keep the modification lasting when using pool.map according to https://stackoverflow.com/questions/15857838/modify-object-in-python-multiprocessing, but we will instead return the modified object
        return(mutobj)

    def worker_seq_context_annot(self,obj):
        mutobj = obj[0]
        kwargs = obj[1]
        mutobj.get_sequence_context_at_annotation(**kwargs)
        # note: it is not possible to keep the modification lasting when using pool.map according to https://stackoverflow.com/questions/15857838/modify-object-in-python-multiprocessing, but we will instead return the modified object
        return(mutobj)


    def get_seq_context(self,**kwargs):
        # pool
        if 'n_jobs' in kwargs:
            n_jobs = int(kwargs['n_jobs'])
        else:
            n_jobs = 10
        print("Running %d jobs at a time in parallel "%n_jobs)
        pool = Pool(n_jobs)
        args = [kwargs] * len(self.mutations)
        input_zip = zip(self.mutations,args)
        new_mutations = pool.map(self.worker_seq_context,
                              input_zip)
        pool.close()
        pool.join()
        # update the set of mutations with the copies of these mutations but with modified wt/mutant sequence contexts
        self.mutations = new_mutations


    def get_seq_context_at_annot(self,**kwargs):
        # a variant of get_seq_context, except it searches the annotation-defined sequence
        # pool
        if 'n_jobs' in kwargs:
            n_jobs = int(kwargs['n_jobs'])
        else:
            n_jobs = 10
        print("Running %d jobs at a time in parallel "%n_jobs)
        pool = Pool(n_jobs)
        args = [kwargs] * len(self.mutations)
        input_zip = zip(self.mutations,args)
        new_mutations = pool.map(self.worker_seq_context_annot,
                              input_zip)
        pool.close()
        pool.join()

        # update the set of mutations with the copies of these mutations but with modified wt/mutant sequence contexts
        self.mutations = new_mutations


    def get_stats_areas(self,**kwargs):
        all_areas = self.calculate_areas(**kwargs)
        return(max(all_areas),min(all_areas),np.mean(all_areas),np.median(all_areas))
