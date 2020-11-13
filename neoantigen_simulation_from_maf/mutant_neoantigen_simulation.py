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
        
#         right_window = seq[::-1]
        
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
        
        
        print("left_window",left_window,seq)
        print("right_window",right_window,seq)
            
        if self.alt != '-':
            print(left_window,"join with",self.alt,"and",right_window[len(self.ref):])
#             new_mutant_sequence = left_window + self.alt + right_window[1:] # should cover SNVs or truncations
            new_mutant_sequence = left_window + self.alt + right_window[len(self.ref):] # should cover SNVs or truncations
        else:
            # for deletions
            print(left_window,"join with",right_window[len(self.ref):])
#             new_mutant_sequence = left_window + right_window[len(self.alt):] # should cover deletions
            new_mutant_sequence = left_window + right_window[len(self.ref):] # should cover deletions
            
        print("New mutant sequence",new_mutant_sequence)
        
        wtseq = seq
        mutseq = new_mutant_sequence
        return(wtseq,mutseq)
    
    
    
    def get_sequence_context(self,**kwargs):
        ### kwargs to get arguments
        # since we require a reference genome, we will check for it. Else, we exit
        if self.ref_genome is None:
            print("No reference genome specified, so cannot retrieve sequence")
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
        
        # only trigger "equalslop" if it is specified; otherwise, default to left/right slop
        if 'equalslop' in kwargs:
            leftslop = rightslop = int(kwargs['equalslop'])
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
                     self.strand] # because the command requires zero-based...

        print(coord,"coord")
        coord = ' '.join([str(i) for i in coord])
        coord_bt = pbt.BedTool(coord,from_string=True)
        
        # second, set fasta
        fasta = pbt.example_filename(self.ref_genome)
        
        # now, read sequence
        a = coord_bt.sequence(fi=fasta)
        a_fasta_out = open(a.seqfn).read()
        wt_seq = a_fasta_out.split("\n")[1]
        self.wt_sequence_context = wt_seq
        
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
    
    # at adjacent annotations
    def get_adjacent_annotations(self,leftslop=0,rightslop=0,remove_intersect_file_at_end=True):
        # get requisite scripts
        scriptdir = '/n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis/neoantigen_simulation_from_maf/' # replace with a function to get the current script's directory
        adjacent_annot = scriptdir + 'get_adjacent_exon_annotations.sh'
        
        # first, format the mutation position into a bedfile and calculate the intersection
        # coord = [self.contig,
        #          self.pos[0] - leftslop,
        #          self.pos[1] + rightslop + 1,
        #          self.strand]
        coord = [self.contig,
                 self.pos[0] - leftslop - 1,
                 self.pos[1] + rightslop,
                 self.strand] # because the command requires zero-based...

        coordentry = ' '.join([str(i) for i in coord])
        coordname = '_'.join([str(i) for i in coord])
        tempcoordbed = '.temp.'+coordname+'.bed'
        tempcoordbedintersect = '.temp.'+coordname+'.intersect_save.bed'
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        coord_bt.intersect(self.gtf,wb=True).saveas(tempcoordbedintersect)
        
        # second, call the sets of adjacent annotations and read them in
#         print(' '.join([adjacent_annot,tempcoordbedintersect,self.gtf,coordname]))
        bashCommand = ['bash',adjacent_annot,tempcoordbedintersect,self.gtf,coordname]
        output = subprocess.call(bashCommand)
        annot_files = glob.glob('.temp.'+coordname+'.[0-9].bed')
        
        if remove_intersect_file_at_end:
            print("removing %s"%tempcoordbedintersect)
            rmcommand = ['rm',tempcoordbedintersect]
            subprocess.call(rmcommand)
        return(annot_files)
    
    def extract_sequence_at_annotations(self,infile,leftslop,rightslop,remove_infile_at_end=True):
        # infile is the bedfile of the intervals containing the mutation and the adjacent event(s) 
        # steps
        # 1. separate the interval that intersects the mutation from the interval that does not (these should be adjacent)
        # coord
        # coord = [self.contig,
        #          self.pos[0],
        #          self.pos[1] + 1,
        #          self.strand]
        coord = [self.contig,
                 self.pos[0] - 1,
                 self.pos[1],
                 self.strand] # make into zero-based

        print(coord)
        coordentry = ' '.join([str(i) for i in coord])
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        
        # slopped-coord
        # slopped_coord = [self.contig,
        #          self.pos[0] - leftslop,
        #          self.pos[1] + rightslop + 1,
        #          self.strand]
        slopped_coord = [self.contig,
                 self.pos[0] - leftslop - 1,
                 self.pos[1] + rightslop,
                 self.strand] # because the command requires zero-based...

        # annotation
        annot_bt = pbt.BedTool(infile)
        annot_intersect = annot_bt.intersect(coord_bt,wa=True)
        annot_adjacent = annot_bt.intersect(coord_bt,wa=True,v=True)
        
        # reject if contigs mismatch for any reason
#         print(annot_intersect,len(annot_intersect),"annot_intersect")
        
#         if len(annot_adjacent) == 0 or annot_intersect == ' ' or annot_intersect == '':
#             print("===")
#             print(annot_adjacent,"Error: no intersection with annotation found. Defaulting to genome context...")
#             print(len(annot_adjacent))
#             print("===")
#             return
        if len(annot_intersect) == 0 or annot_intersect == ' ' or annot_intersect == '':
            print(annot_intersect,"Error: no intersection with annotation found for %s. Defaulting to genome context..."%self.__repr__())
            print(len(annot_intersect))
            # add the option to get the wt-context from genome here and do the mutation...
            print("===")
            return

        
        if annot_intersect[0][0] != self.contig:
            print("Error: contigs mismatch")
            return
        else:
            print("contigs match %s and %s"%(annot_intersect[0][0],self.contig)) # this works...

        # 2. within the intersecting interval, determine whether the slopped sequence is entirely within the intersecting interval
        get_left_adjacent = 0
        get_right_adjacent = 0
        annot_intersect_strand = annot_intersect[0][4]
        print(annot_intersect_strand,"annot_intersect_strand")
#         self.annot_strand.append(annot_intersect_strand)
        
#         print(slopped_coord,annot_intersect[0],"slop -vs- intersected annot")
        if slopped_coord[1] < int(annot_intersect[0][1]):
            print("slopped interval is left-bounded")
            get_left_adjacent = abs(slopped_coord[1] - int(annot_intersect[0][1]))
            # actually, set the strand to the current mutation...
            # seq = [slopped_coord[0],annot_intersect[0][1],slopped_coord[2],annot_intersect_strand]
            seq = [slopped_coord[0],str(int(annot_intersect[0][1]) - 1),slopped_coord[2],annot_intersect_strand] # turning this into a zero-based coordinate
#             seq = [slopped_coord[0],annot_intersect[0][1],slopped_coord[2],slopped_coord[3]]
        else:
#             print(slopped_coord)
#             seq = slopped_coord[:-1] + [annot_intersect_strand] # essentially, keep the strand
            seq = slopped_coord [:-1] + [annot_intersect_strand] # essentially, keep the strand
#         print(seq,"after controlling for left bounding")
            
        if slopped_coord[2] > int(annot_intersect[0][2]):
#             print("slopped interval is right-bounded")
            get_right_adjacent = abs(slopped_coord[2] - int(annot_intersect[0][2]))
            seq = [seq[0],seq[1],annot_intersect[0][2],seq[3]] # modify the right-bound accordingly
#         print(seq,"after controlling for right bounding")
        
        seq = [str(j) for j in seq]            
        
        # 3. For any "overlaps", look for non-intersecting intervals that can come before or after and take those sequences
        seq1 = ['']
        seq2 = ['']
        
        for i in annot_adjacent:
            
            # BUG IS SOMEWHERE HERE -- NOT SPECIFYING OVERLAP LENGTH CORRECTLY
            # i is the adjacent exon
#             print(i,slopped_coord,"compare annot with slopped coordinate")
            if int(i[1]) < slopped_coord[1] and int(i[2]) < slopped_coord[1] and get_left_adjacent > 0:
                # if the current "adjacent exon" has a start-site before the current slop-coordinate site
                # and if the current "adjacent exon" has an end site before the current slop-coordinate site...
                
                print("Exon to left-extend into:",i)
                
                # this looks for the upstream adjacent exon if we need that overlap...
#                 seq1 = [i[0],int(i[2]) - get_left_adjacent - 1,i[2]]
                # seq1 = [i[0],int(i[2]) - get_left_adjacent,i[2],i[3]]
                seq1 = [i[0],int(i[2]) - get_left_adjacent + 1,i[2],i[3]] # does this fix the one-off error?

                print(seq1,"seq1")
                seq1 = [str(j) for j in seq1]

            # BUG IS HERE...
            if int(i[1]) > slopped_coord[2] and int(i[2]) > slopped_coord[1] and get_right_adjacent > 0: # this bug!
                # this looks for the upstream adjacent exon if we need that overlap...
                print("Exon to right-extend into:",i)
                
#                 seq2 = [i[0],int(i[1]) - 1,int(i[1]) + get_right_adjacent]
                seq2 = [i[0],int(i[1]) - 1,int(i[1]) + get_right_adjacent,i[3]]

                print(seq2,"seq2")

                seq2 = [str(j) for j in seq2]
        
        # 4. Format the collection of sequences into a single sequence and set as the wt_sequence_context and mutant sequences
        seq_total = [seq1,seq,seq2]
        # get the strand of the annotations...
        self.annot_strand += list(set([i[-1] for i in seq_total if i[-1] != '.'] or i[-1] != ''))
        self.annot_strand = [i for i in self.annot_strand if i != '']

        seq_total = [' '.join(i) for i in seq_total if i != ['']]
        print(seq_total)
        seq_total_fa = ''
        for i in seq_total:
            # coord bedpt
            coord_bt = pbt.BedTool(i,from_string=True)
            # set fasta
            fasta = pbt.example_filename(self.ref_genome)
            # now, read sequence
            a = coord_bt.sequence(fi=fasta)
            a_fasta_out = open(a.seqfn).read()
            wt_seq = a_fasta_out.split("\n")[1]
            # mutant sequence       
            seq_total_fa += wt_seq

        # Finally, we need to clean up the temporary bed files to prevent over-cluttering if specified
        if remove_infile_at_end:
            print("removing %s"%infile)
            rmcommand = ['rm',infile]
            subprocess.call(rmcommand)
        return(seq_total_fa)
        
        
        # the parent function 'get_sequence_context_at_annotation' will then run the mutate_sequence operation to simulate mutant
        
    
    
    def get_sequence_context_at_annotation(self,**kwargs):
        # NOTE that this only works for the first adjacent feature on either side
        # since we require a GTF, we will check for it. Else, we exit
        if self.gtf is None:
            print("No GTF specified, so cannot retrieve sequence at annotation. Perhaps try get_sequence_context to get the sequence context without annotations")
            return
        if self.ref_genome is None:
            print("No reference genome specified, so cannot retrieve sequence")
            return
        
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
        adjacent_annot = self.get_adjacent_annotations(leftslop,rightslop)
        print(adjacent_annot,"adjacent_annot") # THIS WORKS
        
        # 1 and 2. get the adjacent sequences...
        wt_running = ''
        mut_running = ''
        seq_contexts = []
        for i in adjacent_annot:
            inseq = self.extract_sequence_at_annotations(i,leftslop=leftslop,rightslop=rightslop) # CHECK THIS...
            print(inseq)
            if inseq is None:
                print(self,"no annotation sequence found...")
            # 3. fish out the wt sequence and the mutated wt
            i_wt,i_mut = self.mutate_sequence(seq=inseq,
                                              left_slop=leftslop,
                                              right_slop=rightslop)
            seq_contexts.append((i_wt,i_mut))
        # get unique sequences
        print(seq_contexts)
        seq_contexts = set(seq_contexts)
        print(seq_contexts)
        for i,j in seq_contexts:
            wt_running += i_wt + ","
            mut_running += i_mut + ","
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