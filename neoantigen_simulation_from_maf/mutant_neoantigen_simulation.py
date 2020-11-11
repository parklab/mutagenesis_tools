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
        
#         right_window = self.wt_sequence_context[::-1]
        right_window = seq[::-1]

        if self.ref != '-':
#             left_window = self.wt_sequence_context[:(left_slop)]
            left_window = seq[:(left_slop)]
            right_window = right_window[:right_slop]
        else:
            # the case of insertions
#             left_window = self.wt_sequence_context[:(left_slop+1)]
            left_window = seq[:(left_slop+1)]
            right_window = right_window[:(right_slop+1)]
        right_window = right_window[::-1]
            
        if self.alt != '-':
            new_mutant_sequence = left_window + self.alt + right_window # should cover SNVs or truncations
        else:
            new_mutant_sequence = left_window + right_window # should cover deletions
        
#         wtseq = self.wt_sequence_context
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
            coord = [self.contig,
                     self.pos[0] - leftslop,
                     self.pos[1] + rightslop + 1,
                     self.strand]
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
        coord = [self.contig,
                 self.pos[0] - leftslop,
                 self.pos[1] + rightslop + 1,
                 self.strand]
        coordentry = ' '.join([str(i) for i in coord])
        coordname = '_'.join([str(i) for i in coord])
        tempcoordbed = '.temp.'+coordname+'.bed'
        tempcoordbedintersect = '.temp.'+coordname+'.intersect_save.bed'
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        coord_bt.intersect(self.gtf,wb=True).saveas(tempcoordbedintersect)
        
        # second, call the sets of adjacent annotations and read them in
        print(' '.join([adjacent_annot,tempcoordbedintersect,self.gtf,coordname]))
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
        coord = [self.contig,
                 self.pos[0],
                 self.pos[1] + 1,
                 self.strand]
        coordentry = ' '.join([str(i) for i in coord])
        coord_bt = pbt.BedTool(coordentry,from_string=True)
        
        # slopped-coord
        slopped_coord = [self.contig,
                 self.pos[0] - leftslop,
                 self.pos[1] + rightslop + 1,
                 self.strand]

        # annotation
        annot_bt = pbt.BedTool(infile)
        annot_intersect = annot_bt.intersect(coord_bt,wa=True)
        annot_adjacent = annot_bt.intersect(coord_bt,wa=True,v=True)
        
        # reject if contigs mismatch for any reason
        if annot_intersect[0][0] != self.contig:
            print("Error: contigs mismatch")
            return
        else:
            print("contigs match %s and %s"%(annot_intersect[0][0],self.contig))

        # 2. within the intersecting interval, determine whether the slopped sequence is entirely within the intersecting interval
        get_left_adjacent = 0
        get_right_adjacent = 0
        annot_intersect_strand = annot_intersect[0][4]
        if slopped_coord[1] < int(annot_intersect[0][1]):
            print("slopped interval is left-bounded")
            get_left_adjacent = abs(slopped_coord[1] - int(annot_intersect[0][1]))
            seq = [slopped_coord[0],annot_intersect[0][1],slopped_coord[2],annot_intersect_strand]
        else:
            seq = slopped_coord[:-1] + [annot_intersect_strand]
            
        if slopped_coord[2] > int(annot_intersect[0][2]):
            print("slopped interval is right-bounded")
            get_right_adjacent = abs(slopped_coord[2] - int(annot_intersect[0][2]))
            seq = [seq[0],seq[1],annot_intersect[0][2],seq[3]]
        
        seq = [str(j) for j in seq]            
        
        # 3. For any "overlaps", look for non-intersecting intervals that can come before or after and take those sequences
        for i in annot_adjacent:
            if int(i[1]) < slopped_coord[1] and int(i[2]) < slopped_coord[1] and get_left_adjacent > 0:
                # this looks for the upstream adjacent exon if we need that overlap...
                seq1 = [i[0],int(i[2]) - get_left_adjacent,i[2]]
                seq1 = [str(j) for j in seq1]
            else:
                seq1 = ['']
            if int(i[1]) > slopped_coord[2] and int(i[2]) < slopped_coord[1] and get_right_adjacent > 0:
                # this looks for the upstream adjacent exon if we need that overlap...
                seq2 = [i[0],i[1],int(i[1]) + get_right_adjacent]
                seq2 = [str(j) for j in seq2]
            else:
                seq2 = ['']
        
        # 4. Format the collection of sequences into a single sequence and set as the wt_sequence_context and mutant sequences
        seq_total = [seq1,seq,seq2]
        seq_total = [' '.join(i) for i in seq_total if i != ['']]
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
        
        # 1 and 2. get the adjacent sequences...
        wt_running = ''
        mut_running = ''
        seq_contexts = []
        for i in adjacent_annot:
            inseq = self.extract_sequence_at_annotations(i,leftslop=leftslop,rightslop=rightslop)
            # 3. fish out the wt sequence and the mutated wt
            i_wt,i_mut = self.mutate_sequence(seq=inseq,
                                              left_slop=leftslop,
                                              right_slop=rightslop)
            seq_contexts.append((i_wt,i_mut))
        # get unique sequences
        seq_contexts = set(seq_contexts)
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
            
    # need to add an option to generate peptides from the sequences
    def generate_translation(shear_translations=True,
                             compare_mutant_wt=True,
                             clip_at_stop_codons=False,
                             reset_mutant_vs_wt_pairs=True,
                             **kwargs):
        # mutant sequence
        mut_peptides,mut_frame_pos,mut_frames = tilepep.return_peptides(inseq=obj.mutant_sequence_context,**kwargs)
    
        # wt sequence
        wt_peptides,wt_frame_pos,wt_frames = tilepep.return_peptides(inseq=obj.wt_sequence_context,**kwargs)
        
        # do we store the full frames or the sheared translations
        if shear_translations:
            obj.mutant_residue_context = mut_peptides
            obj.wt_residue_context = wt_peptides
        else:
            obj.mutant_residue_context = mut_frames
            obj.wt_residue_context = wt_frames
                    
        if compare_mutant_wt:
            if reset_mutant_vs_wt_pairs:
                obj.mutant_vs_wt_pairs = []
            # assemble pairs
            for i,j in zip(mut_peptides,wt_peptides):
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
                    obj.mutant_vs_wt_pairs.append((a1,b1,a1==b1,n_differences))
####################
    
        
class mutation_set():
    
    # this option will store a list of mutation objects.
    # test if we can paralllelize operations involving mutation objects (instantiation of mutation objects, mutagenesis, etc)
    # future options will
    # 1. Conduct mutational signature analysis
    # 2. Look for coverage stats across all mutations
    # 3. Do annotation-burden tests
    # 4. Eventually accomodate diverse mutations type (CNs?)
    
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
    
    ## workers
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

    # deployment functions
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
        
    # translation
    ## workers
    def worker_generate_translations(self,obj):
        mutobj = obj[0]
        kwargs = obj[1]
        mutobj.generate_translation(**kwargs)
        # note: it is not possible to keep the modification lasting when using pool.map according to https://stackoverflow.com/questions/15857838/modify-object-in-python-multiprocessing, but we will instead return the modified object
        return(mutobj)

    ## deployment
    def get_translations_at_context(self,**kwargs):
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
        new_mutations = pool.map(self.worker_generate_translations, 
                              input_zip)
        pool.close()
        pool.join()
        
        # update the set of mutations with the copies of these mutations but with modified wt/mutant sequence contexts
        self.mutations = new_mutations

    
########