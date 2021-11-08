import os,sys,regex
import numpy,math,random
import itertools
import numpy
from collections import defaultdict
import argparse


def import_sequence(infa,concat_seq=True,assess_rc=False):
	inseqname = []
	inseq = []
	print(infa)
    # for i in infile:
	current_seq = ''
	with open(infa,'r') as infafile:
		for i in infafile:
			if i[0] != ">":
				current_seq += i.strip('\n')
			else:
				# add the header
				inseqname.append(i.strip("\n")[1:])
				if current_seq != '':
					inseq.append(current_seq)
					current_seq = ''
		inseq.append(current_seq)
	if concat_seq:
		seq_full = ''.join(inseq)
		seqname_full = ';'.join(inseqname)
	else:
		seq_full = inseq
		seqname_full = inseqname
    # if assess_rc, then also get the RC of each sequence
	if assess_rc:
		inseq_rc = []
		inseqname_rc = []
		for i,j in zip(inseqname,inseq):
			inseqname_rc.append("reverse_complement_"+i)
			inseq_rc.append(get_rc_internal(j))
		inseq += inseq_rc
		inseqname += inseqname_rc
	return(seq_full,seqname_full)


def reverse_seq(inseq):
	revseq = [i for i in inseq]
	revseq.reverse()
	revseq = ''.join(revseq)
	return(revseq)

def complement_seq(inseq):
	compl = {'A':'T','C':'G','T':'A','G':'C'}
	compseq = [compl[i.upper()] for i in inseq]
	return(compseq)

def get_new_substitutions (inseq,make_deletions=False):

	bases = ['a','c','g','t']
	if make_deletions:
		bases += '-'
	inds_bases = range(len(bases))
	kmers = []

	#for i in itertools.combinations_with_replacement(bases, len(inseq)):
	for i in itertools.product(inds_bases, repeat = len(inseq)):
		combination = [bases[k] for k in i]
		kmers.append(''.join(combination))

	substitutions = list(set(kmers) - set([inseq.lower()]))
	return(substitutions)


def tile_sequence(inseq,k,overlap=False):
	# elegant solution: https://stackoverflow.com/questions/2485669/consecutive-overlapping-subsets-of-array-numpy-python?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
	as_strided = numpy.lib.stride_tricks.as_strided
	#b = as_strided(a, (11,4), a.strides*2)

	ntiles = max(len(inseq) - k + 1,1) # if the sequence is too small to be tiled, keep the sequence

	seq = numpy.asarray([i for i in inseq])
	#inds = numpy.arange(0,len(seq))
	tiles = as_strided(seq, (ntiles,k), seq.strides*2) # how to keep regions that aren't the same size as the tile?
	# print(len(tiles))

	if not overlap:
		pick_tiles = range(0,len(tiles),k)
		# print(pick_tiles)
		new_tiles = [tiles[i] for i in pick_tiles]
		# excludes the bases that did not end up in the tile

		# if there are extra bases that don't fit in the tile...
	# 	if pick_tiles[-1] + k > len(tiles):
	# 		extra_bases = inseq[pick_tiles[-1]+1:len(tiles)]
	# 		new_tiles.append(extra_bases.split())
	# 	new_tiles = numpy.asarray(new_tiles)
	# 	print(len(new_tiles))
	else:
		new_tiles = tiles
	# we should add any "leftover" sequences
	return(new_tiles)
	# return(tiles)

# def build_sequences (seq_before,seq_after,mut_seqs):
# 	tiles = []
# 	return()

def mutate_gene (inseq,k,indel=0):
	# generate all non-overlapping 'k'-base tiles of the sequence 'inseq'
	# preallocate the size of the sequences. NOTE that we want to generate the sequence with only the "tile" being the mutant
	# for every tile, generate the "mutated" portion and the "constant" portion. Make the mutant sequences and ligate to the constant portion
	tiles = tile_sequence(inseq,k)
	indices = tile_sequence(numpy.arange(0,len(inseq)),k)
	seqind = numpy.arange(0,len(inseq))

	# list for mutated sequences
	mutated_sequences = []
	# go through the tiles
	for i,j in zip(tiles,indices):
		indices_before = [inseq[k] for k in list(seqind) if k not in list(j) and k < j[0]]
		seq_before = ''.join(indices_before)
		indices_after = [inseq[k] for k in list(seqind) if k not in list(j) and k > j[-1]]
		seq_after = ''.join(indices_after)

		# get substitutions
		substitutions = get_new_substitutions(''.join(i),make_deletions=indel > 0)
		for k in substitutions:
			new_seq = seq_before + k + seq_after
			new_seq = new_seq.replace('-','')
			mutated_sequences.append(new_seq)


	# should we make indels?
	if indel > 0:
		print("TO MAKE INDELS")
		indel_sequences = []
		# make all possible insertions



	return(mutated_sequences)

def count_mutation_types(wt,probes):
	syn = 0
	nsyn = 0
	nonsense = 0
	for i in probes:
		if "*" in i:
			nonsense += 1
		elif i == wt:
			syn += 1
		elif i != wt:
			nsyn += 1
	print(("synonymous","nonsynonymous","nonsense"))
	return((syn,nsyn,nonsense))

def codon_table():
	bases = ['t', 'c', 'a', 'g']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))
	return(codon_table)

def reverse_codon_table ():
	in_codon_table = codon_table()
	# set the dictionary of codons
	rev_codon_table = defaultdict(list)
	for i in in_codon_table:
		aa = in_codon_table[i]
		rev_codon_table[aa].append(i)
	return(rev_codon_table)

def translation(inseq):
	codons = codon_table()
	aas = []
	for i in range(0,len(inseq),3):
		test_codon = ''.join(inseq[i:i+3]).lower()
		if len(test_codon) == 3:
			aas.append(codons[test_codon])
	return(''.join(aas))

def get_codons(AA):
	rct = reverse_codon_table()
	return(rct[AA.upper()])

def mass_translations(seqs):
	# this function is unique in that it will produce a unique set of translated sequences without translating every single sequence
	translations = [None] * len(seqs)
	nt_seqs = [None] * len(seqs)
	for n,i in enumerate(seqs):
		tl = translation(i)
		if tl not in translations:
			translations[n] = tl
			nt_seqs[n] = i
	translations = [i for i in translations if i is not None]
	nt_seqs = [i for i in nt_seqs if i is not None]
	return(translations,nt_seqs)

def get_gc_content(seq):
	#regex
	gcre = regex.compile('[GCgc]')
	gcpctg = float(len(gcre.findall(seq)))/float(len(seq))
	return(gcpctg)

def get_repeat_content(seq):
	# in this function, check how many times a particular sequence repeats
	return()

def get_edit_distance(seq,substitution_dict):
	# based on the substitution matrix, determine how different the codons are.
	return()

def get_codon_set_with_gclim(seq,gclim,markovian=True):
	# print("Input sequence",seq)
	# initialize the set of codons to use
	codons = [None] * len(range(0,len(seq),3))
	# the codon sequence
	codon_seq = []
	# the starting position: currpos = 0
	currpos = 0
	# start codons are all those that could produce the AA at the starting positions
	start_codons = get_codons(seq[currpos])
	# for each codon, check whether it is close to the gc lim
	currpos = 0 # initialize the current position
	for currpos in range(0,len(seq)):
		# get all the codons at the current position
		curr_codons = get_codons(seq[currpos])
		# the centerpiece of the method for translation selection:
		# get the codon sequence with the smallest difference in GC content with the previous codon
		# the problem with this method is that it gets influenced by the earlier codons
		# if the function is specified to be "markovian=True," then we will compare just the current codon's GC content to the GC limit
		# if the function is specified to be "markovian=False," then we will compare the current codon's GC content and all of the previous sequence to the GC limit
		if markovian:
			gc_at_diff = list(numpy.absolute([get_gc_content(i) - gclim for i in curr_codons]))
		else:
			compare_seq = [''.join(codon_seq+[i]) for i in curr_codons]
			# print('compare seq',compare_seq)
			gc_at_diff = list(numpy.absolute([get_gc_content(i) - gclim for i in compare_seq]))

		selected_codon = curr_codons[gc_at_diff.index(min(gc_at_diff))]
		codon_seq+=[selected_codon]

	return(''.join(codon_seq))


def get_codon_set(seq,markovian=True):
	"""
	# this will generate all possible combinations of codons that produce the amino acid sequence
	# then, paths that go through the codon sets will be chosen that
	# a. maximize the edit distance between each codon
	# b. maintain the GC-AT percentage at relatively even levels (give or take a certain tolerance)
	# c. avoid STOP codons
	# select a path using a dynamic programming algorithm: pick a random codon to start with, and select the path that meets the criteria
	# repeat for each starting codon and cycle through for n steps until you find a common path through the codons
	# pick the "common" path
	# 1. generate the codon table
	# reverse_codon_table
	# 2. pick a random starting codon. Say the first one
	# start_codon = codons[0][0]
	# 3. the algorithm:
	# for i in amino acid sequence:
	# 	from present position, calculate the edit distance to each codon in the subsequent base and the combined GC-AT content
	# 	find the codon with the greatest edit distance and the smallest difference in GC-AT content
	#	if you find ties, pick one path at random
	#	STORE: present AA and next AA
	# continue the for-loop until you reach the END of the AA sequence
	"""
	# initialize the set of codons to use
	codons = [None] * len(range(0,len(seq),3))
	# the codon sequence
	codon_seq = []
	# the starting position: currpos = 0
	currpos = 0
	# start codons are all those that could produce the AA at the starting positions
	start_codons = get_codons(seq[currpos])
	# pick a random codon to begin
	# for now, pick the first of these codons as the next codon
	#curr_codon = start_codons[random.randint(0,len(start_codons))]
	curr_codon = start_codons[0]
	codon_seq.append(curr_codon)
	# the next AA position that we will go to is at 1
	nextpos = 1
	while nextpos < len(seq) - 1:
		# get all the codons at the next codon set
		next_codons = get_codons(seq[nextpos])
		#print(next_codons)
		#edit_dist = [edit_dist(i) for i in codons]
		# either use the most recent codon or the full preceding sequence to perform the comparison
		if markovian:
			compare_seq = curr_codon
		else:
			compare_seq = ''.join(codon_seq)
		# gc_at_diff = list(numpy.absolute([get_gc_content(i) - get_gc_content(compare_seq) for i in next_codons]))
		gc_at_diff = list(numpy.absolute([get_gc_content(i) - get_gc_content(compare_seq) for i in next_codons]))
		#print(gc_at_diff)
		# get the position "i" with the max edit_dist
		# get the position "j" with the min edit_dist
		# if i = j, then move to the codon
		# codon_seq.append(chosen_codon)

		# the centerpiece of the method for translation selection:
		# get the codon sequence with the smallest difference in GC content with the previous codon
		# the problem with this method is that it gets influenced by the earlier codons
		codon_seq.append(next_codons[gc_at_diff.index(min(gc_at_diff))])

		# the desired method would consider if the possible codon creates a repetitive sequence.
		# if you had a nucleotide transition/transversion matrix, you could compute a probability of observing the previous codon becoming the next (one-to-one correspondence of amino acids)

		# If the codons are also sufficiently distinct from one another in sequence content (edit distance between )

		currpos += 1
		nextpos += 1
		# curr_codon = chosen_codon
	# repeat this function for n steps. Then, get a common codon set
	return(''.join(codon_seq))

def get_codon_set_traversal(seq,seqref,fract_mismatch_tol=0.1,markovian=True):
	"""
	# This function will generate the sets of DNA sequences found within the reference sequence that yield the peptide upon translation.
	# Additional features to allow for (1) sequence composition constraints and (2) mismatch tolerances will be added at a later time
	# STEP 1: INITIALIZATION
	# START: a list of the codons for the starting amino acid. This is the "running sequence"
	# for each amino acid "i" thereafter:
	# 	generate all possible codon sets
	#	propose a combination of a running sequence and a current codon
	#	keep the combination if it is found within the reference sequence
	# terminate once the full peptide as been traversed
	# return the sets of "running sequences"
	# This method is a greedy algorithm and is probably not the most efficient one! Maybe propose random sequences to initialize the set and update its composition based on which sequences are found and which are not?
	"""
	# initialize the set of codons to use
	codons = [None] * len(range(0,len(seq),3))
	# initialize the empty list that will hold the codon sequence
	# codon_seq will be a "list of lists". We will treat the possible codon
	codon_seq = []
	# the starting position: currpos = 0
	# currpos = 0
	# start codons are all those that could produce the AA at the starting positions
	# start_codons = get_codons(seq[currpos])
	#curr_codon = start_codons[random.randint(0,len(start_codons))]
	# curr_codon = start_codons[0]
	# the next AA position that we will go to is at 1
	nextpos = 0
	# nbases_mismatch = int(len(seq) * fract_mismatch_tol)
	while nextpos < len(seq) - 1:
		print("Position %d. Amino acid %s"%(nextpos,seq[nextpos]))
		# we do not want more than half of the base so far to exhibit a mismatch
		nbases_mismatch = int((nextpos + 1) * 3 * fract_mismatch_tol)
		print(nbases_mismatch)
		# get all the codons at the next codon set
		next_codons = get_codons(seq[nextpos])
		print("Next codons:",next_codons)
		# create all possible codon sequences
		if len(codon_seq) > 0:
			proposed_seq = [i+j for j in next_codons for i in codon_seq]
			print("Proposed sequence",proposed_seq)
			# keep the members of proposed_seq if they are found within the seqref
			found_seq = [i for i in proposed_seq if len(regex.findall("("+i+"){s<="+str(nbases_mismatch)+"}",seqref,regex.IGNORECASE)) > 0]
			codon_seq = found_seq
			print("Ongoing codon sequence",codon_seq)
		else:
			codon_seq = next_codons
# 			print("start codon!")
		nextpos += 1
# 		print("---")
		# curr_codon = chosen_codon
	# repeat this function for n steps. Then, get a common codon set
	return((seq,[''.join(i) for i in codon_seq]))

def generate_codon_set (inseq,gclim,steps=1000):
	codons = []
	for i in range(0,steps,1):
		codons.append(get_codon_set(inseq,gclim))

def get_diff_codon_set_simple(aaseq):
	# 1. generate all possible NT sequences of aaset
	# 2. for each sequence, get the
	#	a. substitution distance for the string of codons in frame (assume the same frame) #THIS IS NOT NECESSARY
	# 	b. gc content
	#	c. repeat content
	#	d. cpg content
	#	e. whether the sequence matches the recoded region for shRNA targetting
	# 4. Find the sequence with the greatest edit distance, lowest difference in AT-GC content, lowest repeat content, and lowest CpG content
	# return the sequence
	# this function requires a nucleotide substitution matrix
	return()


def write_fasta(inseq,linelen=50,filename="outfile",otherinfo=''):
	lines = [inseq[i:i+linelen] for i in range(0,len(inseq),linelen)]
	outfile = open(filename+".fa",'w')
	outfile.write(">"+filename+"|"+otherinfo+"\n")
	for i in lines:
		outfile.write(i+"\n")
	outfile.close()


def main():
	parser = argparse.ArgumentParser(description='Mutate a sequence')
	parser.add_argument('-f',
						'--fasta',
						help='FASTA file containing nucleotide sequence(s). REQUIRED',
						type=str)
	parser.add_argument('--indel',
						help='Specify this flag if you also want indels, followed by a number indicating how big of an indel you want at most',
						type=int,
						default=0)
	parser.add_argument('--no_snv',
						dest='snv',
						help='Specify this flag if you do NOT want to produce SNV mutants',
						action='store_false')
	parser.add_argument('-t',
						'--translate',
						help='Specify this flag if you want to translate the mutant sequences',
						action='store_false')
	parser.add_argument('-k',
						'--kmer',
						help='K-mer size. This will make all possible SNVs within each overlapping k-sized subsequence per sequence.',
						type=int,
						default=1)
	parser.add_argument('-r',
						'--revcomp',
						help='Specify this flag if you want to also mutate the reverse complement of your sequences.',
						action='store_true')
	parser.add_argument('-o',
						'--outfile',
						help='Name of the output file. FASTA extension will be automatically appended',
						type=str,
						default='outfile')
	parser.add_argument('-s',
						'--strand',
						help='What is the strand of your sequence? DEFAULT: + strand',
						type=str,
						default='+')
	parser.add_argument('--reduce_representation',
						help="specify this flag if you want to constrain the representation of your sequences",
						action='store_true')
	parser.add_argument('--markovian',
						help='If you only want to constrain sequence representation based on the state of the most recent sequences',
						action='store_true')
	#
	args = parser.parse_args()
	print(args)

	# default arguments
	infasta = args.fasta # sys.argv[1]
	kmer_size = args.kmer # 3
	strand = args.strand # "+"
	markovian = args.markovian # sys.argv[4]

	seq_name = '>test'
	outfile = args.outfile+'.fa'

	with open(outfile,'w') as g:
		g.write('')

	sequences,sequence_names = import_sequence(infasta,concat_seq=False)

	for in_sequence,sequence_name in zip(sequences,sequence_names):
		# strand conversion
		if strand == "-":
			print("negative strand")
			in_sequence = reverse_seq(complement_seq(in_sequence))

		print("Sequence to mutate: %s"%in_sequence)
		# generate mutants
		print("Generating mutants")
		mutants = mutate_gene(in_sequence,kmer_size,args.indel)
		# write as FASTA
		with open(outfile,'a') as g:
			for i,j in enumerate(mutants):
				g.write('>'+sequence_name+"_mutant_"+str(i)+'\n')
				g.write(j+'\n')

if __name__ == "__main__":
	main()
else:
	print("Imported \'make_all_nmer_substitutions\'")
