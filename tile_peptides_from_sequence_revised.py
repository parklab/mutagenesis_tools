import os,sys,re
import math,random
import itertools
import numpy as np
from collections import defaultdict
import argparse

# sys.path.insert(")
# sys.path.insert(0,"/n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis")
sys.path.insert(0,"/n/data1/hms/dbmi/park/vinay/pipelines/mutagenesis")
import make_all_nmer_substitutions_copy_dev as nmersub
# from make_all_nmer_substitutions import nmersub.reverse_seq,nmersub.translation,codon_table,tile_sequence,complement_seq

def get_rc(inseq):
    complements = {'A':'T','T':'A','C':'G','G':'C'}
    return(''.join([complements[i.upper()] if i.upper() == i else complements[i.upper()].lower() for i in inseq[::-1]]))

def import_sequence(infa,concat_seq=True,assess_rc=False):
    inseqname = []
    inseq = []
    # for i in infile:
    current_seq = ''
    with open(infa) as infafile:
        for i in infafile:
            if i[0] != ">":
                # NOTE: we may have multiple fasta lines
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

    print(inseqname)
    print(inseq)

    # if assess_rc, then also get the RC of each sequence
    if assess_rc:
        print("generating reverse complement of input sequences")
        inseq_rc = []
        inseqname_rc = []
        for i,j in zip(inseqname,inseq):
            # print(i,j)
            inseqname_rc.append("reverse_complement_"+i)
            inseq_rc.append(get_rc(j))

        inseq += inseq_rc
        inseqname += inseqname_rc

    print(inseqname)
    print(inseq)

    return(seq_full,seqname_full)


# def import_sequence(infa,concat_seq=True):
#     inseqname = []
#     inseq = []
#     # infile = open(infa)
#     # for i in infile:
#     with open(infa) as infafile:
#         for i in infafile:
#             if i[0] != ">":
#                 inseq.append(i.strip("\n"))
#             else:
#                 inseqname.append(i.strip("\n")[1:])
#     if concat_seq:
#         seq_full = ''.join(inseq)
#         seqname_full = ';'.join(inseqname)
#     else:
#         seq_full = inseq
#         seqname_full = inseqname
#     return(seq_full,seqname_full)
#
def codon_table():
	bases = ['t', 'c', 'a', 'g']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	out_codon_table = dict(zip(codons, amino_acids))
	return(out_codon_table)

def reverse_complement(inseq):
    return(nmersub.reverse_seq(nmersub.complement_seq(inseq)))

def translate_sequence(inseq):
    return(nmersub.translation(inseq))

def start_stop_codons():
    codon_table_in = codon_table()
    start_codon = [i for i in codon_table_in if codon_table_in[i] == "M"]
    stop_codon = [i for i in codon_table_in if codon_table_in[i] == "*"]
    return(start_codon,stop_codon)

def define_frames(inseq,search_start_codons=True,search_stop_codons=True,start_seq_site_lim=3):
    # april 14, 2021: debug this
    # we will save the start and stop sites
    start_seq_sites = []
    end_seq_sites = []
    frames = []
    start_codon,stop_codons = start_stop_codons()
    # print(start_codon,stop_codons)
    ## codons = nmersub.nmersub.nmersub.tile_sequence(inseq = inseq,k=3,overlap=False)
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
    # print(start_seq_sites,stop_seq_sites)
    # print("start-stop")


    # note that if we aren't looking for frames but would like to find an in-frame nmersub.translation, adjust the start and stop sites accordingly so that we can get appropriate frames
    if not search_start_codons and not search_stop_codons:
        # in the case of direct translation
        # other_start = ((stop_seq_sites[0] - start_seq_sites[1]) % 3) - 1
        other_start = ((stop_seq_sites[0] - start_seq_sites[0]) % 3) - 1
        print(other_start)
        # other_start = max(0,other_start) # make sure it starts at 0!
        other_start = min(0,abs(other_start)) # make sure it starts at 0!
        # other_end = stop_seq_sites[0] + ((stop_seq_sites[0] - start_seq_sites[1]) % 3) - 1
        other_end = stop_seq_sites[0] + ((stop_seq_sites[0] - start_seq_sites[0]) % 3) - 1
        print(other_start,other_end)
        # start_seq_sites = start_seq_sites + [other_start,other_start+1,other_start+2]
        # stop_seq_sites = stop_seq_sites + [other_end-2,other_end-1,other_end]
        start_seq_sites = [other_start,other_start+1,other_start+2]
        print(start_seq_sites)
        stop_seq_sites = [other_end-2,other_end-1,other_end]
        print(stop_seq_sites)

    # print(start_seq_sites,stop_seq_sites)

    # print(start_seq_sites,stop_seq_sites,"repeat")


    # search every combination of start and stop codons; keep frames that are divisible by 3.
    frame_positions = [(i,j) for i in start_seq_sites for j in stop_seq_sites if (j - i) % 3 == 0 and j > i]
    # print(frame_positions,"frame")
    # print(frame_positions)
    frames = [inseq[i[0]:i[1]] for i in frame_positions]
    # remove empty frames
    frames = [i for i in frames if len(i) > 0]
    # print(frames,"frames")

    return(frames,frame_positions)

def write_peptides(peptide_array,name,clip_after_stop):
    # outfile = open(os.getcwd()+"/"+name+".txt",'w')
    print("clip after stop:",clip_after_stop)
    outfile = open(name+".txt",'w')
    for i in peptide_array:
        for j in i:
            if clip_after_stop and j[0] != ">":
                # if we specify to clip a peptide after a "*" codon...
                j_clipped = j.split('*')[0]
                outfile.write(j_clipped+"\n")
            else:
                outfile.write(j+"\n")
    outfile.close()

def length_fiter(peptide,minlen=8):
    return(len(peptide) >= minlen)

def write_decoy(in_peptide,randomize=False):
    if randomize:
        return(''.join(random.sample(in_peptide,len(in_peptide))))
    else:
        return(in_peptide[::-1])

def write_peptide_frames_old(peptide_array,name,filename,clip_after_stop,lenfilter=10,write_decoy_peptide=True):
    # write as a FASTA file
    # outfile = open(filename+".fa",'w')
    print(filename)
    outfile = open(filename+".fa",'a') # modified to append
    count = 0
    for i in peptide_array:
        # outfile.write(">"+name+"_"+str(count)+"\n")
        for j in i:
            if clip_after_stop and j[0] != ">":
                j_clipped = j.split('*')[0]
                if len(j_clipped) >= lenfilter:
                    outfile.write(">"+name+"_"+str(count)+"\n")
                    outfile.write(j_clipped+"\n")
                    # if write a decoy, then also provide a decoy peptide
                    if write_decoy_peptide:
                        outfile.write(">decoy_"+name+"_"+str(count)+"\n")
                        decoy_j = write_decoy(j_clipped)
                        outfile.write(decoy_j+"\n")
            else:
                # ADD OPTION HERE TO BREAK THE SEQUENCE INTO SMALLER FRAMES AND NAME EACH ONE "SUBFRAME"
                # HOW MANY STOPS?
                # nstops=j.count('*')
                for k in enumerate(j.split('*')):
                    print(k)
                    k_subframe_count = k[0]
                    k_subframe_str = k[1]
                    if len(k_subframe_str) >= lenfilter:
                        outfile.write(">"+name+"_"+str(count)+"_subframe_"+str(k_subframe_count)+"\n")
                        outfile.write(k_subframe_str+"\n")
                        if write_decoy_peptide:
                            outfile.write(">decoy_"+name+"_"+str(count)+"_subframe_"+str(k_subframe_count)+"\n")
                            # outfile.write(">decoy_"+name+"_"+str(count)+"\n")
                            decoy_k = write_decoy(k_subframe_str)
                            outfile.write(decoy_k+"\n")
            count += 1
    outfile.close()

def write_peptide_frames(peptide_array,name,lenfilter=10):
    # create a separate folder
    if not os.path.exists(name):
        print("writing to new folder")
        os.makedirs(name)
    # write as a FASTA file
    outfile = open(name+"/"+name+".fa",'w')
    count = 0
    for i in peptide_array:
        outfile.write(">"+name+"_"+str(count)+"\n")
        if len(j) >= lenfilter:
            outfile.write(j+"\n")
        count += 1
    outfile.close()


def define_peptides(inseq,name,filename="sequence_translations",peptide_lengths=[8,9,10,11],search_start_codons=False,search_stop_codons=False,search_alt_strand=False,digest_peptides=True,omit_termin=True,clip_after_stop=True):
    # if search_start_codons=True, then get all possible start sites--any position with an M. Else, start from the beginning of the sequence
    # if search_stop_codons=True, then get all possible end sites--any codon with a stop codon. Else, end at the beginning of the sequence
    # if search_alt_strand=True, then get the reverse complement of the sequence
    if search_alt_strand:
        inseq_use = reverse_complement(inseq)
    else:
        inseq_use = inseq
    inseq_frames,inseq_frames_coords = define_frames(inseq=inseq_use,search_start_codons=search_start_codons,search_stop_codons=search_stop_codons) # define trames

    # inseq_translated_frames = [translatstopgae_sequence(i) for i in inseq_frames]

    # April 15, 2021: error is at this line 215 here:
    # print(inseq_frames)
    # print(inseq_frames[0])
    # for i in inseq_frames:
    #     print(i,"frame 1")
    inseq_translated_frames = [translate_sequence(i) for i in inseq_frames]


    # omit polypeptides with more than two stop codons
    if omit_termin:
        inseq_translated_frames = [i.strip("*") for i in inseq_translated_frames if len(re.findall(r'\*',i)) == 1]
    # print("\t%d frames found"%len(inseq_translated_frames)) # suppress
    # print(inseq_translated_frames)
    # inseq_translated_frames = [] # some way to save this as a nested array?
    # print("Emitting full frames") # suppress

    print("clip after stop",clip_after_stop)
    write_peptide_frames_old([inseq_translated_frames],name=name+"_peptide_frames",filename=filename,lenfilter=8,clip_after_stop=clip_after_stop) # at least 8 AA for digestion... # modify to write to the same file!
    # write_peptide_frames(inseq_translated_frames,name=name+"_peptide_frames",lenfilter=0)
    if digest_peptides:
        inseq_translated_frames_peptides = []
        print("\tobtaining peptides ..."),
        for i in peptide_lengths:
            # I should indicate the coordinates for each frame in the future...
            a = [[''.join(i_tile) for i_tile in nmersub.tile_sequence(inseq=j,k=i,overlap=False)] for j in inseq_translated_frames] # this line is throwing an error
            # inseq_translated_frames_peptduleraides.append(a)
            inseq_translated_frames_peptides.append(a)
            # need to allow the fill peptide if len(peptide) < tile
            # print(inseq_translated_frames_peptides)
            write_peptides(a,name=name+"_peptides_"+str(i)+"aa",clip_after_stop=clip_after_stop) # modify to write to the same file!
            print("\t %d %d-mers ..."%(len(a),i)),
        # return(inseq_translated_frames_peptides)
    # else:
    return(inseq_translated_frames)

def return_peptides(inseq,peptide_lengths=[8,9,10,11],search_start_codons=False,search_stop_codons=False,search_alt_strand=False,digest_peptides=True,omit_termin=False,clip_after_stop=True):
    # if search_start_codons=True, then get all possible start sites--any position with an M. Else, start from the beginning of the sequence
    # if search_stop_codons=True, then get all possible end sites--any codon with a stop codon. Else, end at the beginning of the sequence
    # if search_alt_strand=True, then get the reverse complement of the sequence
    if search_alt_strand:
        inseq_use = reverse_complement(inseq)
    else:
        inseq_use = inseq
    # print(inseq_use)
    inseq_frames,inseq_frames_pos = define_frames(inseq=inseq_use,search_start_codons=search_start_codons,search_stop_codons=search_stop_codons)

    # inseq_translated_frames = [translatstopgae_sequence(i) for i in inseq_frames]

    print(inseq_frames) # this is still sorted -- we want to keep the native frame...
    inseq_translated_frames = [translate_sequence(i) for i in inseq_frames]


    # omit polypeptides with more than two stop codons
    if omit_termin:
        inseq_translated_frames = [i.strip("*") for i in inseq_translated_frames if len(re.findall(r'\*',i)) == 1]

    print("frames")
    print(inseq_translated_frames)
    # print("\t%d frames found"%len(inseq_translated_frames)) # suppress
    # print(inseq_translated_frames)
    # inseq_translated_frames = [] # some way to save this as a nested array?
    # print("Emitting full frames") # suppress

    # write_peptide_frames_old([inseq_translated_frames],name=name+"_peptide_frames",filename=filename,lenfilter=8,clip_after_stop=clip_after_stop) # at least 8 AA for digestion... # modify to write to the same file!
    # write_peptide_frames(inseq_translated_frames,name=name+"_peptide_frames",lenfilter=0)
    if digest_peptides:
        inseq_translated_frames_peptides = []
        print("\tobtaining peptides ..."),
        for i in peptide_lengths:
            # I should indicate the coordinates for each frame in the future...
            a = [[''.join(i_tile) for i_tile in nmersub.tile_sequence(inseq=j,k=i,overlap=False)] for j in inseq_translated_frames] # this line is throwing an error
            # inseq_translated_frames_peptduleraides.append(a)
            inseq_translated_frames_peptides += a
            # need to allow the fill peptide if len(peptide) < tile
            print("\t %d %d-mers ..."%(len(a),i)),
        return(inseq_translated_frames_peptides,inseq_frames_pos,inseq_translated_frames)
    else:
        return(inseq_translated_frames,inseq_frames_pos,inseq_translated_frames)



# def nmersub.tile_sequence(inseq,name,peptide_lengths=[8,9,10,11],search_start_codons=True,search_stop_codons=True,search_alt_strand=False,digest_peptides=True):
#     # given the sequence,
#     return()

"""
# 3/26/2019
Next additions to streamline peptide outputs
1. Several of the frames are overlapping; combine them into a single consensus frame
2. Remove specific frames that are too small (below the peptide length limit)
"""

def main():
    # print out as a list

    # STEP 1: take in input
    parser = argparse.ArgumentParser(description='Produce direct translations of individual sequences')
    parser.add_argument('-f',
                        '--fasta',
                        help='FASTA file containing nucleotide sequence. REQUIRED',
                        type=str)
    parser.add_argument('-r',
                        '--revcomp',
                        help='Specify this flag if you want to also translate the reverse complement of your sequences.',
                        action='store_true')
    parser.add_argument('--clip_after_stop',
                        help='Specify this flag if you want to clip direct translations after a stop codon.',
                        action='store_true')
    parser.add_argument('-n',
                        '--name',
                        help='Job name.',
                        type=str,
                        default="job")
    parser.add_argument('--min',
                        help='Minimum peptide size',
                        type=int,
                        default=7)
    parser.add_argument('--overwrite',
                        help='Specify this flag if you want to overwrite output of the same name.',
                        action='store_true')
    parser.add_argument('--trypsinize',
                        help='Specify this flag if you want to simulate trypsin digestions. \
                        Also provided will be the density of trypsin frames per sequence, as well as a table indicating the probability \
                        that each possible trypsin-generated peptide is uniquely associated with each possible frame. \
                        This metric is essentially 1 - (p_{frame i of protein j maps to other proteins -j}), \
                        and the measurements are provided in a table of # of unique digestion products x # of total frames. \
                        These metrics can help with estimating the false positive rate of detecting a possible peptide \
                        generated from trypsinization of the theoretical protein/polypeptide. \
                        If any frame cannot be trypsinized, it will also be included as long as its size is >= MIN',
                        action='store_true')
    args = parser.parse_args()

    # STEP 1: read in FASTA file
    print("hello")
    infa = args.fasta # sys.argv[1] # infasta
    # filename = sys.argv[2]
    print("Now reading %s ..."%infa)
    inseq,inseqname = import_sequence(infa,concat_seq=False,assess_rc=args.revcomp)
    infilename = os.path.basename(infa).split(".txt")[0] # infa.split(".fa")[0]

    # remove file if it exists
    filename=args.name
    if os.path.exists(filename+".fa"):
        os.remove(filename+".fa")

    # STEP 2: define peptides
    print("\t ... defining peptides ... "),

    # if we have multiple, then we don't concat
    if type(inseq) is list:
        for i,j in zip(inseq,inseqname):
            j_rename = "_".join(j.split(":"))
            print(j_rename)
            # inpeptides = define_peptides(inseq=i,name=j_rename,digest_peptides=False,omit_termin=False)
            inpeptides = define_peptides(inseq=i,
            filename=filename,
            name=j_rename,
            digest_peptides=False,
            omit_termin=False,
            search_start_codons=False,
            search_stop_codons=False,
            clip_after_stop=args.clip_after_stop)
    else:
        inpeptides = define_peptides(inseq=inseq,filename=filename,name=infilename,digest_peptides=False,omit_termin=False)
    print("peptides obtained!")

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
else:
    print("Imported \'propose_nres_peptides\'")

# incorporte a function to just produce a nmersub.translation and not to



# in this script, we will
