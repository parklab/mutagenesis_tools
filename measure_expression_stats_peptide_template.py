# this package will compute all possible template DNA sequences of a peptide and obtain the RNA-seq-level coverage (TPM is default) of the template sequence in one or more RNA-seq files (returns as a template).

## IMPORT STATEMENTS
import os,sys,re
import numpy,math,random
import itertools
import numpy as np
from collections import defaultdict
import argparse
from Bio import SeqIO
import pysam,pybedtools

sys.path.insert(0,'/n/data1/hms/dbmi/park/vinay/pipelines/sv_gene_analysis/')
from make_all_nmer_substitutions import reverse_seq,translation,codon_table,tile_sequence,complement_seq,reverse_codon_table,get_codon_set_with_gclim,get_codon_set
# from propose_nres_peptides import import_sequence,

## FILE HANDLING FUNCTIONS
def extract_crux_results(infile,colnames):
    try:
        data = pd.read_csv(infile, sep="\t",header=None)
        data.columns = colnames
        return(data)
    except:
        return(None)

def extract_crux_summary_stats (inpd,colnmaes,summarize_columns):
    # a = extract_crux_results(infile,colnames)
    a = inpd
    for i in summarize_columns:
        # print(i)
        i_mean = np.mean(a.loc[:,i])
        i_var = np.var(a.loc[:,i])
        i_z = (a.loc[:,i] - i_mean)/i_var
        a[i+'_zscore'] = i_z
    return(a)
        
def filter_crux_results (inpd,lower_bound_filters=None,upper_bound_filters=None):
    # print(lower_bound_filters)
    out_summary = inpd
    # apply the filters
    if lower_bound_filters != None:
        for i in lower_bound_filters:
            out_summary = out_summary.loc[out_summary[i] >= lower_bound_filters[i]]
    if upper_bound_filters != None:
        for i in upper_bound_filters:
            out_summary = out_summary.loc[out_summary[i] <= upper_bound_filters[i]]
    return(out_summary)


def get_files_dir (indir,filekey):
#    os.chdir(indir)
    for file in glob.glob(indir+"**/*"+filekey):
#        filestats = os.stat(file)
        if os.stat(file).st_size != 0:
            yield file

## PEPTIDE FILTERING
# generator object to store the peptide search results
def yield_protein_results(infiles,colnames,lower_filters,upper_filters):
    for i in infiles:
        # print(i)
        in_results = extract_crux_results(i,
                                         sample_match.columns)
        in_results['prop_ions_matched'] = in_results['b/y ions matched']/in_results['b/y ions total']
        with_summary = extract_crux_summary_stats(in_results,
                                                sample_match.columns,
                                                outstats)
        filtered_results = filter_crux_results(inpd=with_summary,
                                lower_bound_filters=lower_bound_filters,
                                upper_bound_filters=upper_bound_filters)
        if filtered_results is None:
            pass
        else:
            yield filtered_results
#             print(a.loc[:,['spectrum precursor m/z_zscore','peptide mass_zscore','sequence','e-value']])

# function to extract all unique peptides
def yield_peptides(peptide_results,peptide_colname):
    uniq_peptides = set()
    for i in peptide_results:
#         print(i)
        for j in i.loc[:,peptide_colname]:
            if j not in uniq_peptides:
                yield j
                uniq_peptides.add(j)

# consolidate all peptide filtering operations
def filter_peptide_matches(indir,colnames,filekey="fa_peptides_8aa.txt"):
    peptide_matches = get_files_dir(indir,filekey)
    # sample_match = pd.read_csv('/n/data1/hms/dbmi/park/vinay/analysis/tcga_neoA_analysis/sv_fusion_gene_analysis/double_valid_bkpt_analysis/kirc/2939_testcase/1kb_window_breakpoint/20190416_proteome_txome_comparisons/cptac3_ccrcc_search/2019_04_02_crux_comet_ensembl/17CPTAC_CCRCC_W_JHU_20180517_LUMOS_f16/comet.target.txt',header=0,sep="\t")
    sample_match_columns = colnames
    # get the unique peptides
    
    counter = 0
    lower_bound_filters = {'sp score':float(200),
                          'xcorr score':float(2)}
    upper_bound_filters = {'sp rank':float(2)}
    outstats = ['spectrum precursor m/z',
               'spectrum neutral mass',
               'peptide mass',
               'delta_cn',
               'sp score',
               'sp rank',
               'xcorr score',
               'prop_ions_matched',
               'total matches/spectrum']
    
    # obtain peptide results
    all_peptide_results = yield_protein_results(infiles=peptide_matches,
                                                colnames=sample_match_columns,
                                                outstats=outstats,
                                                lower_filters=lower_bound_filters,
                                                upper_filters=upper_bound_filters)
    
    # get unique peptides
    all_peptides = yield_peptides(all_peptide_results,'sequence')
    # get unique modified peptides
    all_mod_peptides = yield_peptides(all_peptide_results,'modified sequence')
    # return
    return(all_peptides,all_mod_peptides,all_peptide_results)

## MATCH PROTEIN TO GENE
def get_peptide_mapping_uniprot(fileobj,crux_results):
    for j in crux_results.iterrows():
        j_prot = j[1]['protein id'].split("|")[1]
        matches = []
        for line in fileobj:
            if re.search(j_prot+"\t",line):
                matches.append(line.strip("\n").split("\t")[-1])
        yield(j,matches)

def prepare_id_with_search(x):
    entry = x[0]
    line = x[1]
    found = search_peptide_entry(entry,line)
    print(entry,line,found)
    if found:
        tx = line.strip("\n").split("\t")
        return(tx[-1])
    else:
        tx = False
        
def prepare_id(x):
    entry = x[0]
    line = x[1]
    tx = line.strip("\n").split("\t")
    return(entry,tx[-1])

def search_peptide_entry(entry,line):
    # prepare the entry
    # search the line
    return(re.search(entry,line))

def get_uniprot_ids(crux_results):
    for j in crux_results.iterrows():
        j_prot = j[1]['protein id'].split("|")
        # print(j_prot)
        yield((j[0],j_prot[1]))

def match_file_proteins_to_ensembl(inresults,infile):
    mappings = get_uniprot_ids(inresults)
    for j in mappings:
        pdentry = j[1]
        pdind = j[0]
            # how to write a cleaner loop
        with open(infile) as datfile:
            result = map(prepare_id,((pdentry,line) for line in datfile if search_peptide_entry(pdentry,line)))
            yield (inresults.loc[pdind,:],';'.join([r[1] for r in result]))
#             yield ([r[1] for r in result])

def add_ensembl_matches(peptide_results,infile,counterlim):
    counter = 0
    for i in filtered_peptide_results:
        mappings = get_uniprot_ids(i)
        mappings, mappings_backup = itertools.tee(mappings)
        inmatches = []
        for j in mappings:
            pdentry = j[1]
            pdind = j[0]
            # how to write a cleaner loop
            with open(infile) as datfile:
                result = map(prepare_id,((pdentry,line) for line in datfile if search_peptide_entry(pdentry,line)))
                inmatches.append([r for r in result])
        yield ([m for m in mappings_backup],inmatches)
#         counter += 1
#         if counter > counterlim:
#             break

def add_ensembl_matches_faster(peptide_results,infile,counterlim):
    counter = 0
    for i in filtered_peptide_results:
        map_results = match_file_proteins_to_ensembl(i,infile)
        df1_test = []
        for line in map_results:
            newline = pd.Series.tolist(line[0])
            newline.append(line[1])
            df1_test.append(newline)
        newdf = pd.DataFrame(df1_test)
        newdf.columns = list(i.columns) + ['gene']
#         yield map_results
        yield newdf


# parse genefasta and yield the matches
def load_gene_fasta (infile):
    # Install the biopython library (if not already installed)

    # Import parts of Biopython
    # File path to your FASTA file
    path_to_file = infile   # <--- substitute by your local path

    # Open file with "with" statement to avoid problems with access 
    # to original file (in case computer hangs
    # or there will be any other problem)
    with open(path_to_file, mode='r') as handle:

        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        counter = 0
        for record in SeqIO.parse(handle, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            # description = record.description
            sequence = record.seq
            [gene,coordinates] = identifier.split("::")
            yield gene

def find_genes_in_fasta (infile,genelist):
    # Install the biopython library (if not already installed)

    # Import parts of Biopython
    # File path to your FASTA file
    path_to_file = infile   # <--- substitute by your local path

    # Open file with "with" statement to avoid problems with access 
    # to original file (in case computer hangs
    # or there will be any other problem)
    with open(path_to_file, mode='r') as handle:

        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        counter = 0
        for record in SeqIO.parse(handle, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            [gene,coordinates] = identifier.split("::")
            if gene in genelist:
                yield (gene,coordinates)
    
##### TO DO
## TRANSCRIPT IDENTIFICATION SCRIPTS
def get_all_template_dna (inseq,seqref):
    # for a given peptide sequence, what are all possible DNA sequences that could have translated to the peptide?
    # use the get_codon_set function
    return()

def generate_dna_templates(peptides,seqref):
    for i in peptides:
        yield get_codon_set_traversal(seq=i,
                                      seqref = seqref)

def get_template_dna_at_reference (inseq,inref,ingenome):
    # for a given peptide sequence, what are all possible DNA sequences AT THE MATCHING REFERENCE SEQUENCE that could have translated to the peptide?
    # first, find the 
    codon_set = get_codon_set(seq=inseq)
    return(codon_set)
    return()

def get_template_dna_by_tblastn ():
    # this script runs an external function 'tblastn' from the BLAST package to
    # incorporate a function here to run the BLAST shell script and store the output
    return()
    

### UNTESTED

def get_reference_sequence (refname,inref,ingenome):
    # obtains the reference sequence for a particular reference gene "refname" so that the relevant DNA sequence can be extracted
    return()

def get_tx_coverage_seq (inbam,intervals):
    # using the pybedtools coverage command, obtain the RNA-seq read-level coverage at a peptide
    # step 1: import the intervals from 'intervals'
    # step 2: read in the BAM file from 'inbam'
    # step 3: perform the intersection
    return()

### MAIN: Decant all lines to carry out functions here. 

def main(indir):
    
    # get peptide matches 
    # indir = "./cptac3_ccrcc_search/2019_04_16_test_batch_match_crux_comet_run2/"
    colnames = ['scan', 'charge', 'spectrum precursor m/z', 'spectrum neutral mass',
                'peptide mass', 'delta_cn', 'sp score', 'sp rank', 'xcorr score',
                'xcorr rank', 'b/y ions matched', 'b/y ions total',
                'total matches/spectrum', 'sequence', 'modified sequence',
                'modifications', 'protein id', 'flanking aa', 'e-value']
    outstats = ['spectrum precursor m/z',
               'spectrum neutral mass',
               'peptide mass',
               'delta_cn',
               'sp score',
               'sp rank',
               'xcorr score',
               'prop_ions_matched',
               'total matches/spectrum']
    
    # obtain the peptide matches
    peptides,mod_peptides,filtered_peptide_results = filter_peptide_matches(indir,colnames,filekey="fa_peptides_8aa.txt")
    # match the UniprotID of the observed peptides' source proteins to the gene name of interest
    infile = "/n/data1/hms/dbmi/park/vinay/referenceAnnotation/hg19/HUMAN_9606_idmapping_ensembl_gene.dat"
    a = add_ensembl_matches_faster(filtered_peptide_results,infile,2)
    a,top_peptides_ensembl_matched = itertools.tee(a)
    # now, to mine the sequences from the ensembl FASTA file...
    uniq_genes = yield_peptides(top_peptides_ensembl_matched,'gene')
    # 1. load the FASTA file with the ENSEMBL gene names
    genefasta = "/n/data1/hms/dbmi/park/vinay/referenceAnnotation/hg19/Homo_sapiens.GRCh37.75.gene.fa"
    refgenes = load_gene_fasta(genefasta)
    # 2. Search for genes...

    

    return()


### END FUNCTIONS

if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   main()
else:
    print("Imported \'measure_expression_stats_peptide_template\'")
