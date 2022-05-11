# mutagenesis_tools

The scripts in this repository can be used to simulate protein sequences from DNA, recode DNA with nonsynonymous or synonymous variants, and in general mutagenize a sequence *in silico*. The ideal applications are for:

1. Simulating SNV/indel and SV-fusion mutant sequences (genomic and transcriptomic), both nucleotide and amino acid.
3. Matching a peptide or protein sequence back to a mutant nucleotide of interest (for example, in matching protein mass spectrometry hits to a nucleotide database)
  a. match a protein or peptide sequence back to nucleotide sequences in a FASTA. If custom sequences or a set of FASTAs are provided, then this script will also determine which sequence(s) are the unique peptide templates

## Under construction
1. employ intron-retention and SJ-prediction tools to determine if a mutation hitting a gene is likely to cause a novel protein variant
2. Efficiently designing nucleotide libraries to cover the full space of protein variants for a given gene
  a. Generate SNVs required to efficiently produce all codon substitutions of a gene, i.e. what sets of SNVs and indels would be made such that you can cover all single-amino-acid substitutions, insertions, and deletions. This module would be useful for ordering sequences for protein library construction
3. Generate full-length protein isoforms from mutational data in a gene and fold them with AlphaFold2
4. Implement the full TESLA platform for neoantigen immunogenicity prediction
5. Implement the "mutant peptide affinity distribution" approach (that is, compare the distribution of MHC-peptide affinities of all peptides covering a site/sequence with the same distribution from a mutant sequence) to estimate the neoantigenic potential of a mutation


# Utilities

1. This repository implements `blast_tools` and can be paired with `sv_gene_analysis` to carry out fusion/chimeric protein sequence analysis and neoantigen identification (for SNVs, indels, and fusions)
2. `netMHCpan` is also implemented if you want to conduct MHC1-peptide affinity calculations
