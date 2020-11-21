# mutagenesis_tools

The scripts in this repository can be used to simulate protein sequences from DNA, recode DNA with nonsynonymous or synonymous variants, and in general mutagenize a sequence *in silico*. The ideal applications are for:

1. Simulating SNV/indel and SV-fusion mutant sequences, both nucleotide and amino acid.
2. Efficiently designing nucleotide libraries to cover the full space of protein variants for a given gene
3. Matching a peptide or protein sequence back to a mutant nucleotide of interest (for example, in matching protein mass spectrometry hits to a nucleotide database)


# Utilities

1. This repository implements `blast_tools` and can be paired with `sv_gene_analysis` to carry out fusion/chimeric protein sequence analysis and neoantigen identification (for SNVs, indels, and fusions)
2. `netMHCpan` is also implemented if you want to conduct MHC1-peptide affinity calculations

