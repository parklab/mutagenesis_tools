# BLAST tools

## Introduction
The scripts in this directory implement BLAST in order to map peptides or sequences back to some reference of interest. Right now, run-time parameters must be specified as a config file. An example is given in 

## Instructions

### For batch-core

The batch-core scripts parallelize the searches for multiple samples as job arrays. The relevant scripts are `./batch_tblastn_crux_comet_peptides_fromlist.sh` and `./core_tblastn_crux_comet_peptides_fromlist.sh`. 

1. Specify the paths of your input files in one file called `infiles.txt` (you can vary this name out). A sample, `test_infiles_first10.txt` is given here.

2. Set your config file. `sample_config_tblastn.txt` is an example of a config file containing all necessary input parameters. 

```
blastdb <A BLAST index that BLAST will use to search the directory. 
        If you would like to build your own BLAST directory, please see the 
        script `core_makeblastdb.sh`>

evalue  <What e-value should you use as a cutoff? The e-value is 
        meant to indicate the significance of the sequence pairings observed 
        in the alignment. The default cutoff is 1e-5. See 
        https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect and https://en.wikipedia.org/wiki/BLAST_(biotechnology)#Algorithm for further reading>

# the `splitlines` and `nsamplequeue` options are necessary for submitting 
your jobs to SLURM as job arrays. 

splitlines  <if you are submitting `infiles.txt`, how many files do you want 
            to run per job array?>

nsamplequeue  <when submitting the job array, how many jobs in each array should
              be allowed to run at the same time>
```

3. Submit your job to SLURM as follows:

```
bash ./batch_tblastn_crux_comet_peptides_fromlist.sh test_infiles_first10.txt sample_config_tblastn.txt
# Submitted job #######
```

For each list you submit, a separate output directory should be created (the directory should carry the same name as your input file. If your input file paths all have the same name but are each located in a different directory, please make sure that the name itself is different; otherwise, you results will be overwritten). 

You can browse the output of a sample run (using the above instructions) here: `/n/data1/hms/dbmi/park/vinay/analysis/tcga_neoA_analysis/sv_fusion_gene_analysis/double_valid_bkpt_analysis/all_samples/pipeline_testruns/20190719_run5`

### If you would like to build a custom BLAST reference

Make sure you have your FASTA file of interest. **WARNING: BLAST will not accept compressed fasta files for building references!**. Then, 

```
sbatch ./core_makeblastdb.sh \
<reference .fa file> \
<the prefix you would like to give the entries> \
<output directory>
```
The BLASTDB directory should be created in your output directory. **NOTE** that the reference is specified to the actual alignment script as `<output directory>/<prefix used to name the file>`. See the `blastdb` option in the config files for reference. 

Two references involving `hg19` are available for use:
* ENSEMBL hg19 whole genome at `/n/data1/hms/dbmi/park/vinay/referenceGenomes/blastdb/ensembl_reference/`. The config file that submits this reference is `test_config.txt`. 
* ENSEMBL hg19 CDS at `/n/data1/hms/dbmi/park/vinay/referenceGenomes/blastdb/ensembl_cds/`. The config file that submits this reference is `test_config_with_cds.txt`. 

Please see `/n/data1/hms/dbmi/park/vinay/analysis/tcga_neoA_analysis/sv_fusion_gene_analysis/double_valid_bkpt_analysis/all_samples/pipeline_testruns/20190719_run6` for an example on how the hg19 CDS run appears. 

