dir1="/n/data1/hms/dbmi/park/vinay/pipelines/blast_tools"
dir2="/n/data1/hms/dbmi/park/vinay/pipelines/sv_gene_analysis/blast_tools/"

# rsync -a $dir1 $dir2

rsync -atuvP $dir1"/" $dir2 && rsync -atuvP $dir2"/" $dir1

# rsync -atuv /n/data1/hms/dbmi/park/vinay/pipelines/blast_tools/ /n/data1/hms/dbmi/park/vinay/pipelines/sv_gene_analysis/blast_tools
