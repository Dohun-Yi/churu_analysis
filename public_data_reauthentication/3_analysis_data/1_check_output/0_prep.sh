root=../../
workd=$root/analysis/1_check_output
bucket=s3://dohun/sra/1697001520
mkdir -p $workd

ln -srf $root/resources/SRA_Accessions_with_expr.uniq.tab $workd/
aws s3 cp $bucket $workd/1697001520/ --recursive
