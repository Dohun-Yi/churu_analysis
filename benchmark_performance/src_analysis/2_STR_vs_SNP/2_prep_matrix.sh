root=../../
workd=$root/analysis/2_STR_vs_SNP/
fetch_from=$root/analysis/1_benchmark/identity
cello=$root/reference/cellosaurus/cellosaurus.txt
datalist=$root/datalist.txt
matrix=$workd/matrix_lineage.txt

gzip -c $workd/`basename $matrix` > $workd/`basename $matrix`.gz
awk -F":" 'FNR!=1 {print $1}' $matrix > $workd/cell_list.txt
