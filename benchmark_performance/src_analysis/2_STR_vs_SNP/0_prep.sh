root=../../
workd=$root/analysis/2_STR_vs_SNP/
fetch_from=$root/analysis/1_benchmark/identity
cello=$root/reference/cellosaurus/cellosaurus.txt
datalist=$root/datalist.txt

mkdir -p $workd/

ln -srf $fetch_from $workd/
ln -srf $cello $workd/
ln -srf $datalist $workd/
