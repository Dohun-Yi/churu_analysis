root=../../
workd=$root/analysis/1_check_output/
path_in=$workd/1697001520.parse.table
path_filt=$workd/1697001520.parse.filt.table
path_fail=$workd/1697001520.parse.fail.table

awk -F"\t" '$15!="-"' $path_in > $path_filt
awk -F"\t" '$15=="-"' $path_in > $path_fail
