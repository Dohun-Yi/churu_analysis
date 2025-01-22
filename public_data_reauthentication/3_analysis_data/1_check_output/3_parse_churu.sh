root=../../
workd=$root/analysis/1_check_output
datalist=$workd/SRA_churu_output_compare.1697001520.txt
dir_churu=$workd/1697001520/output/
path_out=$workd/1697001520.churu.out
while read srr etc; do
    path_churu=$dir_churu/$srr/churu.out
    if [ -f $path_churu ] ; then
        echo $srr
        head -n 5 $path_churu
        echo
    fi
done < $datalist > $path_out
