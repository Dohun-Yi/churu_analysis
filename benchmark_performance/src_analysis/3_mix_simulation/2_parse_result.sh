root=../../
workd=$root/analysis/3_mix_simulation/

dna=churu_release_dna
rna=churu_release_rna

path_out=$workd/final_result_churu.txt

thp1() {
    f=$1
    echo "$(grep "THP1" $f | cut -f7)"
}
hct116() {
    f=$1
    echo "$(grep "HCT116" $f | cut -f7)"
}

echo -e "mix\tt1_d\thc_d\tt1_r\thc_r" > $path_out
for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    condition=${i1}_${i2}
    path_dna=$workd/$dna/$condition/churu.out
    path_rna=$workd/$rna/$condition/churu.out
    t1_d=$(thp1 $path_dna)
    hc_d=$(hct116 $path_dna)
    t1_r=$(thp1 $path_rna)
    hc_r=$(hct116 $path_rna)
    echo -e "$condition\t$t1_d\t$hc_d\t$t1_r\t$hc_r"
done >> $path_out
