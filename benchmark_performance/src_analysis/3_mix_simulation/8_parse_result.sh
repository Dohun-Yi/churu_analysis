root=../../
workd=$root/analysis/3_mix_simulation/

path_out=$workd/final_result_celid.txt

thp1() {
    f=$1
    first=`awk 'FNR==2' $f | sed 's/^[ \t]*//;s/[ \t]*$//'`
    mix=`awk 'FNR==5' $f | sed 's/^[ \t]*//;s/[ \t]*$//'`
    percent1=`awk 'FNR==5' $f | cut -d"=" -f2 | cut -d"%" -f1`
    percent1=$(python -c "print ($percent1)" )
    percent2=$(python -c "print 100 - $percent1" )

    if [[ "$mix" = *"HCT_116"* ]] && [[ "$first" = *"THP.1"* ]]; then
        echo "$percent2%\t$percent1%"
    elif [[ "$mix" = *"THP.1"* ]] && [[ "$first" = *"HCT_116"* ]]; then
        echo "$percent1%\t$percent2%"
    elif [[ "$first" = *"THP.1"* ]]; then
        echo "$percent2%\t-"
    elif [[ "$first" = *"HCT_116"* ]]; then
        echo "-\t$percent2%"
    else
        echo "-\t-"
    fi
}

echo -e "mix\tt1_c\thc_c" > $path_out
for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    condition=${i1}_${i2}
    # echo $condition
    path=$workd/celid/$condition/identity.out
    # echo $path
    t1="$(thp1 $path)"
    echo -e "$condition\t$t1"
done >> $path_out
