root=../../
workd=$root/analysis/1_benchmark/
datalist=$root/datalist.txt

from_star=$root/processed_data/1_mapping_STAR_mm_noquant/
from_gatk=$root/processed_data/2_variant_call_gatk_rna_seq/
from_varscan=$root/processed_data/2_variant_call_varscan/
from_churu=$root/processed_data/3_identification_CHURU/

output=$workd/time.out.txt

get_time_star() {
    dir_sample=$1
    log=$dir_sample/Log.final.out
    tic=`date -d "$(grep 'Started job on' $log | cut -f2)" +%s`
    toc=`date -d "$(grep 'Finished on' $log | cut -f2)" +%s`
    echo $((toc-tic))
}
get_time_gatk() {
    dir_sample=$1
    f1=$(ls $dir_sample/*exe*/*/*/*gtfToCalling*/ex*/script)
    f2=$(ls $dir_sample/*exe*/*/*/*VariantFilt*/ex*/*tbi)
    tic=`stat -c%Z $f1`
    toc=`stat -c%Z $f2`
    echo $((toc-tic))
}
get_time_churu() {
    dir_sample=$1
    f1=$dir_sample/churu.log
    if [ ! -f $f1 ]; then
        echo '-1'
    else
        line=$(tail -n 1 $f1)
        if [ $(echo $line | grep -c 'run time') -eq 1 ]; then
            echo "$(echo $line | cut -d' ' -f5)"
        else
            echo "-1"
        fi
    fi
}
get_time_varscan() {
    dir_sample=$1
    log=$dir_sample/varscan_time.log
    if [ -f $log ] && [ $(wc -l $log | cut -d' ' -f1) -eq 2 ]; then
        tic=`date -d "$(head -n 1 $log)" +%s`
        toc=`date -d "$(tail -n 1 $log)" +%s`
        echo $((toc-tic))
    else
        echo -1
    fi
}
echo -e "SRA\tSTAR\tGATK\tVARSCAN\tCHURU" > $output
while read sra cell lib layout; do
    line="$sra\t"
    line=$line"$(get_time_star $from_star/$sra)\t"
    line=$line"$(get_time_gatk $from_gatk/$sra)\t"
    line=$line"$(get_time_varscan $from_varscan/$sra)\t"
    line=$line"$(get_time_churu $from_churu/$sra)\t"
    echo -e $line
done < $datalist >> $output
