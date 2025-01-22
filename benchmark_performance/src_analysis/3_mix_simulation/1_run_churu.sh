root=../../
workd=$root/analysis/3_mix_simulation/
from_fq=$workd/fastq_mix/


cpu=64
churu=$root/program/churu/churu
genome=$root/program/churu/reference/hg38.fa
var=$root/program/churu/reference/variants
idx=$root/program/churu/reference/idx
param_common="-g $genome -r $idx -s $var"

# RNA mode (with expression)
mode=release_rna
param="-n rna"
runfile=./1_run_churu_rna.run.sh
to=$workd/churu_$mode/

mkdir -p $to

if [ -f $runfile ]; then
    rm $runfile
fi

for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    fq1=$from_fq/${i1}_${i2}/mix_1.fastq
    fq2=$from_fq/${i1}_${i2}/mix_2.fastq
    outpath=$to/${i1}_${i2}
    log=$outpath.log
    cmd="$churu -f -t $cpu -1 $fq1 -2 $fq2 -o $outpath/ $param $param_common >> $log 2>&1"
    echo $cmd
done > $runfile


# DNA mode (without expression)
mode=release_dna
param="-n dna"
runfile=./1_run_churu_dna.run.sh
to=$workd/churu_$mode/

mkdir -p $to

if [ -f $runfile ]; then
    rm $runfile
fi

for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    fq1=$from_fq/${i1}_${i2}/mix_1.fastq
    fq2=$from_fq/${i1}_${i2}/mix_2.fastq
    outpath=$to/${i1}_${i2}
    log=$outpath.log
    cmd="$churu -f -t $cpu -1 $fq1 -2 $fq2 -o $outpath/ $param $param_common >> $log 2>&1"
    echo $cmd
done > $runfile
