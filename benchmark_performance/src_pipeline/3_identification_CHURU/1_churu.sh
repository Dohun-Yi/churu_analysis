root=../../
workd=$root/processed_data/3_identification_CHURU/
from_fq=$root/raw_data/fastq/
datalist=$root/datalist.txt

cpu=32
churu=$root/program/churu
genome=$root/program/churu/reference/hg38.fa
var=$root/program/churu/reference/variants
idx=$root/program/churu/reference/idx

runner=1_churu.run.sh
rm -f $runner
while read sra sample lib layout; do
    mkdir $workd/$sra/ -p
    outpath=$workd/$sra/churu
    log=$workd/$sra/churu.log

    fq1=$from_fq/${sra}_1.fastq
    fq2=$from_fq/${sra}_2.fastq
    fqs=$from_fq/${sra}.fastq
    param="-M LoadAndKeep -t $cpu -g $genome -r $idx -s $var"
    if [ -f $fq1 ] && [ -f $fq2 ] ; then
        cmd="$churu -1 $fq1 -2 $fq2 -o $outpath $param >> $log 2>&1"
    elif [ -f $fqs ] ; then
        cmd="$churu -1 $fqs -o $outpath $param >> $log 2>&1"
    else
        echo "warning, fastq file does not exists for $sra"
        continue
        # exit
    fi
    echo $cmd >> $runner
done < $datalist
