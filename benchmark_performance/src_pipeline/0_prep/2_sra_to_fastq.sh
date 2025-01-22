root=../../
datalist=$root/datalist.txt
workd=$root/raw_data/
to=$workd/fastq
logs=./2_sra_to_fastq.logs
runfile=./2_sra_to_fastq.run.sh
cpu=30

mkdir -p $logs $to
ml sratoolkit/3.0.0
fastqdump=`which fasterq-dump`

while read sra sample lib layout; do
    sra_file=$workd/sra/$sra/$sra.sra
    log_file=$logs/$sra.log
    if [ ! -f $sra_file ]; then
        sra_file=${sra_file}lite
        if [ ! -f $sra_file ]; then
            echo "#warning! file does not exists $sra"
            continue
        fi
    fi

    if [ $layout == "PAIRED" ]; then
        cmd="$fastqdump $sra_file -e $cpu --split-3 -O $to/ > $log_file 2>&1"
    elif [ $layout == "SINGLE" ]; then
        cmd="$fastqdump $sra_file -e $cpu -O $to/ > $log_file 2>&1"
    else
        echo "error! not single or paired"
        exit
    fi
    echo $cmd
done < $datalist > $runfile

