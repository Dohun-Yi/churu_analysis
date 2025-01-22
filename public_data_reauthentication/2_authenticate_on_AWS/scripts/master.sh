srafile=/data1/instance_resource/SRA_Accessions.batch.txt
source $(dirname $0)/identity.sh
bucket=s3://dohun/sra/$BATCH/output
cpu=$(nproc)

logfile_master=/data1/log/master.log
logfile_down=/data1/log/SRA_downloaded.log
logfile_decomp=/data1/log/SRA_decomped.log
logfile_run=/data1/log/SRA_finished.log
touch $logfile_master $logfile_down $logfile_decomp $logfile_run

churu=/data1/churu/churu.sh
fasterq=/usr/local/ncbi/sra-tools/bin/fasterq-dump

# use six disks. churu is limited to 2 process due to memory.
maxjob=6
pidInUse=('empty' 'empty' 'empty' 'empty' 'empty' 'empty')
disks=('/data1/' '/data2/' '/data3/' '/data4/' '/data5/' '/data6/')


run_demo() {
    cd /data1/churu/
    rm -rf demo/output
    cmd="./churu.sh -t $cpu -1 demo/demo_1.fastq -2 demo/demo_2.fastq -o demo/output -g LoadAndKeep"
    $cmd 1> demo/log.out 2> demo/log.err
    if [ $? != 0 ]; then
        message="{\"text\":\"$INSTANCENAME failed demo run. kill me\"}"
        # Send kill message to my Slack and exit
        # curl -X POST -H 'Content-type: application/json' --data "$message" {URL}
        bash $0/suicide.sh
        exit 1
    fi
    cd -
}
run_demo

download_and_process() {
    acccession=$1
    disk=$2
    dir_downloading=$disk/downloading/
    dir_downloaded=$disk/downloaded/
    mkdir -p $dir_downloading $dir_downloaded
    mkdir -p $dir_downloading/$accession/

    cmd="aws s3 cp s3://sra-pub-run-odp/sra/$accession/$accession $dir_downloading/$accession/$accession --no-sign-request --no-progress"
    tic=$(date +%s)
    echo $cmd; $cmd
    lastexit=$?
    toc=$(date +%s)
    tictoc=$((toc-tic))s
    echo "[$(date)] download($lastexit) $tictoc" >> $dir_downloading/$accession/time
    if [ $lastexit -eq 0 ]; then
        mv $dir_downloading/$accession $dir_downloaded/$accession
        echo "[$(date)] $accession successfully downloaded, $tictoc" >> $logfile_down
        echo >&2 "[$(date)] $accession successfully downloaded, $tictoc"
    else
        rm -rf $dir_downloading/$accession
        echo "[$(date)] $accession failed to download, $tictoc - exit code $lastexit" >> $logfile_down
        echo >&2 "[$(date)] $accession failed to download, $tictoc - exit code $lastexit"
    fi
    return $lastexit
}

decomp() {
    acccession=$1
    disk=$2
    dir_downloaded=$disk/downloaded/
    dir_decomped=$disk/decomped/
    dir_tmp=$disk/tmp/
    mkdir -p $dir_downloaded $dir_decomped

    tic=$(date +%s)
    if [ $(grep -c "$accession failed" $logfile_down) -gt 0 ]; then
        echo >&2 "[$(date)] $accession failed - download failed"
        echo "[$(date)] $accession failed - download failed" >> $logfile
        return 1
    elif [ $(grep -c "$accession success" $logfile_down) -gt 0 ]; then
        tic=$(date +%s)
        while [ ! -f $dir_downloaded/$accession/$accession ]; do
            toc=$(date +%s)
            if [ $((toc-tic)) -gt 60 ]; then
                echo >&2 "[$(date)] $accession failed - transfer wait timeout"
                echo "[$(date)] $accession failed - transfer wait timeout" >> $logfile_decomp
                return 1
            fi
            echo >&2 "[$(date)] waiting for $accession... (transfer) $((toc-tic))s"
            sleep 5
        done
    fi
    ###
    cmd="$fasterq $dir_downloaded/$accession/$accession -O $dir_downloaded/$accession/ -t $dir_tmp"
    tic=$(date +%s)
    echo $cmd; $cmd
    lastexit=$?
    toc=$(date +%s)
    tictoc=$((toc-tic))s
    echo "[$(date)] decomp($lastexit) $tictoc" >> $dir_downloaded/$accession/time
    if [ $lastexit -eq 0 ]; then
        rm $dir_downloaded/$accession/$accession -f
        mv $dir_downloaded/$accession $dir_decomped/
        echo "[$(date)] $accession successfully decomped, $tictoc" >> $logfile_decomp
        echo >&2 "[$(date)] $accession successfully decomped, $tictoc"
    else
        rm -rf $dir_downloaded/$accession
        echo "[$(date)] $accession failed to decomp, $tictoc - exit code $lastexit" >> $logfile_decomp
        echo >&2 "[$(date)] $accession failed to decomp, $tictoc - exit code $lastexit"
    fi
    return $lastexit
}

identify() {
    accession=$1
    disk=$2
    libsource=$3
    layout=$4
    dir_decomped=$disk/decomped/
    dir_processing=$disk/processing/
    outdir=$dir_processing/$accession/churu/
    logfile_churu=$dir_processing/$accession/$accession.churu.log
    errfile_churu=$dir_processing/$accession/$accession.churu.err
    mkdir -p $dir_decomped $dir_processing

    # Check if file exists
    if [ ! -d $dir_decomped/$accession ]; then
        echo "[$(date)] $accession failed - file does not exists.. skip" >> $logfile_run
        echo >&2 "[$(date)] $accession failed - file does not exists.. skip"
        return 1
    fi
    # mv the directory
    mv $dir_decomped/$accession $dir_processing/

    # Check layout
    if [ $layout = "PAIRED" ]; then
        fq1=$dir_processing/$accession/${accession}_1.fastq
        fq2=$dir_processing/$accession/${accession}_2.fastq
        input="-1 $fq1 -2 $fq2"
    elif [ $layout = "SINGLE" ]; then
        fq1=$dir_processing/$accession/${accession}.fastq
        input="-1 $fq1"
    fi

    # Check library source
    if [ $libsource = "GENOMIC" ]; then
        nucl=dna
    elif [ $libsource = "TRANSCRIPTOMIC" ]; then
        nucl=rna
    else
        echo "[$(date)] $accession failed - has unknown libsource. skip." >> $logfile_run
        echo >&2 "[$(date)] $accession failed - has unknown libsource. skip."
        return
    fi
    # Run CHURU
    cmd="$churu -t $cpu $input -o $outdir -g LoadAndKeep -n $nucl"
    echo $cmd
    tic=$(date +%s)
    $cmd 1> $logfile_churu 2> $errfile_churu
    lastexit=$?
    toc=$(date +%s)
    tictoc=$((toc-tic))s
    echo "[$(date)] identified($lastexit) $tictoc" >> $dir_processing/$accession/time

    # Prep files for upload
    find $dir_processing/$accession -exec du -h {} \; 2>&1 > $dir_processing/$accession/sizes
    head -n 100 $dir_processing/$accession/*q > $dir_processing/$accession/head_fq

    # Upload logs
    aws s3 cp $logfile_churu $bucket/$accession/ --no-progress
    aws s3 cp $errfile_churu $bucket/$accession/ --no-progress
    aws s3 cp $dir_processing/$accession/time $bucket/$accession/ --no-progress
    aws s3 cp $dir_processing/$accession/sizes $bucket/$accession/ --no-progress
    aws s3 cp $dir_processing/$accession/head_fq $bucket/$accession/ --no-progress
    aws s3 cp $outdir/churu.out $bucket/$accession/ --no-progress
    aws s3 cp $outdir/output.pileup $bucket/$accession/ --no-progress
    aws s3 cp $outdir/output.pileup.vcf $bucket/$accession/ --no-progress
    # Leave log
    if [ $lastexit -eq 0 ]; then
        echo "[$(date)] $accession successfully identified, $tictoc" >> $logfile_run
        echo >&2 "[$(date)] $accession successfully identified, $tictoc"
    else
        echo "[$(date)] $accession failed, $tictoc - exit code $lastexit" >> $logfile_run
        echo >&2 "[$(date)] $accession failed, $tictoc - exit code $lastexit"
    fi
    # remove files
    rm -rf $dir_processing/$accession
    return $lastexit
}

run() {
    accession=$1
    disk=$2
    libsource=$3
    layout=$4
    tic=$(date +%s)
    # Clear disk
    rm -rf $disk/downloading/ $disk/downloaded/ $disk/decomped/ $disk/tmp/ $disk/decomped/ $disk/processing/

    echo "[$(date)] $accession started on $disk..." | tee -a $logfile_master

    # download
    download_and_process $accession $disk 1>>/data1/log/down.out 2>>/data1/log/down.err
    lastexit=$?
    if [ ! $lastexit -eq 0 ]; then
        echo "[$(date)] $accession failed downloding on $disk." | tee -a $logfile_master
        echo "[$(date)] $accession failed downloding on $disk." >> $logfile_decomp
        echo "[$(date)] $accession failed downloding on $disk." >> $logfile_run
        exit $lastexit
    fi

    # decompress
    decomp $accession $disk  1>>/data1/log/decomp.out 2>>/data1/log/decomp.err
    lastexit=$?
    if [ ! $lastexit -eq 0 ]; then
        echo "[$(date)] $accession failed decompressing on $disk." | tee -a $logfile_master
        echo "[$(date)] $accession failed decompressing on $disk." >> $logfile_run
        exit $lastexit
    fi

    # identify
    while [ $(pgrep -f 'churu.sh' | wc -l) -ge 2 ]; do
	sleep 5
    done
    identify $accession $disk $libsource $layout 1>>/data1/log/identify.out 2>>/data1/log/identify.err
    lastexit=$?
    if [ ! $lastexit -eq 0 ]; then
        echo "[$(date)] $accession failed identifying on $disk." | tee -a $logfile_master
        exit $lastexit
    fi
    toc=$(date +%s)
    echo "[$(date)] $accession successfully ended on $disk... took $((toc-tic))s" | tee -a $logfile_master
    return 0
}

while read accession libsource layout; do
    # check empty disk
    hasEmpty="false"
    while [ $hasEmpty = "false" ]; do
        for i in $(seq 0 $((maxjob-1))); do
            if [ "${pidInUse[$i]}" == "empty" ] || \
                $(! kill -0 ${pidInUse[$i]} 2>/dev/null) ; then 
                hasEmpty="true"
            fi
        done
        sleep 1
    done
    for i in $(seq 0 $((maxjob-1))); do
        # skip if disk is occupied
        if [ "${pidInUse[$i]}" != "empty" ] && \
            $(kill -0 ${pidInUse[$i]} 2>/dev/null) ; then 
            continue
        else
	# if available, start job and break loop
            run $accession ${disks[$i]} $libsource $layout &
            pidInUse[$i]=$!
            break
        fi
    done
    sleep 1
done < <(cut -f1,22,24 $srafile)
wait
