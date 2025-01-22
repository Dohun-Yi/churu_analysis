resources=/data1/instance_resource/
srafile=$resources/SRA_Accessions.batch.txt
dir_log=/data1/log/
logfile_down=$dir_log/SRA_downloaded.log
logfile_decomp=$dir_log/SRA_decomped.log
logfile_run=$dir_log/SRA_finished.log
logfile_sar=$dir_log/sar.out
source $(dirname $0)/identity.sh
bucket=s3://dohun/sra/$BATCH/logs/$INDEX

mkdir -p $dir_log
touch $logfile_down $logfile_decomp $logfile_run $dir_log/progress.log $dir_log/progress.hourly.log

check_file() {
    curtime_s=$1 # in sec
    file=$2
    n=0; s=0; f=0
    while read line; do
	thetime=$(echo $line | cut -f1 -d"|")
	logline=$(echo $line | cut -f2 -d"|")
    thetime_s=$(date -d "$thetime" +%s)
	diff=$((curtime_s - thetime_s))
	#echo >&2 $curtime_s, $thetime_s, $diff, $thetime, $logline
        if [ $diff -lt 3600 ]; then
            n=$((n+1))
            if [[ "${logline,,}" =~ "fail" ]] ; then
                f=$((f+1))
            fi
        fi
    done < <(sed 's/\[//g' $file | sed 's/\]/|/g')
    echo "$n (fail=$f)"
}

getline(){
    file=$1
    echo "$(wc -l $file | cut -f1 -d' ')"
}

format()
{
  MESSAGE="$1"
  cat <<EOF
{
  "text":"$MESSAGE"
}
EOF
}

post() {
  MESSAGE=$(format "$1")
  # send progress information to my Slack
  # curl -X POST -H 'Content-type: application/json' --data "$MESSAGE" {URL}
}


curtime=$(date)
curtime_s=$(date -d "$curtime" +%s)
n_down=$(check_file $curtime_s $logfile_down)
n_decomp=$(check_file $curtime_s $logfile_decomp)
n_run=$(check_file $curtime_s $logfile_run)
n_down_t=$(check_file 0 $logfile_down)
n_decomp_t=$(check_file 0 $logfile_decomp)
n_run_t=$(check_file 0 $logfile_run)
msg="[$curtime] this hour: Download $n_down, Decomp $n_decomp, Run $n_run."
msg="$msg  in total: Download $n_down_t, Decomp $n_decomp_t, Run $n_run_t"
post "$INSTANCENAME: $msg" 2>/dev/null
echo $msg >> $dir_log/progress.hourly.log
echo "DOWN:    $(getline $logfile_down) / $(getline $srafile)" > $dir_log/progress.log
echo "DECOMP:  $(getline $logfile_decomp) / $(getline $srafile)" >>$dir_log/progress.log
echo "RUN:     $(getline $logfile_run) / $(getline $srafile)" >> $dir_log/progress.log
aws s3 cp $dir_log $bucket --recursive
