log=./tracker.log
last_n=0
total=203916
while true; do
    d=$(date)
    n=$(aws s3 ls s3://dohun/sra/1697001520/output/ 2>/dev/null | wc -l)
    percent=$(printf %.3f%% "$((10**9 * $n/$total))e-7")
    echo -e "$d\t$(($n-$last_n))\t$n\t$percent" | tee -a $log
    last_n=$n
    sleep 300
done
