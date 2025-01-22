root=../../
workd=$root/analysis/3_mix_simulation/
mode=celid
workd=$workd/$mode/

runner=./7_CELID.run.sh

dir_celid=$root/programs/CeL-ID/
identifier=$dir_celid/cellAuthentication_args.R

for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    sra=${i1}_${i2}
    mkdir $workd/$sra/ -p

    DP=$workd/$sra/output.pileup_DP.txt
    FREQ=$workd/$sra/output.pileup_FREQ.txt
    lst=$workd/$sra/samplelist.txt
    identity=$workd/$sra/identity.out
    cmd="Rscript $identifier $DP $FREQ > $identity"
    echo $cmd

done > $runner
