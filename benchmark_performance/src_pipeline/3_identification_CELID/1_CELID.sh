root=../../
workd=$root/processed_data/3_identification_CELID/
datalist=$root/datalist.txt

dir_celid=$root/programs/CeL-ID/
identifier=$dir_celid/cellAuthentication_args.R

runner=./1_CELID.run.sh
while read sra sample lib layout; do
    mkdir $workd/$sra/ -p
    DP=$workd/$sra/output.pileup_DP.txt
    FREQ=$workd/$sra/output.pileup_FREQ.txt
    lst=$workd/$sra/samplelist.txt
    identity=$workd/$sra/identity.out
    cmd="Rscript $identifier $DP $FREQ > $identity"
    echo $cmd

done < $datalist > $runner
