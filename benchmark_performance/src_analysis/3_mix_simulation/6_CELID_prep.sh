root=../../
workd=$root/analysis/3_mix_simulation/
mode=celid
from=$workd/$mode/
to=$from

runner=./6_CELID_prep.run.sh

dir_celid=$root/programs/CeL-ID/
parser=$dir_celid/vcfFilter_DP_Freq.pl

for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    sra=${i1}_${i2}
    mkdir $to/$sra/ -p

    vcf=$to/$sra/output.pileup.vcf
    lst=$to/$sra/samplelist.txt
    echo 'Sample1' > $lst
    cmd1="sed 's/^chr//g' $vcf -i"
    cmd2="perl $parser --input $vcf --list $lst --DP=30 --AD=5 --FREQ=25"
    
    echo "$cmd1; $cmd2"
done > $runner
