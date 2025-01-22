root=../../
from=$root/processed_data/2_variant_call_varscan/
to=$root/processed_data/3_identification_CELID/
datalist=$root/datalist.txt

dir_celid=$root/programs/CeL-ID/
parser=$dir_celid/vcfFilter_DP_Freq.pl

runner=./0_prep.run.sh
while read sra sample lib layout; do
    mkdir $to/$sra/ -p
    vcf_ori=$from/$sra/output.pileup.vcf
    vcf=$to/$sra/`basename $vcf_ori`
    lst=$to/$sra/samplelist.txt
    echo 'Sample1' > $lst
    rm -f $vcf
    cp $vcf_ori $vcf
    cmd1="sed 's/^chr//g' $vcf -i"
    cmd2="perl $parser --input $vcf --list $lst --DP=30 --AD=5 --FREQ=25"
    
    echo "$cmd1; $cmd2"
done < $datalist > $runner
