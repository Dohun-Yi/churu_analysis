root=$PWD/../../
workd=$root/analysis/3_mix_simulation/
mode=celid
from=$workd/$mode/
to=$from
ref_fa=$root/reference/hg38/hg38.fa

ml java/jdk_13.0.1

java=java
jar=$root/programs/varScan/v2.4.5/VarScan.v2.4.5.jar

# config
runner=5_mpileup_varscan.run.sh
for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    sra=${i1}_${i2}
    mkdir $to/$sra/ -p

    # pileup=$to/$sra/output.pileup
    bam=$to/$sra/Aligned.sortedByCoord.out.bam
    vcf=$to/$sra/output.pileup.vcf
    log=$to/$sra/varscan_time.log
    cmd_log="date >> $log"

    cmd1="samtools mpileup -B -f $ref_fa $bam"
    cmd2="java -jar $jar mpileup2snp --output-vcf 1"
    cmd="$cmd1 | $cmd2 > $vcf"
    cmd_full="$cmd_log; $cmd1 | $cmd2 > $vcf; $cmd_log"

    if [ -f $log ] && [ $(wc -l $log | cut -d" " -f1) -eq 2 ]; then
        continue
    else
        rm -f $vcf $log
        echo $cmd_full
    fi
done > $runner
