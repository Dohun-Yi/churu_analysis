root=$PWD/../../
from=$root/processed_data/1_mapping_STAR_mm_noquant/
to=$root/processed_data/2_variant_call_varscan/
ref_fa=$root/reference/hg38/hg38.fa
datalist=$root/datalist.txt

ml java/jdk_13.0.1 # is a system default (20230815)
java=java
jar=$root/programs/varScan/v2.4.5/VarScan.v2.4.5.jar

runner=./1_mpileup_varscan_again_for_time.run.sh
while read sra sample lib layout; do
    mkdir $to/$sra/ -p
    bam_ori=$from/$sra/Aligned.sortedByCoord.out.bam
    bam=$to/$sra/Aligned.sortedByCoord.out.bam
    ln -srf $bam_ori $to/$sra/
    ln -srf $bam_ori.bai $to/$sra/
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
done < $datalist > $runner
