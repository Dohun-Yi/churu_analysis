root=../../
workd=$root/analysis/3_mix_simulation/
from_fq=$workd/fastq_mix/
mode=celid
to=$workd/$mode/
runfile=3_mapping_STAR_mm_noquant.run.sh

ml STAR/2.5.4

# config
cpu=80
ref_star=$root/reference/hg38/refseq/GRCh38_latest_genomic.gtf_STAR

for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    sra=${i1}_${i2}
    mkdir $to/$sra/ -p

    fq1=$from_fq/$sra/mix_1.fastq
    fq2=$from_fq/$sra/mix_2.fastq
    outlog=$to/$sra/Log.final.out
    cmd="STAR --runMode alignReads --runThreadN $cpu --genomeDir $ref_star --readFilesIn $fq1 $fq2 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outTmpDir $to/${sra}/STARtmp --outFileNamePrefix $to/$sra/ --limitBAMsortRAM 800000000000"
    if [ -f $outlog ]; then
        continue
    fi
    echo $cmd
done > $runfile
