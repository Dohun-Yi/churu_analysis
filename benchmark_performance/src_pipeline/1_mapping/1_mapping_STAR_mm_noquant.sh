root=../../
from=$root/raw_data/fastq/
to=$root/processed_data/1_mapping_STAR_mm_noquant/
datalist=$root/datalist.txt

# config
cpu=80
ref_star=$root/reference/hg38/refseq/GRCh38_latest_genomic.gtf_STAR
ml STAR/2.5.4

while read sra sample lib layout; do
    mkdir $to/$sra/ -p
    if [ $layout == "SINGLE" ]; then
        fq1=$from/${sra}.fastq
        cmd="STAR --runMode alignReads --runThreadN $cpu --genomeDir $ref_star --readFilesIn $fq1 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outTmpDir $to/${sra}/STARtmp --outFileNamePrefix $to/$sra/ --limitBAMsortRAM 800000000000"
    elif [ $layout == "PAIRED" ]; then
        fq1=$from/${sra}_1.fastq
        fq2=$from/${sra}_2.fastq
        cmd="STAR --runMode alignReads --runThreadN $cpu --genomeDir $ref_star --readFilesIn $fq1 $fq2 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outTmpDir $to/${sra}/STARtmp --outFileNamePrefix $to/$sra/ --limitBAMsortRAM 800000000000"
    else
        echo "error"
        exit
    fi
    outlog=$to/$sra/Log.final.out
    if [ -f $outlog ]; then
        continue
    fi
    echo $cmd
done < $datalist
