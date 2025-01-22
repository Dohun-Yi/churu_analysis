root=$PWD/../../
from=$root/processed_data/1_mapping_STAR_mm_noquant/
to=$root/processed_data/2_variant_call_gatk_rna_seq/
datalist=$root/datalist.txt

ml java/jdk_13.0.1
java=java
jar=$root/programs/cromwell/cromwell-85.jar
conf=$root/programs/cromwell/use_singularity.conf
wdl=$PWD/gatk4-rna-best-practices.wdl
json_t=./template.json

# config
runner=./2_cromwell_gatk_variant_call.run.sh
while read sra sample lib layout; do
    mkdir $to/$sra/ -p
    bam=$to/$sra/Aligned.sortedByCoord.out.bam
    json=$to/$sra/input.json
    sed "s#BAMFILE#$bam#g" $json_t > $json

    vcf=$(ls $to/$sra/cromwell-executions/RNAseq/*/call-VariantFiltration/execution/Aligned.sortedByCoord.out.variant_filtered.vcf.gz.tbi 2>/dev/null) # expected output
    if [ -f $vcf ] && [ ! -z $vcf ] ; then
        continue
    fi
    cmd="cd $to/$sra/ && java -Dconfig.file=$conf -jar $jar run $wdl -i $json"

    echo $cmd
done < $datalist > $runner
