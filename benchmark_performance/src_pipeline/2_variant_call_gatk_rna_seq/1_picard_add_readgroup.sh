root=$PWD/../../
from=$root/processed_data/1_mapping_STAR_mm_noquant/
to=$root/processed_data/2_variant_call_gatk_rna_seq/
datalist=$root/datalist.txt

ml java/jdk_13.0.1 
java=java
jar=$root/programs/picard/2.21.4-3/picard.jar

runner=./1_picard_add_readgroup.run.sh
while read sra sample lib layout; do
    mkdir $to/$sra/ -p
    bam_ori=$from/$sra/Aligned.sortedByCoord.out.bam
    bam=$to/$sra/Aligned.sortedByCoord.out.bam

    cmd="java -jar $jar AddOrReplaceReadGroups I=$bam_ori O=$bam RGID=$sra RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sra"
    echo $cmd
done < $datalist > $runner
