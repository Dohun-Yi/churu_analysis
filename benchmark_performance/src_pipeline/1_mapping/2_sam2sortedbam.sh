root=../../
from=$root/processed_data/1_mapping_STAR_mm_noquant
cpu=30

for bam in $(ls $from/*/*sorted*bam)
do
    if [ -f $bam.bai ]; then
        continue
    fi
    echo "samtools index $bam -@ $cpu"
done
