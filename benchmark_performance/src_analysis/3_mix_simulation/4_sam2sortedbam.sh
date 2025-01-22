root=../../
from=$root/analysis/3_mix_simulation/celid/
cpu=30
runfile=4_sam2sortedbam.run.sh

for bam in $(ls $from/*/*sorted*bam)
do
    if [ -f $bam.bai ]; then
        continue
    fi
    echo "samtools index $bam -@ $cpu"
done > $runfile
