root=../../
from=$root/processed_data/2_variant_call_gatk_rna_seq/
workd=$root/processed_data/3_identification_uniquorn
datalist=$root/datalist.txt

while read sra sample lib layout; do
    path1=$from/$sra/cromwell-executions/RNAseq/
    mkdir -p $workd/$sra

    for path in $(ls -d $path1/*); do
        vcf=$path/call-VariantFiltration/execution/Aligned.sortedByCoord.out.variant_filtered.vcf.gz
        if [ ! -f $vcf ]; then
            echo "Warning! vcf file does not exists on $sra"
            continue
        fi
        ln -srf $vcf $workd/$sra/$sra.vcf.gz
    done
done < $datalist
