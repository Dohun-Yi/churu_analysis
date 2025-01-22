root=../../
workd=$root/analysis/3_mix_simulation/
path_fastq=$root/raw_data/fastq/
datalist=$root/datalist.txt

mkdir -p $workd/
ln -srf $root/reference/CCLE/22Q4/Model.csv $workd/
ln -srf $root/reference/CCLE/22Q4/OmicsSomaticMutations.csv $workd/
ln -srf $root/reference/CCLE/22Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv $workd/

mix1=SRR8616091 # THP-1
mix2=SRR8615282 # HCT116
mkdir -p $workd/fastq_original $workd/fastq_shuffled $workd/fastq_mix
ln -srf $path_fastq/$mix1* $workd/fastq_original/
ln -srf $path_fastq/$mix2* $workd/fastq_original/

for mix in $mix1 $mix2; do
    echo "shuffling $mix..."
    fq1_i=$workd/fastq_original/${mix}_1.fastq
    fq2_i=$workd/fastq_original/${mix}_2.fastq
    fq1_o=$workd/fastq_shuffled/${mix}_1.fastq
    fq2_o=$workd/fastq_shuffled/${mix}_2.fastq
    paste <(cat $fq1_i) <(cat $fq2_i) | paste - - - - | shuf | awk -F'\t' -v f1=$fq1_o -v f2=$fq2_o '{OFS="\n"; print $1,$3,$5,$7 > f1; print $2,$4,$6,$8 > f2}'
done

fq1_1=$workd/fastq_shuffled/${mix1}_1.fastq
fq1_2=$workd/fastq_shuffled/${mix1}_2.fastq
fq2_1=$workd/fastq_shuffled/${mix2}_1.fastq
fq2_2=$workd/fastq_shuffled/${mix2}_2.fastq
maxread=80000000

# Divide
for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    n_line1=$(expr $maxread \* 4 \* $i1 / 10)
    n_line2=$(expr $maxread \* 4 \* $i2 / 10)
    echo $n_line1 $n_line2
    folder=$workd/fastq_mix/${i1}_${i2}/
    mkdir -p $folder
    head -n $n_line1 $fq1_1 > $folder/`basename $fq1_1` &
    head -n $n_line1 $fq1_2 > $folder/`basename $fq1_2` &
    head -n $n_line2 $fq2_1 > $folder/`basename $fq2_1` &
    head -n $n_line2 $fq2_2 > $folder/`basename $fq2_2` &
done
wait

# Concat
for i1 in $(seq 0 10); do
    i2=$(expr 10 - $i1 )
    echo $n_line1 $n_line2
    folder=$workd/fastq_mix/${i1}_${i2}/
    mkdir -p $folder
    cat $folder/`basename $fq1_1` $folder/`basename $fq2_1` > $folder/mix_1.fastq &
    cat $folder/`basename $fq1_2` $folder/`basename $fq2_2` > $folder/mix_2.fastq &
done
wait
