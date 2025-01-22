root=../
to=./instance_resource
# target1=$root/resources/SRA_Accessions_with_expr.uniq.tab.filter_f25.libs.selection # downsampled data for pilot test
target1=$root/resources/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup # full data
target2=./scripts
bucket=s3://dohun/
n_instance=40

### prepare sample list and scripts and upload to bucket ###
if [ ! -f $target1 ]; then
    echo "error: required file not found $target1"
    exit 1
fi

if [ ! -d $target2 ]; then
    echo "error: required directory not found $target2"
    exit 1
fi

rm -rf $to/
aws s3 rm $bucket/instance_resource --recursive

mkdir -p $to/
cp $target1 $to/SRA_Accessions.txt
cp -Lrf $target2 $to/
splitter $to/SRA_Accessions.txt $n_instance

aws s3 cp --recursive $to $bucket/instance_resource

### prepare churu code and required reference files and upload to bucket ###
churu=./churu
resources=$churu/resources/

if [ -f $resources/hg38.fa ]; then
    echo "genome file not found: $resources/hg38.fa"; exit 1
fi

for f in Model.csv OmicsExpressionProteinCodingGenesTPMLogp1.csv OmicsSomaticMutations.csv; do
    if [ -f $resources/$f ]; then
        echo "CCLE annotation file not found: $resources/$f"; exit 1
    fi
done

if [ -f $resources/hg38.fa.sa ]; then
    echo "genome index for bwa not found: $resources/hg38.fa.sa"; exit 1
fi

if [ -d $resources/hg38.fa_star ]; then
    echo "genome index for star not found: $resources/hg38.fa_star"; exit 1
fi

aws s3 cp --recursive $churu $bucket/churu
