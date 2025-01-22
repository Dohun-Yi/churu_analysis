# beta release version. WGS not supported

# TODO: license

# Default parameters
cpu=8 # number of processors
k=31 # size of k-mer
K=5 # minimal prevalence of k-mer for sequencing error filtration
mpileup_cut=5 # minimum read required for SNP call 
to=churu # output directory
nucl=rna
fq1=
fq2=
bam=
shared=NoSharedMemory

# TODO: modify the name of opts
usage() {
    echo "Usage: $(basename $0) [-1 fq1] [-2 fq2] [-o output]"
    echo "Options:"
    echo "  -1 FILE   first fastq file"
    echo "  -2 FILE   second fastq file"
    echo "  -b FILE   BAM file"
    echo "  -o DIR    output directory (default: $to)"
    echo "  -t INT    number of threads (default: $cpu)"
    echo "  -k INT    size of k-mer (default: $k)"
    echo "  -K INT    k-mers below this count will be considered sequencing error (default: $K)"
    echo "  -c INT    minimum depth for SNP call (default: $mpileup_cut)"
    echo "  -n STR    type of nucleotide. {dna|rna} (default: $nucl)"
    echo "  -g STR    use shared memory for STAR aligner {NoSharedMemory|LoadAndKeep} (default: $shared)"
    echo "Examples:"
    echo "  Case 1: $(basename $0) -1 fq1 -2 fq2  # paired-end input"
    echo "  Case 2: $(basename $0) -1 fq1  # single-end input"
    echo "  Case 3: $(basename $0) -b bam  # bam input"
    fq1=demo/demo_1.fastq
    fq2=demo/demo_2.fastq
    to=demo/output
    echo "Demo:"
    echo "  $(basename $0) -1 $fq1 -2 $fq2 -o $to"
    echo "Requiremetns:"
    echo "  ./resources/hg38.fa"
    echo "  ./resources/kmerdb"
    echo "  ./resources/reftable"
    echo "  ./resources/model.csv"
    echo "  ./resources/tpm.csv"
    echo "  ./resources/variants.csv"
    exit 1
}



# TODO: k is dependent on reference. further correction needed
# TODO: make an indexing code
# TODO: use minimum resources
while getopts ":b:1:2:k:K:c:t:o:n:g:" opt; do
  case $opt in
    t)
      cpu="$OPTARG"
      ;;
    c)
      mpileup_cut="$OPTARG"
      ;;
    K)
      K="$OPTARG"
      ;;
    k)
      k="$OPTARG"
      ;;
    b)
      bam="$OPTARG"
      ;;
    1)
      fq1="$OPTARG"
      ;;
    2)
      fq2="$OPTARG"
      ;;
    o)
      to="$OPTARG"
      ;;
    g)
      shared="$OPTARG"
      ;;
    n)
      nucl=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]')
      echo $nucl
      ;;
    h)
      usage
      ;;
    \?)
      echo >&2 "Invalid option: -$OPTARG"
      usage
      ;;
    :)
      echo >&2 "Option -$OPTARG requires an argument."
      usage
      ;;
    *)
      echo >&2 "Unknown option: $1"
      usage
      ;;
  esac
done

shift $((OPTIND-1))
if [ $# -ne 0 ]; then
    echo >&2 "Error: Unexpected positional argument(s): $*"
    usage
    exit 1
fi

if [[ ! -z "$fq1" ]] && [ -f "$fq1" ]; then
    # FASTQ ipnut
    if [[ ! -z "$fq2" ]] && [ -f "$fq2" ]; then
        # Paired-end
        :
    else
        # Single-end
        fq2=NONE
    fi
    if [[ ! -z "$bam" ]]; then
        echo >&2 "WARNING: FASTQ and BAM is both specified. using FASTQ as input.."
        bam=NONE
    fi
    mode=FASTQ
elif [[ ! -z "$bam" ]] && [ -f "$bam" ]; then
    # BAM input
    if [ ! -f $bam.bai ]; then
        echo >&2 "ERROR: bam index is not found"
        echo >&2 "  expected index name: $bam.bai"
    fi
    mode=BAM
    ref_bed=$resources/common.kmerdb.pos.bed
elif [[ -z "$fq1" ]] || [[ -z "$fq2" ]] || [[ -z "$bam" ]]; then
    echo >&2 "ERROR: input is required."
    echo >&2 ""
    usage
fi

# Check if nucl is set coorrectly
if [ $nucl != "dna" ] && [ $nucl != "rna" ]; then
    echo >&2 "ERROR: -n must be either DNA or RNA"
    usage
fi

if [ $shared != "LoadAndKeep" ] && [ $shared != "NoSharedMemory" ]; then
    echo >&2 "ERROR: -g must be either LoadAndKeep or NoSharedMemory"
    usage
fi

if [ -d $to ]; then
    echo >&2 "ERROR: the output directory already exists. please specify other name"
    exit 1
fi


# path of input
rootdir=$(dirname $0)
resources=$rootdir/resources/
ref_fa=$resources/hg38.fa
# ref_star=$resources/GRCh38_latest_genomic.gtf_STAR
ref_star=$resources/hg38.fa_star
# reftable=$resources/OmicsSomaticMutations.tx.noSJ.span0.twin.k$k.pos.txt # pos
reftable=$resources/rna.kmerdb.bed # pos
# kmerdb=$resources/OmicsSomaticMutations.tx.noSJ.span0.refsnp.k$k.fasta.kmerdb # kmerDB 8
kmerdb=$resources/rna.kmerdb.fa.kmc # kmerDB 8
csv=$resources/OmicsSomaticMutations.csv
model=$resources/Model.csv
expression=$resources/OmicsExpressionProteinCodingGenesTPMLogp1.csv


# path of programs
check_dependency() {
    local program_name=$1
    local program_path=$(command -v "$program_name" 2> /dev/null)
    local rootdir=$(dirname $0)
    if [ -n "$program_path" ]; then
        echo >&2 "checking $1... OK"
        echo "$program_path"
        return 0
    elif [ -f "$rootdir/$program_name" ]; then
        echo >&2 "checking $1... OK"
        echo "$rootdir/$program_name"
        return 0
    else
        echo >&2 "checking $1... ERROR: Dependency '$program_name' not found"
        return 1
    fi
}

# Display parameter values in a professional format
echo >&2 "--------------------------------------------"
echo >&2 "            CHURU Parameters                "
echo >&2 "--------------------------------------------"
echo >&2 "fq1:        $fq1"
echo >&2 "fq2:        $fq2"
echo >&2 "bam:        $bam"
echo >&2 "output dir: $to/"
echo >&2 "mode:       $mode"
echo >&2 "cpu:        $cpu"
echo >&2 "--------------------------------------------"
echo >&2 "nucl:       $nucl"
echo >&2 "c:          $mpileup_cut"
echo >&2 "k:          $k"
echo >&2 "K:          $K"
echo >&2 "g:          $shared"
echo >&2 "--------------------------------------------"
samtools=$(check_dependency samtools)>&1 || exit 1
varscan=$(check_dependency VarScan.v2.4.5.jar)>&1 || exit 1
kmc=$(check_dependency kmc)>&1 || exit 1
kmctools=$(check_dependency kmc_tools)>&1 || exit 1
filterbyAF=$(check_dependency filter_by_AF.py)>&1 || exit 1
whoami=$(check_dependency identify.py)>&1 || exit 1
star=$(check_dependency STAR)>&1 || { [ "$nucl" == "rna" ] && exit 1; }
bwa=$(check_dependency "bwa-mem2" || check_dependency "bwa")>&1 || { [ "$nucl" == "dna" ] && exit 1; }
parallel=$(check_dependency parallel)>&1 || exit 1
java=$(check_dependency java)>&1 || exit 1
awk=$(check_dependency awk)>&1 || exit 1
cut=$(check_dependency cut)>&1 || exit 1
split=$(check_dependency split)>&1 || exit 1
python=$(check_dependency python)>&1 || exit 1
echo >&2 "--------------------------------------------"
$python -c "import scipy, numpy, argparse, csv, math" || { [ ! $? -eq 0 ] && echo "ERROR: python module not found" && exit 1;}
# STAR, parallel, java, awk, cut, split is used, but not packed here


# path of output
fq1f=$to/filtered_1.fastq # filtered read 1
fq2f=$to/filtered_2.fastq # filtered read 2
bam=$to/align/Aligned.sortedByCoord.out.bam # mapped read
bam_tmp=$to/align/Aligned.out.bam # unsorted mapped read
bed_list=$to/bedlist.txt # bed spllit file list
bed_prefix=$to/tmp/refsnppos.split # bed split file prefix
pileup=$to/output.pileup # output of mpileup
vcf=$to/output.pileup.vcf # output of varscan
mkdir -p $to/ $to/align/ $to/tmp


time_start=`date`


if [ "$mode" == "FASTQ" ]; then

    echo "## 1. READ FILTRATION ##"
    echo "# filter start:" `date`
    # 1-1. get sample kmer
    if [ -f $fq2 ]; then # paired
        echo -e "$fq1\n$fq2" > $to/fqlist
    else # single
        echo -e "$fq1" > $to/fqlist
    fi
    cmd="$kmc -t$cpu -k$k -cs100000 -ci$K @$to/fqlist $to/samplekmer $to/" # CS10000, for DB8
    if [ ! -f $to/samplekmer.kmc_pre ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi
    echo $cmd

    # 1-2. intersect with snp kmer
    cmd="$kmctools -t$cpu simple $to/samplekmer $kmerdb intersect $to/kmers_intersect -ocleft"
    if [ ! -f $to/kmers_intersect.kmc_pre ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi

    # 1-3. dump to table
    cmd="$kmctools -t$cpu transform $to/kmers_intersect dump $to/kmers_intersect.txt"
    if [ ! -f $to/kmers_intersect.txt ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi

    # 1-4. prep snp & reference pair kmer
    refsnpkmer=$to/refsnpkmer.txt
    refsnppos=$to/refsnppos.bed
    refsnplog=$to/refsnp.log
    cmd="$python $filterbyAF $to/kmers_intersect.txt $reftable $refsnpkmer $refsnppos > $refsnplog"
    if [ ! -f $to/refsnpkmer.txt ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi

    # 1-5. refsnpkmer to fasta
    cmd="$awk '{print \">\"FNR\"\\n\"\$1}' $refsnpkmer > $refsnpkmer.fa"
    if [ ! -f $refsnpkmer.fa ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi

    # 1-6. fasta to kmc
    cmd="$kmc -t$cpu -k$k -fm -ci1 $refsnpkmer.fa $to/refsnpkmer $to/"
    if [ ! -f $refsnpkmer.kmc_pre ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi

    # 1-7. read collect (KMC)
    time_start_filter=`date`
    fqf=$to/filtered.kmc.fastq
    cmd="$kmctools -t$cpu filter $to/refsnpkmer @$to/fqlist -ci1 $fqf"
    if [ ! -f $fqf ]; then echo $cmd && eval "$cmd"; fi
    if [ $? != 0 ]; then exit $?; fi
    echo "#kmc filter: $time_start_filter ~ "`date`


    echo "## 2. READ MAPPING ##"
    echo "# mapping start:" `date`
    if [[ "$nucl" == "rna" ]]; then
        cmd="$star --genomeLoad $shared --runMode alignReads --runThreadN $cpu --genomeDir $ref_star --readFilesIn $fqf --outSAMattributes All --outSAMtype BAM Unsorted --outTmpDir $to/align/STARtmp --outFileNamePrefix $to/align/ && $samtools sort -@ $cpu $bam_tmp > $bam"
    elif [[ "$nucl" == "dna" ]]; then
        cmd="$bwa mem -t $cpu $ref_fa $fqf | $samtools view -Sbh - -@ $cpu | $samtools sort -@ $cpu > $bam"
    fi
    if [ ! -f $bam ]; then echo $cmd && eval $cmd; fi
    if [ $? != 0 ]; then exit $?; fi
    echo "# mapping end:" `date`

    echo "# indexing start:" `date`
    cmd="$samtools index $bam -@ $cpu"
    if [ ! -f $bam.bai ]; then echo $cmd && eval $cmd; fi
    if [ $? != 0 ]; then exit $?; fi
    echo "# indexing end:" `date`

elif [ "$mode" == "BAM" ]; then
    refsnppos=$ref_bed
fi

echo "## 3. PILEUP ##"
# prep bed file to sequentially splitted files
n_lines=`wc -l $refsnppos | cut -f1 -d" "` # total line of bed file
n=3 # hardcoded
n_line_per_file=`echo "( $n_lines / ( $cpu * $n ) ) + 1" | bc` # line for each split file

echo "# bed prep start:" `date`
cmd="$split -d -a 3 -l $n_line_per_file $refsnppos $bed_prefix" # split
echo $cmd && eval $cmd
if [ $? != 0 ]; then exit $?; fi

cmd="ls "$bed_prefix*" > $bed_list" # save the list of splitted files
eval $cmd
if [ $? != 0 ]; then exit $?; fi
echo "# bed prep end:" `date`

# TODO: consider use of -B parameter (https://varscan.sourceforge.net/germline-calling.html)
echo "# mpileup start:" `date`
cmd="$parallel -j $cpu --colsep '\t' $samtools mpileup $bam -l {1} -f $ref_fa '|' awk \'\\\$4 \>= $mpileup_cut\' :::: $bed_list > $pileup 2> /dev/null"
if [ ! -f $pileup ]; then echo $cmd && eval $cmd; fi
if [ $? != 0 ]; then exit $?; fi
echo "# mpileup end:" `date`


echo "## 4. VARIANT CALL ##"
echo "# varscan start:" `date`
cmd="$java -jar $varscan mpileup2snp $pileup --output-vcf 1 > $vcf"
if [ ! -f $vcf ]; then echo $cmd && eval $cmd; fi
if [ $? != 0 ]; then exit $?; fi
echo "# varscan end:" `date`


echo "## 5. WHO AM I ##"
echo "# identify start:" `date`
cell_pval=$to/churu.out
cmd="$python $whoami $vcf $csv --nucl $nucl --ccle-model-file $model --ccle-expression-file $expression > $cell_pval"
if [ ! -f $cell_pval ]; then echo $cmd && eval $cmd; fi
if [ $? != 0 ]; then exit $?; fi
echo "# identify end:" `date`

time_end=`date`
tic=$(date -d "$time_start" +%s)
toc=$(date -d "$time_end" +%s)
runtime=$((toc - tic))
echo "## start    : $time_start"
echo "## end      : $time_end"
echo "## run time : $runtime sec"
