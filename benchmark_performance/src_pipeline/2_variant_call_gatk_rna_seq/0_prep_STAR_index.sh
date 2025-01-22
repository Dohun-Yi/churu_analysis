root=../../
ref_fasta=$root/reference/hg38/hg38.fa
annotations_gtf=$root/reference/hg38/refseq/GRCh38_latest_genomic.gtf
genome_dir=$root/resources/STAR2_5
genome_tar=$root/resources/STAR2_5.tar.gz
threads=128

cmd="STAR --runMode genomeGenerate --genomeDir ${genome_dir} --genomeFastaFiles ${ref_fasta} --sjdbGTFfile ${annotations_gtf} --runThreadN ${threads} && cd ${genome_dir}/../ && tar -zcvf $(basename $genome_tar) $(basename $genome_dir)"

echo $cmd
