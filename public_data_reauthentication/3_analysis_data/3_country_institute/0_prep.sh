root=../../
workd=$root/analysis/3_country_institute
resources=$root/resources

mkdir -p $workd

ln -srf $resources/ROR/table $workd/
ln -srf $resources/GEO/parsed_GSE.txt $workd/
ln -srf $resources/biosample_set.clean.geo.txt $workd/ # sample - geo
ln -srf $resources/biosample_set.clean.summary.interest.txt $workd/ # sample - institutes by Owner
ln -srf $resources/bioproject.gse.interest.txt $workd/ # project - gse
ln -srf $resources/bioproject.journals.interest.txt $workd/ # project - journal
ln -srf $resources/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup $workd/ # sra - target
ln -srf $resources/ENA/ena_sample_human.xml $workd/ # ENA - center name

cd $workd
xtract -input ena_samples_human.xml -pattern SAMPLE -element @accession @center_name
cd -
