root=../../
workd=$root/analysis/7_snp_profile_of_nonccle2ccle/
match=$root/analysis/4_check_misidentified
institute=$root/analysis/3_country_institute

mkdir -p $workd

ln -srf $match/1697001520 $workd/
ln -srf $match/1697001520.compare_result.final.txt $workd/
ln -srf $match/1697001520.compare_result.final.representative.txt $workd/
ln -srf $match/SRA_churu_output_compare.1697001520.meta $workd/
ln -srf $institute/institute_selected.parsed.corrected.txt $workd/
