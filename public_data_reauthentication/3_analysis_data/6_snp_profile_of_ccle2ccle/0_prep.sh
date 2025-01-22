root=../../
workd=$root/analysis/6_snp_profile_of_ccle2ccle
churu=$root/analysis/1_check_output
match=$root/analysis/4_check_misidentified
ccle_snp=$root/reference/CCLE/22Q4/OmicsSomaticMutations.csv
ccle_expr=$root/reference/CCLE/22Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv

mkdir -p $workd

ln -srf $churu/1697001520 $workd/
ln -srf $match/1697001520.compare_result.final.txt $workd/
ln -srf $match/1697001520.compare_result.final.representative.txt $workd/
ln -srf $match/cello2ccle.table $workd/
ln -srf $match/SRA_churu_output_compare.1697001520.meta $workd/
ln -srf $match/Model.tsv $workd/
ln -srf $ccle_snp $workd/
ln -srf $ccle_expr $workd/
