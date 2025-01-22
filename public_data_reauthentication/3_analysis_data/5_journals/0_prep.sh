# NOTE: authentication rates for journal is downloaded from
# https://doi.org/10.1016/j.isci.2020.101698

root=../../
workd=$root/analysis/5_journals/

mkdir -p $workd

ln -srf $root/analysis/4_check_misidentified/1697001520.compare_result.final.txt $workd
