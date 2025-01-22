root=../../
workd=$root/analysis/9_doubling_time
cello=$root/reference/cellosaurus/cellosaurus.txt
data=$root/analysis/4_check_misidentified/1697001520.compare_result.final.txt
institute=$root/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt
finaltable=$root/analysis/8_final_tables/final_table.txt

mkdir -p $workd/
ln -srf $cello $workd/
ln -srf $data $workd/
ln -srf $institute $workd/
ln -srf $finaltable $workd/
