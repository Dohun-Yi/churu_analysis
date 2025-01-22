root=../../
workd=$root/analysis/4_check_misidentified
cell_name=$root/analysis/2_cell_name_guess
churu=$root/analysis/1_check_output

mkdir -p $workd

ln -srf $churu/1697001520 $workd/
ln -srf $churu/1697001520.parse.table $workd/
ln -srf $churu/1697001520.parse.filt.table $workd/
ln -srf $churu/SRA_churu_output_compare.1697001520.meta $workd/
ln -srf $churu/cellosaurus.table.reform.acc.cvcl $workd/
ln -srf $cell_name/cellname_by_access_number.txt $workd/
ln -srf $cell_name/cellname_guesse.txt $workd/
ln -srf $cell_name/cellname_guesse_full.txt $workd/
ln -srf $cell_name/cellname_guess_merge.txt $workd/
ln -srf $root/reference/CCLE/22Q4/Model.tsv $workd/
ln -srf $root/reference/CCLE/22Q4/OmicsSomaticMutations.csv $workd/
ln -srf $root/reference/cellosaurus/cellosaurus.txt $workd/
ln -srf $root/analysis/3_country_institute/institute_selected.parsed.corrected.txt $workd/
