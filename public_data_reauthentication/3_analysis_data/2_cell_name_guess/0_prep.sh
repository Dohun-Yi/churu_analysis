root=../../
workd=$root/analysis/2_cell_name_guess
source1=$root/analysis/1_check_output
source2=$root/resources/

mkdir -p $workd $workd/guess/
ln -srf $source1/cellosaurus.table.reform.acc $workd/
ln -srf $source2/biosample_set.clean.titles.filt.txt $workd/
ln -srf $source2/biosample_set.clean.summary.interest.filt.txt $workd/
