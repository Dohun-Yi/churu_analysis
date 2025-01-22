root=../../
workd=$root/analysis/3_country_institute/
path_in=$workd/institute_selected.parsed.txt
path_out1=$workd/institute_error1_unidentified.txt
path_out2=$workd/institute_error2_multiple_candidates.txt

awk -F"\t" '$4=="-"' $path_in | cut -f2 | sort | uniq -c | sort -k1,1nr > $path_out1
cut -f2- $path_in | sort | uniq | cut -f1 | sort | uniq -c | awk '$1!=1' | sort -k1,1nr > $path_out2
