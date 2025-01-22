root=../../
workd=$root/analysis/7_snp_profile_of_nonccle2ccle/
input=$workd/nonccle_records.full.txt

output1=$workd/nonccle_records.full.contam_unknown.txt
output2=$workd/nonccle_records.full.others.txt

awk -F"\t" '$5 ~ "contam-unknown-nonccle"' $input > $output1
awk -F"\t" '$5 !~ "contam-unknown-nonccle"' $input > $output2

