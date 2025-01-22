root=../../
workd=$root/analysis/7_snp_profile_of_nonccle2ccle/
match=$workd/1697001520.compare_result.final.representative.txt

claims=$workd/nonccle_claims
records=$workd/nonccle_records
grep "contam-unknown-nonccle" $match | cut -f3 | sort | uniq > $claims
cat $claims | xargs -I{} awk -F"\t" -v x={} '$3==x && $5 != "match"' $match > $records
