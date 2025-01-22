root=../../
workd=$root/analysis/3_country_institute
output=$workd/institute_selected.txt
output_full=$workd/institute_selected_full.txt

cat $workd/select_gpt/institute_cellname_guesse.txt_* | sort -k1,1V > $output

# concat multiple json to array
echo "[" > $output_full
for f in $(ls $workd/select_gpt/institute_selected_full.txt_* | sort -V); do
    cat $f
    echo ","
done >> $output_full
echo "{}]" >> $output_full
