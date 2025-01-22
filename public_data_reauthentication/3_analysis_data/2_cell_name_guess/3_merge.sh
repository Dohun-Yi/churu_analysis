root=../../
workd=$root/analysis/2_cell_name_guess
output1=$workd/cellname_guesse.txt
output2=$workd/cellname_guesse_full.txt

cat $workd/guess/cellname_guesse.txt_* | sort -k1,1V > $output1

# concat multiple json to array
echo "[" > $output2
for f in $(ls $workd/guess/*full* | sort -V); do
    cat $f
    echo ","
done >> $output2
echo "{}]" >> $output2
