root=../
workd=$root/resources/
target=$workd/SRA_Accessions_with_expr.uniq.tab

fields=$(seq 21 25)
for i in $fields; do
    out=$target.f$i
    echo "cut -f$i $target | sort | uniq -c | sort -k1,1gr > $out &"
done
