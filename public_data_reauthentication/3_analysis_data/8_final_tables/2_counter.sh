root=../../
workd=$root/analysis/8_final_tables/

modes="filt1 filt2"
for mode in $modes; do
    echo "> $mode"
    # 1. cell guess
    input1=$workd/$mode/cellname_by_access_number.txt
    n=$(awk -F'\t' '$2!=""' $input1 | wc -l)
    echo "guess - access number: $n"

    input2=$workd/$mode/cellname_guesse.txt
    n=$(awk -F'\t' '$2!="none of these"' $input2 | wc -l)
    echo "guess - chatGPT: $n"

    n=$(cat $input1 $input2 | awk -F'\t' '$2!="none of these" && $2!=""' | cut -f1 | sort | uniq | wc -l)
    echo "guess - any: $n"

    # 2. n of institutes
    input=$workd/$mode/institute_selected.parsed.corrected.add_time.txt
    n=$(cut -f3 $input | sort | uniq | grep -v 'DISCARD' | wc -l)
    echo "institutes: $n"

    n=$(cut -f3 $input | grep -v 'DISCARD' | wc -l)
    echo "sample with institutes: $n"

    # 3. n of journals and others
    input=$workd/$mode/sam2journal.txt
    n=$(cut -f2 $input | sort | uniq | wc -l)
    echo "projects: $n"

    n=$(awk -F"\t" '$3!=""' $input | cut -f3 | sort | uniq | wc -l)
    echo "publications: $n"

    n=$(awk -F"\t" '$3!=""' $input | wc -l)
    echo "sample with publications: $n"

    n=$(awk -F"\t" '$6!=""' $input | cut -f6 | sort | uniq | wc -l)
    echo "journals: $n (selected one for each sample)"
done
