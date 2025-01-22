# Select ~1000 samples for pilot test on AWS
root=../
workd=$root/resources/
target=$workd/SRA_Accessions_with_expr.uniq.tab.filter_f25
libs=$target.libs
selection=$target.libs.selection

# if [ ! -f $libs ]; then
#     cut -f21- $target | sort | uniq -c | sort -k1,1gr | awk '$1>30' | cut -f2- > $libs
# fi
cut -f21- $target | sort | uniq -c | sort -k1,1gr | awk '$1>30' | cut -f2- > $libs

IFS=$'\n'
while read line; do
    cmd="shuf $target | grep -F \"$line\" | head -n 2"
    echo -e $cmd 1>&2
    eval $cmd
done < $libs > $selection
