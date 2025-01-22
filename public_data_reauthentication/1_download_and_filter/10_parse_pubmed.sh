root=../
workd=$root/resources/

runner2=11_parse_pubmed.run2.sh
rm -f $runner1 $runner2
for i in $(seq 1 1530); do
    i0=`printf "%04i" $i`
    input=$workd/PUBMED/pubmed23n$i0.xml
    middle=$workd/PUBMED/pubmed23n$i0.clean.xml
    output=$workd/PUBMED/pubmed23n$i0.table
    cmd2="xtract -strict -input $input -pattern PubmedArticle -tab '\t' -sep ',,' -def '-' -element PMID ELocationID ISSN Journal/Title PubDate/Year > $output"
    echo $cmd2 >> $runner2
done

input=bioproject.journals.interest.txt
output=bioproject.journals.interest.pubmedids.txt
cmd="awk -F\"\t\" '\$2!=\"-\" {print \$2}' $input  | sed 's/,,/\n/g' | sort | uniq > $output"
echo $cmd

input=$output
output=bioproject.journals.interest.pubmedids.name.txt
cmd="cat $input | parallel -j 256 grep {} PUBMED/*table > $output"
echo $cmd
