root=../
workd=$root/resources/

# Extract Journals From BioProject
input=$workd/bioproject.xml
output=$workd/bioproject.journals.txt
cmd="xtract -input $input -pattern Package -tab '\t' -sep ',,' -def '-' -element ArchiveID@accession Publication@id Journal/JournalTitle Journal/Year > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi

# Extract GSE access number From BioProject
input=$workd/bioproject.xml
output=$workd/bioproject.gse.txt
cmd="xtract -input $input -pattern Package -tab '\t' -element ArchiveID@accession -block dbXREF -if @db -equals GEO -element ID > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi
