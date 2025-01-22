root=../
workd=$root/resources/
xml=$workd/biosample_set.xml

# Step 0. Clean data (remove &#13; and line breaks if not ended with > character)
input=$xml
output=$workd/biosample_set.clean.xml
cmd="sed '/>$/!N;s/\n/ /;s/&#13;//;' $xml > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi

# Step 1. extract project summary data
input=$output
output=$workd/biosample_set.clean.summary.txt
# long delimiter and seperator is because user specified attribute name and contents 
# doesn't really have any consensus
cmd="xtract -input $input -pattern BioSample -tab '\t' -sep ',' -def '-' -element  BioSample@accession Organism@taxonomy_id BioSample@submission_date Owner/Name -block Attribute -tab ';;DELIMITER;;' -sep '==EQUAL==' -element Attribute@attribute_name,Attribute  > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi

# Step 2. Filter Project of interests
# this file occationally contain TAB in field 4, which is problematic
input=$output
output=$workd//biosample_set.clean.summary.interest.txt
_species=9606 # human
_pattern="cell.line=="
cmd="grep -E $_pattern $input | awk '\$2==$_species' > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi

## Additional steps to parse title and alias of sample
# Step 3. 
input=$workd/biosample_set.clean.xml
output=$workd/biosample_set.clean.titles.txt
cmd="xtract -input $input -pattern BioSample -tab '\t' -sep ',' -def '-' -element BioSample@accession Title -block Id -if @db_label -equals 'Sample name' -element Id > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi

## now the GEO sample data is required (added 20240725)
# Step 4. 
input=$workd/biosample_set.clean.xml
output=$workd/biosample_set.clean.geo.txt
cmd="xtract -input $input -pattern BioSample -element BioSample@accession -block Links -if @type -equals url -and @label -contains GEO -element Link Link@label > $output"
if [ ! -f $output ]; then echo $cmd && eval $cmd; fi
