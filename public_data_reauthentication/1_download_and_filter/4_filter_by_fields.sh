root=../
workd=$root/resources/
target=$workd/SRA_Accessions_with_expr.uniq.tab

# Field 21: library_strategy
in=$target
out=$target.filter_f21
negative='^(AMPLICON|NULL|OTHER|miRNA-Seq|POOLCLONE|other|FL-cDNA|Synthetic-Long-Read|CLONE|EST|CLONEEND|CTS|SELEX)$'
awk -F"\t" -v x="$negative" 'tolower($21) !~ tolower(x)' $in > $out

# Field 22: library_source
in=$out
out=$target.filter_f22
positive='^(GENOMIC|TRANSCRIPTOMIC)$'
awk -F"\t" -v x="$positive" 'tolower($22) ~ tolower(x)' $in > $out

# Field 23: library_selection
in=$out
out=$target.filter_f23
negative="^(NULL|CAGE|RACE|HMPR|repeat fractionation|padlock probes capture method|MF)$"
awk -F"\t" -v x="$negative" 'tolower($23) !~ tolower(x)' $in > $out

# Field 24: library_layout
in=$out
out=$target.filter_f24
negative="^(NULL)$"
awk -F"\t" -v x="$negative" 'tolower($24) !~ tolower(x)' $in > $out

# Field 24: library_layout
in=$out
out=$target.filter_f25
#negative="NULL|MinION|BGISEQ-500|Ion Torrent Proton|PacBio RS II|"
negative="^(NULL|AB SOLiD System 2.0|AB 5500xl Genetic Analyzer|AB SOLiD 4 System|Helicos HeliScope|AB SOLiD 4 System|AB 3730 Genetic Analyzer|AB 5500 Genetic Analyzer|AB SOLiD System 3.0|AB SOLiD System|AB 3730xL Genetic Analyzer|454 GS Junior|454 GS FLX+|AB 5500xl-W Genetic Analysis System|AB SOLiD 3 Plus System|454 GS 20|454 GS|BGISEQ-50|MinION|PacBio RS II|PacBio RS|PromethION|GridION|Sequel II|Sequel)$"
awk -F"\t" -v x="$negative" 'tolower($25) !~ tolower(x)' $in > $out
