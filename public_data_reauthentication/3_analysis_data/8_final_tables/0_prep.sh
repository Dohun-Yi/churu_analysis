root=../../
workd=$root/analysis/8_final_tables/
resources=$root/resources/

mkdir -p $workd $workd/raw $workd/manual $workd/bulk

# BioSample - 203,916
ln -srf $resources/SRA_Accessions_with_expr.uniq.tab.filter_f25_dedup $workd/raw/

# BioSample - 203,331 (Download / decompress fail removed)
ln -srf $workd/../1_check_output/1697001520.parse.table $workd/raw/

# cell name guess - 203,916
ln -srf $workd/../2_cell_name_guess/cellname_by_access_number.txt $workd/raw/
ln -srf $workd/../2_cell_name_guess/cellname_guesse.txt $workd/raw/
ln -srf $workd/../2_cell_name_guess/cellname_guess_all.txt $workd/raw/

# institute name - corrected
ln -srf $workd/../3_country_institute/institute_selected.parsed.corrected.add_time.txt $workd/raw
ln -srf $workd/../3_country_institute/institute_error_parsed_for_manual_correction.xlsx $workd/manual/

# verification
ln -srf $workd/../4_check_misidentified/cello2ccle.table $workd/raw/
ln -srf $workd/../4_check_misidentified/1697001520.compare_result.final.txt $workd/raw
ln -srf $workd/../4_check_misidentified/1697001520.compare_result.final.representative.txt $workd/raw
ln -srf $workd/../4_check_misidentified/manual_inspection_result.2_last_verification.xlsx $workd/manual/

# Journal
ln -srf $workd/../5_journals/sam2journal.txt $workd/raw/

# Bulk data
ln -srf $workd/../2_cell_name_guess/cellname_guesse_full.txt $workd/bulk/
ln -srf $workd/../3_country_institute/institute_selected_full.txt $workd/bulk/
ln -srf $workd/../1_check_output/1697001520.churu.out $workd/bulk/
