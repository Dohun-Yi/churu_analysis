# Global Cell Line Contamination: Regional Variations and Escalating Risks
# Dohun Yi and Jin-Wu Nam

# This repository contains codes used for all analyses of the CHURU study


# Below directory contains the code responsible for benchmarking analylsis on speed and accuracy
benchmark_performance:
    src_pipeline # processing of raw data to variant calling and identification
    src_analysis # comparison of time and accuracy, and mixture simulation


# Below directory contains the code responsible for re-authentication of public data and related analysis
public_data_reauthentication:
    1_download_and_filter # acquisition of public data and parsing, data selelction
    2_authenticate_on_AWS # main re-authentication codes for Amazon Web Service
    3_analysis_data # all related contamination analysis and figures
