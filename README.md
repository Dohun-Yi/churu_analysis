# Global cell line contamination: regional variations and escalating risks
### Dohun Yi and Jin-Wu Nam

This repository contains code used for all analyses of the CHURU study.

## Directory Structure

---

### `benchmark_performance`
Contains code responsible for benchmarking analysis on speed and accuracy:
- **`src_pipeline`**  
  Processing of raw data to variant calling and identification.
- **`src_analysis`**  
  Comparison of time and accuracy, and mixture simulation.

---

### `public_data_reauthentication`
Contains code responsible for re-authentication of public data and related analysis:
1. **`1_download_and_filter`**  
   Acquisition of public data, parsing, and data selection.
2. **`2_authenticate_on_AWS`**  
   Main re-authentication code for Amazon Web Service.
3. **`3_analysis_data`**  
   All related contamination analysis and figure generation.
