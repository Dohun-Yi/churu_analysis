library("Uniquorn")
library("tidyverse")

# CAUTION:: referencce different!!!!!!! (22Q2)
# reference build is required only once in a lifetime
# ccle_mut = '../../reference/CCLE/22Q2/CCLE_mutations.csv'
# ccle_sam = '../../reference/CCLE/22Q2/sample_info.csv'
# initiate_canonical_databases(ccle_file=ccle_mut, ccle_sample_file=ccle_sam, ref_gen='GRCH38')

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript script.R <vcf_file> <file_path>\n")
  quit(status = 1)
}

vcf_file <- args[1]
file_path <- args[2]

ident_result = as_tibble(vcf_file %>% identify_vcf_file(  ref_gen = "GRCH38"))
write.table(ident_result, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
