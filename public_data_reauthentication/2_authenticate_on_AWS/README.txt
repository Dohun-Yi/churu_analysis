# !! This code involves launching multiple AWS EC2 instances and data transfer to S3 storage bucket
# !! Beforee start, users are expected to have proper AWS setup, and all user-specific parameters in each script should be modified accordingly
# !! Related parameters includes access key, security-group, key-name, bucket and others

# 1. prepare required resources by running the script "0_prep_resource.sh"
#    it will initialize bucket instance resource
#    and upload list of sample access number and related scripts to bucket ("./instance_resource")
#    and full churu (AWS version) code and related references to bucket ("./churu")
#    prerequisites:
#    - the list of samples must be prepared using the scripts in directory "1_download_and_filter"
#    - the references files like genome fasta (hg38.fa), annotation (gencodev34 and CCLE), index file for STAR and BWA-MEM and
#      churu reference built by churu_build.py are expected to be placed on ("./churu/resources/")
#      prepared "./churu/resources/" will look like below.
#      (the example contains index files for bwa-mem2, but index for bwa-mem is permitted)
#      ./churu/resources/:    
#        hg38.fa
#        hg38.fa.0123
#        hg38.fa.amb
#        hg38.fa.ann
#        hg38.fa.bwt
#        hg38.fa.bwt.2bit.64
#        hg38.fa.fai
#        hg38.fa.pac
#        hg38.fa.sa
#        hg38.fa_star
#        Model.csv
#        OmicsExpressionProteinCodingGenesTPMLogp1.csv
#        OmicsSomaticMutations.csv
#        common.kmerdb.pos.bed
#        rna.kmerdb.bed
#        rna.kmerdb.fa.kmc.kmc_pre
#        rna.kmerdb.fa.kmc.kmc_suf

bash 0_prep_resource.sh

# 2. to start analysis on AWS, go to "instance_manager" directory and run "1_start_instance_vol6_250.sh" script
#    again, the account specific informations like access key, security-group, key-name, bucket address and others should be modified
#    file list on "./instance_manager/"
#        1_start_instance_vol6_250.sh: launch AWS instances and start analysis
#        2_stop_all_instance.sh: terminate all instances
#        tracker.sh: monitor progress
#        util_scan_instance.sh: check instance status

cd ./instance_manager/
bash 1_start_instance_vol6_250.sh
