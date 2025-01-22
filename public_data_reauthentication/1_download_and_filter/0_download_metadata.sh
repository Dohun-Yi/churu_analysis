# Date of download
# 2023/08/22 17:13 # biosample and bioproject
# 2023/11/02 15:30 # pubmed
root=../resources/

cd $root/
mkdir -p $root $root/PUBMED
wget https://ftp.ncbi.nlm.nih.gov/biosample/biosample_set.xml.gz
wget https://ftp.ncbi.nlm.nih.gov/bioproject/bioproject.xml
wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Datalist/NCBI_SRA_Datalist_20230821.gz
wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_20230821.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20230821.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
wget https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Run_Members.tab

gzip -d biosample_set.xml.gz
gzip -d NCBI_SRA_Datalist_20230821.gz
tar zxvf NCBI_SRA_Metadata_20230821.tar.gz
tar zxvf NCBI_SRA_Metadata_Full_20230821.tar.gz
cd -

cd $root/PUBMED/
seq 1 1166 | xargs -I{} wget https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed23n{}.xml.gz &
seq 1167 1530 | xargs -I{} wget https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/pubmed23n{}.xml.gz &
cd -

echo "Go to scimago and manually download csv file: https://www.scimagojr.com/"
