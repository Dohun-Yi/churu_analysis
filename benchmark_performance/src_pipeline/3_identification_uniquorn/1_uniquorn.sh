root=../../
workd=$root/processed_data/3_identification_uniquorn/
datalist=$root/datalist.txt

identifier=./uniquorn.R

runner=./1_uniquorn.run.sh
while read sra sample lib layout; do
    mkdir $workd/$sra/ -p
    vcf=$workd/$sra/$sra.vcf.gz
    identity=$workd/$sra/identity.out
    cmd="Rscript $identifier $vcf $identity"
    echo $cmd

done < $datalist > $runner
