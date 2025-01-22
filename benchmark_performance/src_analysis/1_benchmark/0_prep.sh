root=../../
workd=$root/analysis/1_benchmark
datalist=$root/datalist.txt

model=$root/reference/CCLE/22Q4/Model.tsv
cell_list_celid=./cellids.celid.txt
cell_list_uniqu=./cellids.uniquorn.txt
cell_list_churu=./cellids.churu.txt

dir_id=$workd/identity/
mkdir -p $dir_id

ln -srf $model $workd/
ln -srf $cell_list_celid $workd/cell_list_celid.txt
ln -srf $cell_list_uniqu $workd/cell_list_uniqu.txt
ln -srf $cell_list_churu $workd/cell_list_churu.txt

while read sra cell lib layout; do
    out_churu=$root/processed_data/3_identification_CHURU/$sra/churu/churu.out
    out_celid=$root/processed_data/3_identification_CELID/$sra/identity.out
    out_uniqu=$root/processed_data/3_identification_uniquorn/$sra/identity.out

    ln -srf $out_churu $dir_id/$sra.churu.out
    ln -srf $out_celid $dir_id/$sra.celid.out
    ln -srf $out_uniqu $dir_id/$sra.uniqu.out
done < $datalist
