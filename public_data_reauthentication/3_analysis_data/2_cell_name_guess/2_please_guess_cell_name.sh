batch=2000 # 100 per batch
src=2_please_guess_cell_name.py
runner=2_please_guess_cell_name.run.sh
for i in $(seq 0 $batch); do
	echo "python $src $i $batch"
done > $runner
