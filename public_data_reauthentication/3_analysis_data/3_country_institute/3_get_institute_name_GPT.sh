batch=2000 # 100 per batch
src=3_get_institute_name_GPT.py
runner=3_get_institute_name_GPT.run.sh
for i in $(seq 0 $batch); do
        echo "python $src $i $batch"
done > $runner
