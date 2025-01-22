src=6_make_group_per_pattern.py
runner=6_make_group_per_pattern.run.sh
for i in $(seq 0 178); do
    echo "python3 $src $i"
done > $runner
