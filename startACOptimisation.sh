#!/usr/bin/bash

python ACOptimization.py

#mpirun -np 4 pw.x -inp ACOptimization/in/a_const_000.in > ACOptimization/out/a_const_000.out
for i in ACOptimization/in/*.in; do
    filename=$(basename "$i")
    start_time=$(date '+%Y-%m-%d %H:%M:%S')
    start_seconds=$(date +%s)

    echo "File: $filename"
    echo "Start time: $start_time"

    mpirun -np 4 pw.x -inp "$i" > "ACOptimization/out/${filename%.in}.out"

    end_time=$(date '+%Y-%m-%d %H:%M:%S')
    end_seconds=$(date +%s)
    duration=$((end_seconds - start_seconds))
    echo "Finish time: $end_time"
    echo "Elapsed time: $duration sec"
    echo ""
done