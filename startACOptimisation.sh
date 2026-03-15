#!/usr/bin/bash

python ACOptimization.py

for i in ACOptimization/in/*.in; do
    mpirun -np 4 pw.x -inp "$i" > "ACOptimization/out/$(basename "${i%.in}.out")"
done