#!/usr/bin/bash

python ACOptimization.py

mpirun -np 4 pw.x -inp ACOptimization/in/a_const_000.in > ACOptimization/out/a_const_000.out
#for i in ACOptimization/in/*.in; do
#    mpirun -np 4 pw.x -inp "$i" > "ACOptimization/out/$(basename "${i%.in}.out")"
#done