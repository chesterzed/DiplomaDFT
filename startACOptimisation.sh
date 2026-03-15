#!/usr/bin/bash

python ACOptimization.py

mpirun -np 4 pw.x -inp ACOptimization/in/laves_000.in > ACOptimization/out/laves_000.out