import numpy as np
import matplotlib.pyplot as plt
import subprocess
from datetime import datetime
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.pwscf import PWInput
from ase.io import write

import os
import sys
import glob
from math import sqrt
from scipy.optimize import minimize_scalar


PARAMS = ['iso', 'a_const', 'c_const']
folderName = 'ACOptimization'
input_dir = f"{folderName}/in"
output_dir = f"{folderName}/out"
diagram_output_dir = f'{folderName}/plots'
log_dir = 'logs'

Ry_to_J = 2.179874099E-18
A_to_m = 10E-10
RyA_to_Jm = 2179.874099