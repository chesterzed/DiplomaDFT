import subprocess
import os
from datetime import datetime
from math import sqrt
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.pwscf import PWInput

import glob


def run_qe_calculation(input_file, output_file, np=4):
    cmd = f"mpirun -np {np} pw.x -inp {input_file} > {output_file}"

    print(f"Start time: {datetime.now()}")
    print(f"Filename: {input_file}")
    print(cmd)
    try:
        result = subprocess.run(cmd,
                                shell=True,
                                capture_output=True,
                                text=True,
                                timeout=18000)

        print(f"Finish time: {datetime.now()}")
        print(f"Exit code: {result.returncode}")

        if result.returncode != 0:
            print(f"Error code: {result.stderr}")
            return False
        return True

    except subprocess.TimeoutExpired:
        print("Timeout: more than 5 hours")
        return False
    except Exception as e:
        print(f"Exception: {e}")
        return False


def V(a, c):
    return (sqrt(3)/2) * a**2 * c


def generateConfig(a, c, crds, spcs, path, idx, name="laves"):
    print('lattice', end=', ')
    lattice = Lattice.hexagonal(a, c)
    print('structure', end=', ')
    structure = Structure.from_spacegroup(
        "P6_3/mmc",
        lattice,
        spcs,
        crds,
        coords_are_cartesian=False
    )
    print('atoms', end=', ')
    atoms = AseAtomsAdaptor.get_atoms(structure)
    os.makedirs(f"{path}", exist_ok=True)
    os.makedirs(f"{path}/in", exist_ok=True)
    os.makedirs(f"{path}/out", exist_ok=True)

    print('pmg_struct', end=', ')
    pmg_struct = Structure(
        lattice=atoms.get_cell(),
        species=[a.symbol for a in atoms],
        coords=atoms.get_scaled_positions()
    )
    print('pwinput', end=', ')
    pwinput = PWInput(
        structure=pmg_struct,
        pseudo={
            "Mg": "Mg.rel-pbesol-spnl-kjpaw_psl.1.0.0.UPF",
            "Ni": "Ni.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF"},
        control={
            "title": 'Optimisation A C for MgNi2',
            "calculation": "relax",
            "forc_conv_thr": 0.001,
            "etot_conv_thr": 0.00001,
            "prefix": f"{name}",
            "outdir": f"./{path}/out",
            "pseudo_dir": "./pseudo",
            # "tstress": True,
            # "tprnfor": True,
            "wf_collect": False},
        system={
            "ibrav": 4,
            "celldm(1)": a,
            "celldm(3)": c/a,
            "ecutwfc": 50,
            "ecutrho": 400,
            "occupations": "smearing",
            "smearing": "m-v",
            "degauss": 0.02,
            "lspinorb": True,
            "noncolin": True,},
        electrons={"conv_thr": 1e-8,
                   "electron_maxstep": 120,
                   "mixing_beta": 0.4,
                   "mixing_ndim": 15,
                   "mixing_mode": "plain",
                   "diagonalization": 'david'},
        kpoints_grid=(3, 3, 1),
        ions={
            "ion_dynamics":'bfgs',
              },
        cell=None,
    )
    pwinput.cell_parameters = None
    pwinput.cell = None
    print('done')
    filename = f'{path}/in/{name}_{idx:03d}.in'
    pwinput.write_file(filename)
    print("Создан QE input:", filename)


print('base vars init')
x1_ni = 1/6
species = ["Mg", "Mg", "Ni", "Ni", "Ni"]
coords = [
    [0, 0, 0.094],            # Mg1 (4e)
    [1/3, 2/3, 0.8442],       # Mg2 (4f)
    [1/3, 2/3, 0.12514],      # Ni1 (4f)
    [x1_ni, 2*x1_ni, 0.25],   # Ni2 (6h)
    [0.5, 0, 0]               # Ni3 (6g)
]


input_dir = "ACOptimization/in"
output_dir = "ACOptimization/out"
folderName = "ACOptimization"
a0 = 4.824
c0 = 15.826
k_list = [i/100 for i in range(90, 111, 2)]
P = 1000000
r3 = sqrt(3)
a2 = a0*a0

while P > 0.1:
    print('loop started')
    print('lists generation')
    V0 = V(a0, c0)
    a_list = [sqrt((2*V0*k)/(r3 * c0)) for k in k_list]
    c_list = [(2*V0*k)/(r3 * a2) for k in k_list]
    f = c0 / a0
    iso_a_list = [np.cbrt((2*V0*k)/(r3*f)) for k in k_list]
    iso_list = [(float(a), float(iso_a_list[i] * f)) for i, a in enumerate(iso_a_list)]

    print(iso_list)

    print('config generation')
    print('c const, a list')
    for i, a in enumerate(a_list):
        print(k_list[i], end=' ')
        generateConfig(a, c0, coords, species, folderName, i, "a_const")
        break
    print('\nc list, a const')
    # for i, c in enumerate(c_list):
    #     print(k_list[i], end=' ')
    #     generateConfig(a0, c, coords, species, folderName, i, "c_const")
    print('\nc, a iso')
    # for i, pair in enumerate(iso_list):
    #     print(k_list[i], end=' ')
    #     generateConfig(pair[0], pair[1], coords, species, folderName, i, "iso")

    P = 0
    break

    print('DFT energies')
    for in_file in glob.glob(f"{input_dir}/*.in"):
        base_name = os.path.basename(in_file).replace('.in', '.out')
        out_file = f"{output_dir}/{base_name}"
        run_qe_calculation(in_file, out_file, np=4)
        break

    print('getting energies from file')
    print("drawing diagrams")
    print("saving diagrams")

    print("making polynomials")
    print("looking for minimum")
    print("getting derivative")

    print(f"Comparing ") # todo: add comparisons into output
    print("Updating P")
    print(f"Looking for new average point")


