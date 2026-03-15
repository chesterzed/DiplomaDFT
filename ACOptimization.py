import os
from math import sqrt
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.pwscf import PWInput


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
            "calculation": "vc-relax",
            "forc_conv_thr": 0.001,
            "prefix": f"{name}",
            "outdir": f"./{path}/out",
            "pseudo_dir": "./pseudo",
            "tstress": ".true.",
            "tprnfor": ".true.",
            "wf_collect": False},
        system={"ecutwfc": 50,
                "ecutrho": 400,
                "occupations": "smearing",
                "smearing": "m-v",
                "degauss": 0.02,
                "lspinorb": True,
                "noncolin": True,},
        electrons={"conv_thr": 1e-8,
                   "electron_maxstep": 120,
                   "mixing_beta": 0.1,
                   "mixing_ndim": 15,
                   "mixing_mode": "local-TF",
                   "diagonalization": 'david'},
        kpoints_grid=(4, 4, 1),
        ions={
            "ion_dynamics":'bfgs',
              },
        cell={
            "cell_dynamics": 'bfgs',
            "cell_dofree": 'all',
            },
    )
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

folderName = "ACOptimization"
a0 = 4.824
c0 = 15.826
k_list = [i/100 for i in range(80, 121, 2)]
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
    print('\nc list, a const')
    for i, c in enumerate(c_list):
        print(k_list[i], end=' ')
        generateConfig(a0, c, coords, species, folderName, i, "c_const")
    print('\nc, a iso')
    for i, pair in enumerate(iso_list):
        print(k_list[i], end=' ')
        generateConfig(pair[0], pair[1], coords, species, folderName, i, "iso")
    P = 0





