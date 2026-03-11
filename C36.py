import os
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.pwscf import PWInput
from ase.io import write
from subfunctions import *

# параметры решётки
a = 4.9
c = 15.9

lattice = Lattice.hexagonal(a, c)

x1_ni = 1/6
species = ["Mg", "Mg", "Ni", "Ni", "Ni"]
coords = [
    [0, 0, 0.094],        # Mg1 (4e)
    [1/3, 2/3, 0.8442],    # Mg2 (4f)
    [1/3, 2/3, 0.12514],    # Ni1 (4f)
    [x1_ni, 2*x1_ni, 0.25], # Ni2 (6h)
    [0.5, 0, 0]           # Ni3 (6g)
]

# Mg1 Mg 0 0 0.094 1 0.0
# Mg2 Mg 0.3333 0.6667 0.8442 1 0.0
# Ni1 Ni 0.3333 0.6667 0.12514 1 0.0
# Ni3 Ni 0.16429 0.32858 0.25 1 0.0
# Ni2 Ni 0.5 0 0 1 0.0

structure = Structure.from_spacegroup(
    "P6_3/mmc",
    lattice,
    species,
    coords,
    coords_are_cartesian=False
)

print(structure)
print("Unit cell atoms:", len(structure))

demonstrate(structure, name="C36 Laves phase (basic cell)", colors={"Mg":"red","Ni":"blue"})

# суперячейка
supercell = structure.copy()
supercell.make_supercell([2,2,2])

atoms = AseAtomsAdaptor.get_atoms(supercell)
print("Supercell atoms:", len(atoms))
print("Supercell atoms:", len(supercell))
demonstrate(supercell, name="C36 Laves phase (supercell)", colors={"Mg":"red","Ni":"blue"})


os.makedirs("C36", exist_ok=True)
os.makedirs("C36/in", exist_ok=True)
os.makedirs("C36/xyz", exist_ok=True)

structures = []

for i in range(1200):
    config = random_displacement(atoms)
    config = strain_cell(config)
    config = random_swap(config)
    config = vacancy_defect(config, 0.01)
    structures.append(config)


print("Generated configurations:", len(structures))

for i in range(1200):
    config = random_displacement(atoms)
    config = strain_cell(config)
    config = random_swap(config)
    config = vacancy_defect(config, 0.01)
    structures.append(config)

for idx, atoms_cfg in enumerate(structures):
    pmg_struct = Structure(
        lattice=atoms_cfg.get_cell(),
        species=[a.symbol for a in atoms_cfg],
        coords=atoms_cfg.get_scaled_positions()
    )
    pwinput = PWInput(
        structure=pmg_struct,
        pseudo={
            "Mg":"Mg.rel-pbesol-spnl-kjpaw_psl.1.0.0.UPF", 
            "Ni": "Ni.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF"},
        control={
            "calculation":"scf", 
            "prefix":"laves_C36", 
            "outdir":"./C36_out", 
            "pseudo_dir":"./pseudo", 
            "tstress": ".true.", 
            "tprnfor": ".true."},
        system={"ecutwfc":50, 
                "ecutrho":400, 
                "occupations": "smearing", 
                "smearing": "mp", 
                "degauss": 0.02,
                "lspinorb": True,
                "noncolin": True},
        electrons={"conv_thr":1e-8,
                   "electron_maxstep": 200,
                   "mixing_beta": 0.4,
                   "mixing_mode": "plain",
                   "diagonalization": 'david'},
        kpoints_grid=(4,4,4),
        ions=None,
        cell=None
    )

    filename = f'C36/in/laves_{idx:03d}.in'
    pwinput.write_file(filename)
    print("Создан QE input:", filename)


for i, s in enumerate(structures):
    write(f"C36/xyz/structure_{i}.xyz", s)
