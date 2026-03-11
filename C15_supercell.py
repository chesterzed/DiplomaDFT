import os
from pymatgen.core import Structure, Lattice
from pymatgen.io.pwscf import PWInput
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
from subfunctions import *


# Генерация C15
a = 7.1
lattice = Lattice.cubic(a)

species = ["Mg", "Cu"]
coords = [
    [0,0,0],
    [0.625,0.625,0.625]
]

structure = Structure.from_spacegroup(
    "Fd-3m",
    lattice,
    species,
    coords
)

# Отрисовка
demonstrate(structure=structure, name="C15 structure", colors={"Mg":"red","Cu":"blue"})

# Создание супер ячейки
supercell = structure.copy()
supercell.make_supercell([2,2,2])
print("super cell done")
print(supercell)

# Перевод в другой формат
atoms = AseAtomsAdaptor.get_atoms(supercell)
print("got atoms")
print(atoms)
demonstrate(structure=atoms, name="C15 structure (supercell)", colors={"Mg":"red","Cu":"blue"})

# Создание конфигов (входного файла)
structures = []

for i in range(1200):
    config = random_displacement(atoms)
    config = strain_cell(config)
    config = random_swap(config)
    config = vacancy_defect(config, 0.01)
    structures.append(config)

print("config")

os.makedirs("C15", exist_ok=True)
os.makedirs("C15/in", exist_ok=True)
os.makedirs("C15/xyz", exist_ok=True)


for idx, atoms_cfg in enumerate(structures):
    pmg_struct = Structure(
        lattice=atoms_cfg.get_cell(),
        species=[a.symbol for a in atoms_cfg],
        coords=atoms_cfg.get_scaled_positions()
    )
    pwinput = PWInput(
        structure=pmg_struct,
        pseudo={"Mg":"Mg.rel-pbesol-spnl-kjpaw_psl.1.0.0.UPF","Cu":"Cu.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF"},
        control={"calculation":"scf", "prefix":"laves_C15", "outdir":"./С15_out", "pseudo_dir":"./pseudo", "tstress": ".true.", "tprnfor": ".true."},
        system={"ecutwfc":50, "ecutrho":400, "occupations": "smearing", "smearing": "mp", "degauss": 0.02, "lspinorb": True, "noncolin": True},
        electrons={"conv_thr":1e-8, "electron_maxstep": 200, "mixing_beta": 0.4, "mixing_mode": "plain", "diagonalization": 'david'},
        kpoints_grid=(4,4,4),
        ions=None,
        cell=None
    )

    filename = f'C15/in/laves_{idx:03d}.in'
    pwinput.write_file(filename)
    print("Создан QE input:", filename)

for i, s in enumerate(structures):
    print(i, s)
    write(f"C15/xyz/structure_{i}.xyz", s)