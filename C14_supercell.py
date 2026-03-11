import os
from pymatgen.core import Structure, Lattice
from pymatgen.io.pwscf import PWInput
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
from subfunctions import *


# параметры решётки
a = 4.85
c = 7.85
lattice = Lattice.hexagonal(a, c)

# асимметричная часть C14
species = [
    "Mg",      # 4f
    "Zn",      # 2a
    "Zn"       # 6h
]

coords = [
    [1/3, 2/3, 0.0625],      # 4f (A)
    [0, 0, 0],               # 2a (B)
    [0.8333, 0.6666, 0.25]   # 6h (B)
]

# генерация структуры
structure = Structure.from_spacegroup(
    "P6_3/mmc",
    lattice,
    species,
    coords
)
demonstrate(structure, name="C14 Laves phase (12 atoms)", colors={"Mg":"red","Zn":"blue"})

structure.make_supercell([2, 2, 2])
print("Atoms in unit cell:", len(structure))
atoms = AseAtomsAdaptor.get_atoms(structure)

# визуализация
demonstrate(atoms, name="C14 Laves phase (12 atoms)", colors={"Mg":"red","Zn":"blue"})

# Сохранение структуры в файлы
os.makedirs("C14", exist_ok=True)
os.makedirs("C14/in", exist_ok=True)
os.makedirs("C14/xyz", exist_ok=True)

structures = []

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
        pseudo={"Mg":"Mg.rel-pbe-spnl-kjpaw_psl.1.0.0.UPF", "Zn":"Zn.rel-pbe-dnl-kjpaw_psl.1.0.0.UPF"},
        control={"calculation":"scf", "prefix":"laves_C14", "outdir":"./C14_out", "pseudo_dir":"./pseudo", "tstress": True, "tprnfor": True},
        system={"ecutwfc":40, "ecutrho":240, "occupations": "smearing", "smearing": "mp", "degauss": 0.02, "lspinorb": True, "noncolin": True},
        electrons={"conv_thr":1e-8, "electron_maxstep": 200, "mixing_beta": 0.4, "mixing_mode": "plain", "diagonalization": 'david'},
        kpoints_grid=(3,3,3),
        ions=None,
        cell=None
    )

    filename = f'C14/in/laves_{idx:03d}.in'
    pwinput.write_file(filename)
    print("Создан QE input:", filename)


for i, s in enumerate(structures):
    write(f"C14/xyz/structure_{i}.xyz", s)
