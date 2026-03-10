import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write

import matplotlib.pyplot as plt


def demonstrate(structure, name="C15 structure"):
    coords = structure.cart_coords
    species = [site.specie.symbol for site in structure]

    colors = {"Mg":"red","Cu":"blue"}

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, element in enumerate(species):
        ax.scatter(
            coords[i][0],
            coords[i][1],
            coords[i][2],
            color=colors[element],
            s=80
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.title(name)
    plt.show()

def random_displacement(atoms, amplitude=0.05):
    new_atoms = atoms.copy()
    displacement = amplitude * np.random.randn(len(atoms),3)
    new_atoms.positions += displacement
    return new_atoms

def strain_cell(atoms, strain=0.02):
    new_atoms = atoms.copy()
    cell = new_atoms.get_cell()
    deformation = np.eye(3) + strain * np.random.randn(3,3)
    new_cell = cell @ deformation
    new_atoms.set_cell(new_cell, scale_atoms=True)
    return new_atoms


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

demonstrate(structure=structure)

supercell = structure.copy()
supercell.make_supercell([2,2,2])
print("super cell done")
print(supercell)

atoms = AseAtomsAdaptor.get_atoms(supercell)
print("got atoms")

structures = []

for i in range(200):
    config = random_displacement(atoms)
    config = strain_cell(config)
    structures.append(config)

print("config")

for i, s in enumerate(structures):
    print(i, s)
    write(f"C15/structure_{i}.xyz", s)