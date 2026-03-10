import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
import matplotlib.pyplot as plt


def demonstrate(structure, name="C36 structure"):

    coords = structure.cart_coords
    species = [site.specie.symbol for site in structure]

    colors = {"Mg":"red","Ni":"blue"}

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')

    for i, element in enumerate(species):
        ax.scatter(
            coords[i][0],
            coords[i][1],
            coords[i][2],
            color=colors[element],
            s=50,
            alpha=0.8
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


# параметры решётки
a = 4.9
c = 15.9

lattice = Lattice.hexagonal(a, c)

# https://www.ctcms.nist.gov/~knc6/jsmol/JVASP-11969.html
# https://www.crystallography.net/cod/2106100.html
x1_ni = 1/6
species = ["Mg", "Mg", "Ni", "Ni", "Ni"]
coords = [
    [0, 0, 0.094],        # Mg1 (4e)
    [1/3, 2/3, 0.8442],    # Mg2 (4f)
    [1/3, 2/3, 12514],    # Ni1 (4f)
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

demonstrate(structure, "C36 Laves phase")


# суперячейка
supercell = structure.copy()
supercell.make_supercell([2,2,2])

atoms = AseAtomsAdaptor.get_atoms(supercell)
print("Supercell atoms:", len(atoms))
print("Supercell atoms:", len(supercell))
demonstrate(supercell, "C36 Laves phase")


structures = []

for i in range(800):

    config = random_displacement(atoms)
    config = strain_cell(config)

    structures.append(config)


print("Generated configurations:", len(structures))


for i, s in enumerate(structures):
    write(f"C36/structure_{i}.xyz", s)