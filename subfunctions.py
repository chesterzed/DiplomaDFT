import numpy as np
import matplotlib.pyplot as plt


def demonstrate(structure, name="C15 structure", colors: dict={"Mg":"red","Cu":"blue"} ):
    if "Structure" in str(type(structure)):
        coords = structure.cart_coords
        species = [site.specie.symbol for site in structure]
    elif "MSONAtoms" in str(type(structure)):
        coords = structure.get_positions()
        species = [atom.symbol for atom in structure]


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, element in enumerate(species):
        ax.scatter(
            coords[i][0],
            coords[i][1],
            coords[i][2],
            color=colors[element],
            s=80,
            alpha=0.8
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.title(name)
    plt.show()


def random_displacement(atoms, amplitude=0.05):
    new_atoms = atoms.copy()
    new_atoms.positions += amplitude * np.random.randn(len(atoms),3)
    return new_atoms

def strain_cell(atoms, strain=0.02):
    new_atoms = atoms.copy()
    cell = new_atoms.get_cell()
    deformation = np.eye(3) + strain * np.random.randn(3,3)
    new_cell = cell @ deformation
    new_atoms.set_cell(new_cell, scale_atoms=True)
    return new_atoms

def random_swap(atoms, swap_prob=0.05):
    new_atoms = atoms.copy()
    for j in range(len(new_atoms)):
        if np.random.rand() < swap_prob:
            # случайный атом для обмена
            k = np.random.randint(0,len(new_atoms))
            new_atoms[j].symbol, new_atoms[k].symbol = new_atoms[k].symbol, new_atoms[j].symbol
    return new_atoms

def vacancy_defect(atoms, vac_prob=0.02):
    new_atoms = atoms.copy()
    mask = np.random.rand(len(new_atoms)) > vac_prob
    new_atoms = new_atoms[mask]
    return new_atoms