from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
import matplotlib.pyplot as plt


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
    [1/3, 2/3, 0.0625],   # 4f (A)
    [0, 0, 0],            # 2a (B)
    [0.8333, 0.6666, 0.25] # 6h (B)
]


# генерация структуры
structure = Structure.from_spacegroup(
    "P6_3/mmc",
    lattice,
    species,
    coords
)
for site in structure:
    print(site)

print("Atoms in unit cell:", len(structure))


# визуализация
atoms = AseAtomsAdaptor.get_atoms(structure)
coords = atoms.get_positions()
species = [atom.symbol for atom in atoms]

colors = {"Mg":"red","Zn":"blue"}

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="3d")

for i,s in enumerate(species):
    ax.scatter(coords[i,0], coords[i,1], coords[i,2],
               color=colors[s], s=120)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.title("C14 Laves phase (12 atoms)")
plt.show()