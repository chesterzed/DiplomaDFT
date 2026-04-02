from ase.data import atomic_masses_legacy

from subfunctions import *

class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()

datestamp = datetime.now().strftime('%yY%mM%dD_%HH%MM%SS')
os.makedirs(log_dir, exist_ok=True)
os.makedirs("imgs", exist_ok=True)
log_file = open(f'{log_dir}/output_{datestamp}.log', 'w')

original_stdout = sys.stdout
sys.stdout = Tee(original_stdout, log_file)


a = 4.7350490214018
c = 15.5341803094331

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

# demonstrate(structure, name="C36 Laves phase (basic cell)", colors={"Mg":"red","Ni":"blue"})

# суперячейка
supercell = structure.copy()
# supercell.make_supercell([2,2,2])
#
atoms = AseAtomsAdaptor.get_atoms(supercell)
# print("Supercell atoms:", len(atoms))
# print("Supercell atoms:", len(supercell))
# demonstrate(supercell, name="C36 Laves phase (supercell)", colors={"Mg":"red","Ni":"blue"})


os.makedirs("C36", exist_ok=True)
os.makedirs("C36/in", exist_ok=True)
os.makedirs('C36/out', exist_ok=True)
# os.makedirs("C36/xyz", exist_ok=True)

structures = [atoms]

print(atoms)

for i in range(100):
    config = random_displacement(atoms, amplitude=0.01)
    config = strain_cell(config, strain=0.005)
    config = random_swap(config, swap_prob=0.01)
    config = vacancy_defect(config, vac_prob=0.009)
    if is_valid_structure(config):
        structures.append(config)
    # print(config)
    # demonstrate(config, name=f"{i}. C36 Laves phase (basic cell)", colors={"Mg":"red","Ni":"blue"}, save=True)

print("Generated configurations:", len(structures))

for idx, atoms_cfg in enumerate(structures):
    pmg_struct = Structure(
        lattice=atoms_cfg.get_cell(),
        species=[a.symbol for a in atoms_cfg],
        coords=atoms_cfg.get_scaled_positions()
    )

    if idx >= 3:
        method = 'vc-relax'
        in_file = f'C36/in/laves_{method}_{idx:03d}.in'
        out_file = f"C36/out/laves_{method}_{idx:03d}.out"
        pwinput = make_pwinput(pmg_struct, method, f"./C36/out_logs/C36_out_{idx}")
        pwinput.write_file(in_file)
        print("Создан QE input:", in_file)
        run_qe_calculation(in_file, out_file, np=1)

        method = 'relax'
        in_file = f'C36/in/laves_{method}_{idx:03d}.in'
        out_file = f"C36/out/laves_{method}_{idx:03d}.out"
        pwinput = make_pwinput(pmg_struct, method, f"./C36/out_logs/C36_out_{idx}")
        pwinput.write_file(in_file)
        print("Создан QE input:", in_file)
        run_qe_calculation(in_file, out_file, np=1)

        method = 'scf'
        in_file = f'C36/in/laves_{method}_{idx:03d}.in'
        out_file = f"C36/out/laves_{method}_{idx:03d}.out"
        pwinput = make_pwinput(pmg_struct, method, f"./C36/out_logs/C36_out_{idx}", kpoints_grid=(8, 8, 3))
        pwinput.write_file(in_file)
        print("Создан QE input:", in_file)
        run_qe_calculation(in_file, out_file, np=1)


sys.stdout = original_stdout
log_file.close()