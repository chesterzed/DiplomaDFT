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

for i in range(250):
    config = random_displacement(atoms, amplitude=0.005)
    config = strain_cell(config, strain=0.001)
    config = random_swap(config, swap_prob=0.005)
    config = vacancy_defect(config, vac_prob=0.01)
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

    if idx >= 123:
        methods = ['scf']
        for i, method in enumerate(methods):
            in_file = f'C36/in/laves_{idx:03d}_{i}_{method}.in'
            out_file = f"C36/out/laves_{idx:03d}_{i}_{method}.out"
            pwinput = make_pwinput(pmg_struct=pmg_struct,
                                   calculation=method,
                                   outdir=f"./C36/out_logs/C36_out_{idx}",
                                   restart_mode="from_scratch" if method == methods[0] else "restart",
                                   kpoints_grid=(3, 3, 1),
                                   )
            pwinput.write_file(in_file)
            print("Создан QE input:", in_file)
            run_qe_calculation(in_file, out_file, np=1)


sys.stdout = original_stdout
log_file.close()