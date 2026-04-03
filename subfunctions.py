from defines import *


def demonstrate(structure, name="C15 structure", colors: dict={"Mg":"red","Cu":"blue"}, save: bool=False):
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
    if save:
        plt.savefig(f"imgs/{name}.png")
    else:
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


################### ACOptimization ###################


def run_qe_calculation(input_file, output_file, np=1):
    cmd = f"""
        export OMPI_MCA_coll=^hcoll ;
        export OMPI_MCA_btl=tcp,self ;
        export OMPI_MCA_btl=^openib
        export PATH=/home/chesterzed/nvidia_hpc_sdk/hpc_sdk/Linux_x86_64/26.1/comm_libs/12.9/hpcx/hpcx-2.25.1/ompi/bin:$PATH ;
        export PATH=/home/chesterzed/nvidia_hpc_sdk/hpc_sdk/Linux_x86_64/26.1/compilers/bin:$PATH ;
        export LIBRARY_PATH=/home/chesterzed/nvidia_hpc_sdk/hpc_sdk/Linux_x86_64/26.1/cuda/13.1/lib64:$LIBRARY_PATH ;
        export LD_LIBRARY_PATH=/home/chesterzed/nvidia_hpc_sdk/hpc_sdk/Linux_x86_64/26.1/comm_libs/12.9/hpcx/hpcx-2.25.1/ompi/lib:$LD_LIBRARY_PATH ;
        mpirun -n {np} ~/Projects/qe-7.2-gpu/bin/pw.x -inp {input_file} > {output_file}
    """

    # cmd = f"mpirun -np {np} pw.x -inp {input_file} > {output_file}"

    print(f"Start time: {datetime.now()}")
    print(f"Filename: {input_file}")
    print(cmd)
    try:
        result = subprocess.run(cmd,
                                shell=True,
                                capture_output=True,
                                executable='/bin/bash',
                                text=True,
                                timeout=360000)

        print(f"Finish time: {datetime.now()}")
        print(f"Exit code: {result.returncode}")

        if result.returncode != 0:
            print(f"Error code: {result.stderr}")
            return False
        return True

    except subprocess.TimeoutExpired:
        print("Timeout: more than 10 hours")
        return False
    except Exception as e:
        print(f"Exception: {e}")
        return False


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
            "Ni": "Ni.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF"
        },
        control={
            "title": 'Optimisation A C for MgNi2',
            "calculation": "relax",
            "forc_conv_thr": 0.001,
            "etot_conv_thr": 0.00001,
            "prefix": f"{name}",
            "outdir": f"./{path}/out",
            "pseudo_dir": "./pseudo",
            # "tstress": True,
            # "tprnfor": True,
            "wf_collect": False
        },
        system={
            "ibrav": 0,
            "ecutwfc": 50,
            "ecutrho": 400,
            "occupations": "smearing",
            "smearing": "m-v",
            "degauss": 0.03,
            "lspinorb": True,
            "noncolin": True,
        },
        electrons={
            "conv_thr": 1e-8,
            "electron_maxstep": 120,
            "mixing_beta": 0.1,
            "mixing_ndim": 20,
            "mixing_mode": "local-TF",
            "diagonalization": 'david'
        },
        kpoints_grid=(3, 3, 1),
        ions={
            "ion_dynamics":'bfgs',
            "upscale": 60.0,
            "trust_radius_max": 0.05,
            "trust_radius_min": 0.001,
        },
        cell={
            "cell_dynamics": 'bfgs',
            "cell_dofree": 'all',
            "press_conv_thr": 0.01,
            "press": 0.0,
        },
    )
    print('done')
    filename = f'{path}/in/{name}_{idx:03d}.in'
    pwinput.write_file(filename)
    print("Создан QE input:", filename)


def get_K_list(param):
    if param == 'iso':
        return np.arange(0.9, 1.11, 0.025).tolist()
    else:
        return [i / 100 for i in range(92, 109, 2)]


def get_param(lnum):
    return PARAMS[lnum % len(PARAMS)]


def get_energies_from_file(file_path):
    search_string = 'Final energy'
    print(file_path)
    with open(file_path, 'r') as f:
        for line in f:
            if search_string in line:
                return line.strip().split()[3]
    return None


def is_valid_structure(atoms, min_dist=1.8):
    dists = atoms.get_all_distances(mic=True)
    return (dists[dists > 0] > min_dist).all()


def make_pwinput(pmg_struct, calculation,
                 outdir="./C36/out_logs/C36_out",
                 kpoints_grid: tuple=(3, 3, 1)
                 ):
    os.makedirs(outdir, exist_ok=True)
    return PWInput(
        structure=pmg_struct,
        pseudo={
            "Mg": "Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF",
            "Ni": "Ni.pbe-spn-kjpaw_psl.1.0.0.UPF"
        },
        control={
            "title": 'C36_MgNi2',
            "calculation": calculation,
            "forc_conv_thr": 0.0001,
            "etot_conv_thr": 0.00001,
            "prefix":"laves_C36",
            "outdir":outdir,
            "pseudo_dir":"./pseudo",
            "tstress": True,
            "tprnfor": True,
            "restart_mode": "from_scratch" if calculation == "vc-relax" else "restart",
            "nstep": 12 if calculation in ['vc-relax', 'relax'] else 1
        },
        system={
            "ibrav": 0,
            "ecutwfc": 80,
            "ecutrho": 480,
            "occupations": "smearing",
            "smearing": "m-v",
            "degauss": 0.03,
            "lspinorb": False,
            "noncolin": False,
        },
        electrons={
            "conv_thr": 1e-5,
            "electron_maxstep": 120,
            "mixing_beta": 0.2,
            "mixing_ndim": 12,
            "mixing_mode": "plain",
            "diagonalization": 'david'
        },
        kpoints_grid=kpoints_grid,
        ions={
            "ion_dynamics": "bfgs",
            "trust_radius_max": 0.05,
        } if calculation in ['vc-relax', 'relax'] else None,
        cell={
            "cell_dynamics": "bfgs"
        } if calculation == "vc-relax" else None,
    )