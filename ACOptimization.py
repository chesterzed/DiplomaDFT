import subprocess
import os
from datetime import datetime
from math import sqrt

from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.pwscf import PWInput
import glob


PARAMS = ['iso', 'a_const', 'c_const']
folderName = 'ACOptimization'
input_dir = f"{folderName}/in"
output_dir = f"{folderName}/out"
diagram_output_dir = f'{folderName}/plots'

Ry_to_J = 2.179874099E-18
A_to_m = 10E-10
RyA_to_Jm = 2179.874099

def run_qe_calculation(input_file, output_file, np=4):
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
                                timeout=36000)

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

##########################################################################################


print('base vars init')

species = ["Mg", "Mg", "Ni", "Ni", "Ni"]
coords = [
    [0, 0, 0.094],            # Mg1 (4e)
    [1/3, 2/3, 0.8442],       # Mg2 (4f)
    [1/3, 2/3, 0.12514],      # Ni1 (4f)
    [1/6, 2/6, 0.25],         # Ni2 (6h)
    [0.5, 0, 0]               # Ni3 (6g)
]

# base (for iso)
# a0 = 4.824
# c0 = 15.826
# loop_number = 0

# for c_const
a0 = 4.735049021244854
c0 = 15.534180308918131
loop_number = 1

P = 1000000
r3 = sqrt(3)
a2 = a0*a0

while P > 0.1:
    print('loop started')
    print('lists generation')
    V0 = V(a0, c0)

    a_list = [sqrt((2*V0*k)/(r3 * c0)) for k in get_K_list('c_const')]
    c_list = [(2*V0*k)/(r3 * a2) for k in get_K_list('a_const')]

    f = c0 / a0
    iso_a_list = [np.cbrt((2*V0*k)/(r3*f)) for k in get_K_list('iso')]
    iso_list = [(float(a), float(a * f)) for a in iso_a_list]


    print('draw diagram')
    print([t[1] / t[0] for t in iso_list])
    plt.scatter(a_list, [c0]*len(a_list))
    plt.scatter([a0] * len(a_list), c_list)
    plt.scatter(iso_a_list, [float(iso_a_list[i] * f) for i, a in enumerate(iso_a_list)])
    plt.show()

    print('config generation')
    k_list = get_K_list(get_param(loop_number))
    match get_param(loop_number):
        case 'iso':
            print('\nc, a iso')
            for i, pair in enumerate(iso_list):
                print(k_list[i], end=' ')
                generateConfig(pair[0], pair[1], coords, species, folderName, i, "iso")
        case 'a_const':
            print('\nc var, a const')
            for i, c in enumerate(c_list):
                print(k_list[i], end=' ')
                generateConfig(a0, c, coords, species, folderName, i, "a_const")
        case 'c_const':
            print('c const, a var')
            for i, a in enumerate(a_list):
                print(k_list[i], end=' ')
                generateConfig(a, c0, coords, species, folderName, i, "c_const")


    print('ordering \'in\' files')
    in_files = sorted(glob.glob(f"{input_dir}/{get_param(loop_number)}*.in"))
    print('DFT energies')
    energies = []
    for in_file in in_files:
        base_name = os.path.basename(in_file).replace('.in', '.out')
        out_file = f"{output_dir}/{base_name}"
        run_qe_calculation(in_file, out_file, np=1)
        print(f'getting energies from file {out_file}: ', end=' ')
        energies.append(float(get_energies_from_file(out_file)))
        print(energies[-1])


    print("making polynomials")
    V_list = [V0*k for k in k_list]
    coefficients = np.polyfit(V_list, energies, 4)
    poly = np.poly1d(coefficients)
    V_poly = np.linspace(min(V_list), max(V_list), 100)
    E_poly = poly(V_poly)
    print(f"\tcoefficients: {coefficients}")
    print(f"\tpolynomial: {poly}")

    print("drawing diagrams")
    plt.plot(V_list, energies, 'o', label='Base points')
    plt.plot(V_poly, E_poly, '-', label='polynomial')
    plt.legend()

    print("saving diagrams...")
    os.makedirs(diagram_output_dir, exist_ok=True)
    plt.savefig(os.path.join(diagram_output_dir, f'{get_param(loop_number)}_{loop_number}.png'), dpi=300, bbox_inches='tight')
    plt.show(block=False)
    plt.pause(2)
    plt.close()

    print("looking for minimum")
    result = minimize_scalar(poly, bounds=(min(V_list), max(V_list)), method='bounded')
    Vmin = result.x
    print(f"\tV = {Vmin}")
    print(f"\tE = {result.fun}")

    print("getting derivative for center of range V")
    deriv = np.polyder(poly)
    new_P = deriv(Vmin)
    print(f"\tP = {' + '.join([f'{x:0.15f}x^{len(deriv) - i}' for i, x in enumerate(list(deriv))])} where x is {Vmin}")
    print(f"\tP = {new_P * RyA_to_Jm}")

    print(f"Comparing ")  # todo: add comparisons into output
    print("Updating P")
    if P > new_P * RyA_to_Jm:
        P = new_P * RyA_to_Jm

    print(f"Calculating new a0 c0")
    match get_param(loop_number):
        case 'iso':
            print('\nc, a iso')
            a0 = np.cbrt((2 * Vmin) / (r3 * f))
            c0 = a0 * f
        case 'a_const':
            print('\nc var, a const')
            c0 = (2*Vmin)/(r3 * a2)
        case 'c_const':
            print('c const, a var')
            a0 = sqrt(2*Vmin/(r3 * c0))

    print(
        f"""Loop {loop_number} finished!
Vmin = {Vmin}
Emin = {result.fun}
new A = {a0}
new C = {c0}
new C/A = {c0 / a0}""")

    loop_number += 1