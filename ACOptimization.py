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


log_file = open('output.log', 'w')
original_stdout = sys.stdout
sys.stdout = Tee(original_stdout, log_file)

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
    print(f"\tcoefficients: {list(coefficients)}")
    print(f"\tP = {' + '.join([f'{x:0.15f}x^{len(poly) - i}' for i, x in enumerate(list(poly))])}")

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

    print(f"Comparing {P} ~ {new_P * RyA_to_Jm}...")
    if P > new_P * RyA_to_Jm:
        print("Updating P")
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

sys.stdout = original_stdout
log_file.close()