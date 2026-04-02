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
log_file = open(f'{log_dir}/output_{datestamp}.log', 'w')

original_stdout = sys.stdout
sys.stdout = Tee(original_stdout, log_file)

###############################################################

# find files
os.makedirs('C36', exist_ok=True)
os.makedirs('C36/out', exist_ok=True)
in_files = sorted(glob.glob(f"C36/in/laves*.in"))[2:]
print(' '.join(in_files))
energies = []
for in_file in in_files:
    base_name = os.path.basename(in_file).replace('.in', '.out')
    out_file = f"C36/out/{base_name}"
    run_qe_calculation(in_file, out_file, np=1)
    # print(f'getting energies from file {out_file}: ', end=' ')
    # energies.append(float(get_energies_from_file(out_file)))
    # print(energies[-1])

###############################################################

sys.stdout = original_stdout
log_file.close()