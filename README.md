# DilomaDFT


## links

#### Pseudo potentials:
- Zn.rel-pbe-dnl-kjpaw_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/zn
- Cu.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/cu
- Ni.rel-pbesol-spn-kjpaw_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/ni
- Mg.rel-pbesol-spnl-kjpaw_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/mg

#### .cif files for (Space-group :P6_3/mmc, 194; C14; MgNi2)
- https://www.ctcms.nist.gov/~knc6/jsmol/JVASP-11969.html
- https://www.crystallography.net/cod/2106100.html

just picture:
- http://img.chem.ucl.ac.uk/sgp/large/194az1.htm

### install Quantum Espresson on linux ubuntu 22.04
- https://qe-dft.readthedocs.io/en/latest/setup/install/

### PW Input documentation
- https://www.quantum-espresso.org/Doc/INPUT_PW.html

### Requirements
- python 3.11.9

### Installation guide of Quantum Espresso with GPU acceleration
- soon...

### Useful commands

- Check amount of correct out files
`ls *0_scf.out | grep -rlE "convergence has been achieved in" . | wc -l`

