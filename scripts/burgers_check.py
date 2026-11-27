"""Burgers path hcp->bcc: DFT free energy vs POD v4 + Te correction, on the same configurations."""
import sys, glob, os, tempfile, numpy as np
from ase.io import read, write
from lammps import lammps
D = sys.argv[1] if len(sys.argv) > 1 else 'potentials/te_snap'
tmp = tempfile.mkdtemp()
def pe(atoms, b):
    write(f'{tmp}/c.data', atoms, format='lammps-data', specorder=['Ti'], masses=True)
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', f'read_data {tmp}/c.data',
        'pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401 snap',
        'pair_coeff * * table potentials/pod_v4/wall.table WALL 1.6',
        'pair_coeff * * pod potentials/pod_v4/Ti_param.pod potentials/pod_v4/Ti_coefficients.pod Ti',
        'pair_coeff * * zbl 22 22',
        f'pair_coeff * * snap {D}/Ti_te.snapcoeff {D}/Ti_te.snapparam Te{b:02d}', 'run 0'])
    e = L.get_thermo('pe') / len(atoms); L.close(); return e
print('Burgers path, F - F(hcp) in meV/atom:  DFT / potential')
for T in (0.1, 0.5, 1.0, 1.5, 2.0):
    fs = sorted(glob.glob(f'dft_te/done/Te{T}__*__B_burgers_*.xyz'))
    if not fs: continue
    at = [read(f) for f in fs]
    dft = np.array([a.info['F'] / len(a) for a in at]); pot = np.array([pe(a, int(round(T * 10))) for a in at])
    dft, pot = 1000 * (dft - dft[0]), 1000 * (pot - pot[0])
    print(f'kTe {T:.1f} eV  ' + '  '.join(f'{d:6.1f}/{p:6.1f}' for d, p in zip(dft, pot)))
