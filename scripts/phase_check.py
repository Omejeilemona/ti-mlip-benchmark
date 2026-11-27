"""F(phase) - F(hcp) at 0.94/1.00/1.06 V0: DFT vs POD v4 + correction. usage: phase_check.py [corr dir]"""
import sys, glob, tempfile, numpy as np
from ase.io import read, write
from lammps import lammps
D = sys.argv[1] if len(sys.argv) > 1 else 'potentials/te_snap_res'
tmp = tempfile.mkdtemp()
def pe(atoms, b):
    write(f'{tmp}/c.data', atoms, format='lammps-data', specorder=['Ti'], masses=True)
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', f'read_data {tmp}/c.data',
        'pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401 snap',
        'pair_coeff * * table potentials/pod_v4/wall.table WALL 1.6',
        'pair_coeff * * pod potentials/pod_v4/Ti_param.pod potentials/pod_v4/Ti_coefficients.pod Ti',
        'pair_coeff * * zbl 22 22', f'pair_coeff * * snap {D}/Ti_te.snapcoeff {D}/Ti_te.snapparam Te{b:02d}', 'run 0'])
    e = L.get_thermo('pe') / len(atoms); L.close(); return e
print(f'F - F(hcp) in meV/atom at 0.94 / 1.00 / 1.06 V0:  DFT / potential  ({D})')
for T in (0.1, 0.5, 1.0, 1.5, 2.0):
    for ph in ('fcc', 'bcc', 'omega'):
        out = []
        for v in ('0.94', '1.00', '1.06'):
            fp, fh = glob.glob(f'dft_te/done/Te{T}__*__A_eos_{ph}_v{v}.xyz'), glob.glob(f'dft_te/done/Te{T}__*__A_eos_hcp_v{v}.xyz')
            if not fp or not fh: continue
            p, h = read(fp[0]), read(fh[0]); b = int(round(T * 10))
            d = 1000 * (p.info['F'] / len(p) - h.info['F'] / len(h)); m = 1000 * (pe(p, b) - pe(h, b))
            out.append(f'{d:6.1f}/{m:6.1f}')
        if out: print(f'kTe {T:.1f}  {ph:5s}  ' + '   '.join(out))
