"""Basal generalised stacking-fault path (hcp -> intrinsic fault I2): DFT vs POD v4 + Te corrections.
gamma(s) = [F(s) - F(0)] / A, A = sqrt(3)/2 a^2 (one fault per cell). Run from the repository root:
    python3 scripts/gsf_check.py"""
import glob, tempfile, numpy as np
from ase.io import read, write
from lammps import lammps
CORR = [('potentials/te_snap', 'first corr.'), ('potentials/te_snap_res', 'refitted corr.')]
EVA2MJ = 16021.766                                       # eV/A^2 -> mJ/m^2
tmp = tempfile.mkdtemp()
def pe(atoms, b, D):
    write(f'{tmp}/c.data', atoms, format='lammps-data', specorder=['Ti'], masses=True)
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', f'read_data {tmp}/c.data',
        'pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401 snap',
        'pair_coeff * * table potentials/pod_v4/wall.table WALL 1.6',
        'pair_coeff * * pod potentials/pod_v4/Ti_param.pod potentials/pod_v4/Ti_coefficients.pod Ti',
        'pair_coeff * * zbl 22 22', f'pair_coeff * * snap {D}/Ti_te.snapcoeff {D}/Ti_te.snapparam Te{b:02d}', 'run 0'])
    e = L.get_thermo('pe'); L.close(); return e
for T in (0.1, 0.5, 1.0):
    fs = sorted(glob.glob(f'dft_te/done/Te{T}__*__K_gsf_basal_*.xyz'))
    if len(fs) < 2: print(f'kTe {T}: {len(fs)} points done, skipped'); continue
    at = [read(f) for f in fs]; s = np.array([int(f[-6:-4]) / 10 for f in fs])
    A = np.linalg.norm(np.cross(at[0].cell[0], at[0].cell[1]))
    g = {'DFT': np.array([a.info['F'] for a in at])}
    for D, lab in CORR: g[lab] = np.array([pe(a, int(round(T * 10)), D) for a in at])
    g = {k: EVA2MJ * (v - v[0]) / A for k, v in g.items()}
    print(f'\nkTe = {T} eV: gamma (mJ/m^2) along the Shockley partial, s = 0 (hcp) ... 1 (intrinsic fault I2)')
    print('   s    ' + ''.join(f'{k:>16s}' for k in g))
    for i in range(len(s)): print(f'  {s[i]:.1f}  ' + ''.join(f'{g[k][i]:16.0f}' for k in g))
    print('  max   ' + ''.join(f'{g[k].max():16.0f}' for k in g) + '   (unstable fault energy)')
    print('  I2    ' + ''.join(f'{g[k][-1]:16.0f}' for k in g) + '   (intrinsic fault energy)')
