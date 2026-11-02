"""Test potentials on new v3 DFT configurations that were left out of the fit.
Every 5th liquid, every 4th hot-solid and every 4th surface configuration of dft_v3/done
is held out (10 configurations); fit POD on the rest with
    python scripts/heldout_newdata.py split       -> dft_v3/train, dft_v3/test
    DONE_DIR=dft_v3/train python scripts/fit_v3.py
then compare potentials on the held-out set:
    python scripts/heldout_newdata.py eval potentials/pod_v1 potentials/pod_v3
GPAW energies are shifted to the reference energy scale with the same constant as in fit_v3.py."""
import sys, os, glob, shutil, json
import numpy as np
from ase.io import read

SRC = 'dft_v3/done'
if sys.argv[1] == 'split':
    fs = sorted(glob.glob(f'{SRC}/*.xyz'))
    for d in ('dft_v3/train', 'dft_v3/test'):
        shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    for f in fs:
        dst = 'dft_v3/train'
        for cat, step in (('liquid', 5), ('hotsolid', 4), ('surface', 4)):
            sel = [g for g in fs if f'_{cat}_' in g]
            if f in sel and sel.index(f) % step == step - 1:
                dst = 'dft_v3/test'
        shutil.copy(f, dst)
    print(len(os.listdir('dft_v3/train')), 'train,', len(os.listdir('dft_v3/test')), 'test')
    sys.exit()

from ase.calculators.lammpslib import LAMMPSlib
shift = json.load(open('results/fit_v3_settings.json'))['shift_eV_per_atom']
for d in sys.argv[2:]:
    calc = LAMMPSlib(lmpcmds=['pair_style pod', f'pair_coeff * * {d}/Ti_param.pod {d}/Ti_coefficients.pod Ti'],
                     atom_types={'Ti': 1}, keep_alive=True, log_file=None)
    res = {}
    for f in sorted(glob.glob('dft_v3/test/*.xyz')):
        a = read(f); cat = os.path.basename(f).split('_')[1]
        ed = a.get_potential_energy() / len(a) + shift; fd = a.get_forces()
        b = a.copy(); b.calc = calc
        res.setdefault(cat, []).append((b.get_potential_energy() / len(b) - ed, (b.get_forces() - fd).ravel()))
    print(d)
    for cat, v in res.items():
        de = np.array([x[0] for x in v]); df = np.concatenate([x[1] for x in v])
        print(f'  {cat:9s} n={len(v)}  E_rmse = {1000*np.sqrt(np.mean(de**2)):5.1f} meV/atom   '
              f'F_rmse = {np.sqrt(np.mean(df**2)):.3f} eV/A')
