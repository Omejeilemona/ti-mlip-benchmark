"""Evaluate a potential on held-out DFT configurations (never used in fitting).
usage: python heldout.py <snap|meam> <potdir> <n_per_group> <out.npz>
Energies compared after removing the mean offset (MEAM has a different energy zero)."""
import sys, numpy as np
from ase import Atoms
from ase.calculators.lammpslib import LAMMPSlib
kind, potdir, n, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
src = sys.argv[5] if len(sys.argv) > 5 else None
D = 'data/dft/261022_AllInitAndBulk_256atoms_with_types_curatedData/'
H = np.load('fit/holdout_ids.npz')
E = np.load(D + 'energy.npy'); F = np.load(D + 'force.npy').reshape(len(E), -1, 3)
C = np.load(D + 'coord.npy').reshape(len(E), -1, 3); B = np.load(D + 'box.npy').reshape(-1, 3, 3)
if kind == 'snap':
    cmds = ['pair_style hybrid/overlay zbl 1.8 2.3 snap', 'pair_coeff * * zbl 22 22',
            f'pair_coeff * * snap {potdir}/Ti_snap.snapcoeff {potdir}/Ti_snap.snapparam Ti']
elif kind == 'pod':
    cmds = ['pair_style pod', f'pair_coeff * * {potdir}/Ti_param.pod {potdir}/Ti_coefficients.pod Ti']
else:
    cmds = ['pair_style meam/spline', f'pair_coeff * * {potdir}/Ti.meam.spline Ti']
calc = LAMMPSlib(lmpcmds=cmds, atom_types={'Ti': 1}, keep_alive=True, log_file=None)
res = {}
rng = np.random.default_rng(1)
for g in ('Init', 'Bulk'):
    ids = rng.choice(H[g], size=min(n, len(H[g])), replace=False)
    ed, ep, fd, fp = [], [], [], []
    for i in ids:
        at = Atoms('Ti%d' % C.shape[1], positions=C[i], cell=B[i], pbc=True); at.calc = calc
        ep.append(at.get_potential_energy() / len(at)); ed.append(E[i] / len(at))
        fp.append(at.get_forces()); fd.append(F[i])
    ed, ep = np.array(ed), np.array(ep); fd, fp = np.concatenate(fd), np.concatenate(fp)
    off = np.mean(ep - ed) if kind == 'meam' else 0.0
    eres = ep - ed - off; fres = (fp - fd).ravel()
    res[g] = dict(n=len(ids), E_rmse_meV=1000 * np.sqrt(np.mean(eres**2)), F_rmse=np.sqrt(np.mean(fres**2)),
                  F_mae=np.mean(np.abs(fres)), F_std_dft=fd.std())
    np.savez(out.replace('.npz', f'_{g}.npz'), ed=ed, ep=ep - off, fd=fd, fp=fp)
    print(g, {k: round(float(v), 4) for k, v in res[g].items()}, flush=True)
