"""Te correction for bulk MD with POD v4: dF_res(R;Te) = F(R;Te) - E_POD(R), fitted as linear SNAP per Te and
splined onto bins Te01..Te20 (bin Te00 = POD v4 alone). Unlike fit_te_snap.py (target F - E0, used with EAM),
this correction also absorbs the error of POD v4 itself on the finite-Te configurations, including the Burgers path.
Run from the repository root:  python3 scripts/fit_te_snap_res.py  [weight of Burgers/EOS energies, default 20]
Writes potentials/te_snap_res/{Ti_te.snapcoeff, Ti_te.snapparam, fit_report.json}."""
import os, sys, glob, json, tempfile
import numpy as np
from ase.io import read, write
from scipy.interpolate import CubicSpline
from lammps import lammps
TE_FIT = (0.1, 0.5, 1.0, 1.5, 2.0)
BINS = np.round(np.arange(0.0, 2.01, 0.1), 2)
SNAP = dict(rcutfac=4.8, twojmax=8, rfac0=0.99363, rmin0=0.0, radelem=0.5, wj=1.0)
W_E, W_F, RIDGE = 1000.0, 1.0, 1e-8
W_PATH = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0   # extra energy weight for Burgers path and EOS cells
EXCLUDE = ('J_',); TEST_EVERY = 7; OUT = 'potentials/te_snap_res'
POD = ['pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401',
       'pair_coeff * * table potentials/pod_v4/wall.table WALL 1.6',
       'pair_coeff * * pod potentials/pod_v4/Ti_param.pod potentials/pod_v4/Ti_coefficients.pod Ti',
       'pair_coeff * * zbl 22 22']
tmp = tempfile.mkdtemp()
def lmp(atoms, pair):
    write(f'{tmp}/c.data', atoms, format='lammps-data', specorder=['Ti'], masses=True)
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none', '-nocite'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', 'atom_modify map array sort 0 0',
                     f'read_data {tmp}/c.data'] + pair)
    return L
def descriptors(atoms):
    L = lmp(atoms, ['pair_style zero 6.0', 'pair_coeff * *',
        f'compute snap all snap {SNAP["rcutfac"]} {SNAP["rfac0"]} {SNAP["twojmax"]} {SNAP["radelem"]} {SNAP["wj"]} '
        f'rmin0 {SNAP["rmin0"]} bzeroflag 1 quadraticflag 0 switchflag 1', 'thermo_style custom step pe', 'run 0'])
    A = L.numpy.extract_compute('snap', 0, 2).copy(); N = L.get_natoms()
    ids = L.numpy.extract_atom('id')[:N].copy(); L.close()
    nc = A.shape[1] - 1
    Fr = A[1:1 + 3 * N, :nc].reshape(N, 3, nc)[np.argsort(ids)].reshape(3 * N, nc)
    return np.r_[N, A[0, :nc]], np.c_[np.zeros(3 * N), Fr]
def pod(atoms):
    L = lmp(atoms, POD + ['run 0']); N = L.get_natoms()
    e = L.get_thermo('pe'); ids = L.numpy.extract_atom('id')[:N].copy(); f = L.numpy.extract_atom('f')[:N].copy(); L.close()
    return e, f[np.argsort(ids)]
data = {}
for f in sorted(glob.glob('dft_te/done/Te*__*.xyz')):
    b = os.path.basename(f)[:-4]; kT = float(b.split('__')[0][2:]); data.setdefault(b.split('__', 2)[2], {})[kT] = read(f)
names = sorted(n for n in data if 0.1 in data[n] and not n.startswith(EXCLUDE))
test = set(names[TEST_EVERY - 1::TEST_EVERY])
print(f'{len(names)} configurations fitted ({len(test)} held out in the test); path/EOS energy weight x{W_PATH:g}')
desc = {n: descriptors(data[n][0.1]) for n in names}
ref = {n: pod(data[n][0.1]) for n in names}
w_of = lambda n: W_PATH if n.startswith(('A_', 'B_')) else 1.0
def system(kT, sub):
    R, y, w = [], [], []
    for n in sub:
        if kT not in data[n]: continue
        a = data[n][kT]; N = len(a); E, F = desc[n]; ep, fp = ref[n]
        R.append(E); y.append(a.info['F'] - ep); w.append(np.sqrt(W_E * w_of(n)) / N)
        d = (a.get_forces() - fp).ravel(); R.extend(F); y.extend(d); w.extend([np.sqrt(W_F)] * len(d))
    return np.array(R), np.array(y), np.array(w)
def fit(kT, sub):
    A, y, w = system(kT, sub); Aw, yw = A * w[:, None], y * w
    return np.linalg.solve(Aw.T @ Aw + RIDGE * np.eye(A.shape[1]), Aw.T @ yw)
def err(beta, kT, sub):
    e, f, e0 = [], [], []
    for n in sub:
        if kT not in data[n]: continue
        a = data[n][kT]; N = len(a); E, F = desc[n]; ep, fp = ref[n]
        t = a.info['F'] - ep; e.append((E @ beta - t) / N); e0.append(t / N)
        f.extend(F @ beta - (a.get_forces() - fp).ravel())
    r = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else float('nan')
    return dict(E_meV=1000 * r(e), F=r(f), target_meV=1000 * r(e0), n=len(e))
rep = {'heldout': {}, 'train': {}, 'W_PATH': W_PATH}
print(' kTe   held-out E rmse (meV/atom) [target rms]   force rmse (eV/A)')
for kT in TE_FIT:
    er = err(fit(kT, [n for n in names if n not in test]), kT, sorted(test)); rep['heldout'][kT] = er
    print(f'{kT:4.1f}   {er["E_meV"]:7.1f}  [{er["target_meV"]:7.1f}]   {er["F"]:.3f}   ({er["n"]} configs)')
betas = {kT: fit(kT, names) for kT in TE_FIT}
for kT in TE_FIT: rep['train'][kT] = err(betas[kT], kT, names)
os.makedirs(OUT, exist_ok=True)
Ts = np.array(TE_FIT); B = np.array([betas[t] for t in Ts])
coef = np.zeros((len(BINS), B.shape[1])); coef[1:] = CubicSpline(Ts, B, axis=0, bc_type='natural')(BINS[1:])
lines = ['# POD v4 residual + electron-temperature correction for Ti, one SNAP set per 0.1 eV bin (Te00 = POD v4 only)',
         f'{len(BINS)} {B.shape[1]}']
for t, c in zip(BINS, coef):
    lines.append(f'Te{int(round(t * 10)):02d} {SNAP["radelem"]} {SNAP["wj"]}   # k_BTe = {t:.1f} eV')
    lines += [f'{x:.12e}' for x in c]
open(f'{OUT}/Ti_te.snapcoeff', 'w').write('\n'.join(lines) + '\n')
open(f'{OUT}/Ti_te.snapparam', 'w').write(f'rcutfac {SNAP["rcutfac"]}\ntwojmax {SNAP["twojmax"]}\nrfac0 {SNAP["rfac0"]}\n'
                                          f'rmin0 {SNAP["rmin0"]}\nbzeroflag 1\nquadraticflag 0\nswitchflag 1\n')
json.dump(rep, open(f'{OUT}/fit_report.json', 'w'), indent=1)
print(f'wrote {OUT}/Ti_te.snapcoeff')
