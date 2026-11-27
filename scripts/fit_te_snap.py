import os, glob, json, tempfile
import numpy as np
from ase.io import read, write
from scipy.interpolate import CubicSpline
from lammps import lammps
TE_FIT = (0.1, 0.5, 1.0, 1.5, 2.0); BINS = np.round(np.arange(0.0, 2.01, 0.1), 2)
SNAP = dict(rcutfac=4.8, twojmax=8, rfac0=0.99363, rmin0=0.0, radelem=0.5, wj=1.0)
W_E, W_F, RIDGE = 1000.0, 1.0, 1e-8; EXCLUDE = ('J_',); TEST_EVERY = 7; OUT = 'potentials/te_snap'
def descriptors(atoms, tmp):
    write(f'{tmp}/c.data', atoms, format='lammps-data', specorder=['Ti'], masses=True)
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none', '-nocite'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', 'atom_modify map array sort 0 0',
        f'read_data {tmp}/c.data', 'pair_style zero 6.0', 'pair_coeff * *',
        f'compute snap all snap {SNAP["rcutfac"]} {SNAP["rfac0"]} {SNAP["twojmax"]} {SNAP["radelem"]} {SNAP["wj"]} '
        f'rmin0 {SNAP["rmin0"]} bzeroflag 1 quadraticflag 0 switchflag 1', 'thermo_style custom step pe', 'run 0'])
    A = L.numpy.extract_compute('snap', 0, 2).copy(); N = L.get_natoms()
    ids = L.numpy.extract_atom('id')[:N].copy(); L.close()
    nc = A.shape[1] - 1; E = np.r_[N, A[0, :nc]]
    Fr = A[1:1 + 3 * N, :nc].reshape(N, 3, nc)[np.argsort(ids)].reshape(3 * N, nc)
    return E, np.c_[np.zeros(3 * N), Fr]
data = {}
for f in sorted(glob.glob('dft_te/done/Te*__*.xyz')):
    b = os.path.basename(f)[:-4]; data.setdefault(b.split('__', 2)[2], {})[float(b.split('__')[0][2:])] = read(f)
names = sorted(n for n in data if 0.1 in data[n])
fitnames = [n for n in names if not n.startswith(EXCLUDE)]
test = set(fitnames[TEST_EVERY - 1::TEST_EVERY])
print(f'{len(names)} configurations; {len(fitnames)} fitted, {len(test)} held out')
tmp = tempfile.mkdtemp(); desc = {n: descriptors(data[n][0.1], tmp) for n in names}
def system(kT, sub):
    rows, y, w = [], [], []
    for n in sub:
        if kT not in data[n]: continue
        a0, a1 = data[n][0.1], data[n][kT]; N = len(a0); E, F = desc[n]
        rows.append(E); y.append(a1.info['F'] - a0.info['E0']); w.append(np.sqrt(W_E) / N)
        d = (a1.get_forces() - a0.get_forces()).ravel(); rows.extend(F); y.extend(d); w.extend([np.sqrt(W_F)] * len(d))
    return np.array(rows), np.array(y), np.array(w)
def fit(kT, sub):
    A, y, w = system(kT, sub); Aw = A * w[:, None]
    return np.linalg.solve(Aw.T @ Aw + RIDGE * np.eye(A.shape[1]), Aw.T @ (y * w))
def errors(beta, kT, sub):
    e, f, e0, f0 = [], [], [], []
    for n in sub:
        if kT not in data[n]: continue
        a0, a1 = data[n][0.1], data[n][kT]; N = len(a0); E, F = desc[n]
        dE = a1.info['F'] - a0.info['E0']; d = (a1.get_forces() - a0.get_forces()).ravel()
        e.append((E @ beta - dE) / N); f.extend(F @ beta - d); e0.append(dE / N); f0.extend(d)
    r = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else float('nan')
    return dict(E_meV=1000 * r(e), F=r(f), dE_meV=1000 * r(e0), dF=r(f0), n=len(e))
rep = {'heldout': {}, 'train': {}, 'small_not_fitted': {}}
print('\nHeld-out: kTe  n  dF rmse meV/atom [size]  force rmse eV/A [size]')
for kT in TE_FIT:
    er = errors(fit(kT, [n for n in fitnames if n not in test]), kT, sorted(test)); rep['heldout'][kT] = er
    print(f'  {kT:4.1f} {er["n"]:2d}  {er["E_meV"]:7.1f} [{er["dE_meV"]:7.1f}]   {er["F"]:6.3f} [{er["dF"]:6.3f}]')
betas = {0.0: np.zeros(len(desc[names[0]][0]))}
print('\nFinal fit (training errors)')
for kT in TE_FIT:
    betas[kT] = fit(kT, fitnames); er = errors(betas[kT], kT, fitnames); rep['train'][kT] = er
    sc = errors(betas[kT], kT, [n for n in names if n.startswith(EXCLUDE)]); rep['small_not_fitted'][kT] = sc
    print(f'  {kT:4.1f}  dF {er["E_meV"]:6.1f} meV/atom, F {er["F"]:.3f} eV/A ({er["n"]})   small clusters (not fitted): {sc["E_meV"]:.0f} meV/atom')
print('\nF(phase)-F(hcp) at V0 (meV/atom), relative to ground state: DFT / correction')
pa = lambda n, kT, b=None: ((data[n][kT].info['F'] - data[n][0.1].info['E0']) if b is None else desc[n][0] @ b) / len(data[n][0.1])
for ph in ('bcc', 'fcc', 'omega'):
    nm, h = f'A_eos_{ph}_v1.00', 'A_eos_hcp_v1.00'
    print(f'  {ph:5s}' + ''.join(f'  {kT}: {1000*(pa(nm,kT)-pa(h,kT)):6.1f}/{1000*(pa(nm,kT,betas[kT])-pa(h,kT,betas[kT])):6.1f}' for kT in TE_FIT))
os.makedirs(OUT, exist_ok=True)
Ts = np.array(sorted(betas)); B = np.array([betas[t] for t in Ts])
coef = CubicSpline(Ts, B, axis=0, bc_type='natural')(BINS); coef[0] = 0.0
lines = ['# dF(R;Te) correction for Ti, one SNAP coefficient set per 0.1 eV bin', f'{len(BINS)} {B.shape[1]}']
for t, c in zip(BINS, coef):
    lines.append(f'Te{int(round(t * 10)):02d} {SNAP["radelem"]} {SNAP["wj"]}'); lines += [f'{x:.12e}' for x in c]
open(f'{OUT}/Ti_te.snapcoeff', 'w').write('\n'.join(lines) + '\n')
open(f'{OUT}/Ti_te.snapparam', 'w').write(f'rcutfac {SNAP["rcutfac"]}\ntwojmax {SNAP["twojmax"]}\nrfac0 {SNAP["rfac0"]}\n'
    f'rmin0 {SNAP["rmin0"]}\nbzeroflag 1\nquadraticflag 0\nswitchflag 1\n')
json.dump(rep, open(f'{OUT}/fit_report.json', 'w'), indent=1)
print(f'\nwrote {OUT}/Ti_te.snapcoeff (bins Te00..Te20 = 0..2 eV) and Ti_te.snapparam')
