"""Group J of the finite-Te campaign: isolated atom, dimer curve, trimers, tetramers and Ti5-Ti8 clusters in a
12 A vacuum box (1 k-point). Needed for evaporation and fragmentation: the first POD fits contain no
low-coordination or close-contact configurations. Written to queue_add/, then moved into dft_te/queue/."""
import os, numpy as np
from ase import Atoms
from ase.io import write
rng = np.random.default_rng(11)
TE = (0.1, 0.5, 1.0, 1.5, 2.0); B = 12.5
cfg = []
def add(name, pos):
    a = Atoms(f'Ti{len(pos)}', positions=pos, cell=[B, B, B], pbc=True); a.center(); cfg.append((name, a))
add('J_atom', [[0, 0, 0]])
for r in (1.7, 1.9, 2.1, 2.3, 2.6, 3.0, 3.5, 4.2):
    add(f'J_dimer_r{r:.1f}', [[0, 0, 0], [r, 0, 0]])
for r in (2.2, 2.5, 2.9):
    add(f'J_trimer_tri_r{r:.1f}', [[0, 0, 0], [r, 0, 0], [r/2, r*np.sqrt(3)/2, 0]])
add('J_trimer_lin_r2.4', [[0, 0, 0], [2.4, 0, 0], [4.8, 0, 0]])
for r in (2.4, 2.8):
    t = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * r / (2*np.sqrt(2)); add(f'J_tetra_r{r:.1f}', t)
def compact(n, d=2.6):                      # random compact cluster grown atom by atom
    p = [np.zeros(3)]
    while len(p) < n:
        q = p[rng.integers(len(p))] + rng.normal(size=3) * 0.0
        v = rng.normal(size=3); q = q + d * v / np.linalg.norm(v)
        if min(np.linalg.norm(q - x) for x in p) > 0.92 * d: p.append(q)
    return np.array(p)
for n in (5, 6, 7, 8):
    add(f'J_cluster_Ti{n}', compact(n))
os.makedirs('queue_add', exist_ok=True)
for te in TE:
    for name, at in cfg:
        write(f'queue_add/Te{te:.1f}__{len(at):03d}__{name}.xyz', at, format='extxyz')
print(len(cfg), 'configurations x', len(TE), '=', len(cfg) * len(TE), 'jobs')
