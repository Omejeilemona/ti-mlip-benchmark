"""Generate the v3 DFT queue: one extxyz file per configuration in queue/.
File names start with a priority number, so workers take the most important first.
Liquid and hot-solid snapshots come from short POD v1 MD runs in LAMMPS.
usage: python scripts/gen_structures.py <pod_v1_dir>"""
import sys, os, numpy as np
from ase import Atoms
from ase.build import bulk, hcp0001, hcp10m10, bcc110
from ase.io import write
pod = sys.argv[1]
rng = np.random.default_rng(2026)
os.makedirs('queue', exist_ok=True)
n = {'k': 0}
def put(prio, cat, at):
    n['k'] += 1
    at = at.copy(); at.pbc = True; at.wrap()
    write(f'queue/{prio:02d}_{cat}_{n["k"]:03d}.xyz', at, format='extxyz')
def omega(a=4.57, c=2.83):
    return Atoms('Ti3', scaled_positions=[(0, 0, 0), (1/3, 2/3, .5), (2/3, 1/3, .5)],
                 cell=[[a, 0, 0], [-a/2, a*3**.5/2, 0], [0, 0, c]], pbc=True)
def rattled(at, s):
    b = at.copy(); b.positions += rng.normal(0, s, b.positions.shape); return b

# 00 alignment references (same settings as everything else)
put(0, 'align', bulk('Ti', 'hcp', a=2.93, c=4.64)); put(0, 'align', bulk('Ti', 'bcc', a=3.25))
put(0, 'align', bulk('Ti', 'fcc', a=4.11)); put(0, 'align', omega())
# 10 omega: volumes and rattled supercells
for sc in (0.94, 0.97, 1.0, 1.03, 1.06): put(10, 'omega', omega(4.57*sc, 2.83*sc))
for s in (0.03, 0.04, 0.05, 0.06, 0.06, 0.07, 0.08, 0.08, 0.09, 0.1): put(11, 'omega', rattled(omega().repeat((1, 1, 2)), s))
# 20 vacancies: hcp 16-1 atoms and bcc 16-1 atoms, rattled
for s in (0.0, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1):
    v = bulk('Ti', 'hcp', a=2.93, c=4.64, orthorhombic=True).repeat((2, 1, 2)); del v[0]; put(20, 'defect', rattled(v, s))
for s in (0.0, 0.05, 0.08, 0.1):
    v = bulk('Ti', 'bcc', a=3.27, cubic=True).repeat((2, 2, 2)); del v[0]; put(21, 'defect', rattled(v, s))
for s in (0.03, 0.05, 0.08):   # octahedral interstitial in hcp (16 + 1 atoms)
    v = bulk('Ti', 'hcp', a=2.93, c=4.64, orthorhombic=True).repeat((2, 1, 2))
    g = np.stack(np.meshgrid(*[np.linspace(0, 1, 13)[:-1]] * 3, indexing='ij'), -1).reshape(-1, 3)
    pts = v.cell.cartesian_positions(g)
    dmin = [np.min(np.linalg.norm(v.cell.cartesian_positions((v.cell.scaled_positions(v.positions - p) + 0.5) % 1 - 0.5), axis=1)) for p in pts]
    v.append('Ti'); v.positions[-1] = pts[int(np.argmax(dmin))]   # largest empty site = octahedral
    put(22, 'defect', rattled(v, s))
# 30 liquid and 40 hot solid: snapshots from POD v1 MD
from lammps import lammps
def md_snaps(lattice_cmds, T, nsnap, every, pre=()):
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p'] + lattice_cmds + ['mass 1 47.867',
        'pair_style pod', f'pair_coeff * * {pod}/Ti_param.pod {pod}/Ti_coefficients.pod Ti',
        'timestep 0.002'] + list(pre) + [ f'velocity all create {T} {int(T)} mom yes rot yes',
        f'fix 1 all npt temp {T} {T} 0.2 iso 0 0 2.0', 'run 3000'])
    out = []
    for k in range(nsnap):
        L.command(f'run {every}')
        x = np.array(L.numpy.extract_atom('x'))[:L.get_natoms()].copy()
        boxlo, boxhi, xy, yz, xz, _, _ = L.extract_box()
        cell = [[boxhi[0]-boxlo[0], 0, 0], [xy, boxhi[1]-boxlo[1], 0], [xz, yz, boxhi[2]-boxlo[2]]]
        out.append(Atoms('Ti%d' % len(x), positions=x - np.array(boxlo), cell=cell, pbc=True))
    L.close(); return out
melt = ['lattice bcc 3.33', 'region b block 0 2 0 2 0 3', 'create_box 1 b', 'create_atoms 1 box']   # 24 atoms
for T in (2000, 2400, 3000, 3500):
    for a in md_snaps(melt, T, 5, 1500, pre=['velocity all create 4000 7 mom yes', 'fix m all nvt temp 4000 4000 0.1', 'run 3000', 'unfix m']):
        put(30, 'liquid', a)
hcp16 = ['lattice custom 2.95 a1 1 0 0 a2 0 1.7320508 0 a3 0 0 1.588 basis 0 0 0 basis 0.5 0.5 0 basis 0.5 0.8333333 0.5 basis 0 0.3333333 0.5',
         'region b prism 0 2 0 1 0 2 0 0 0', 'create_box 1 b', 'create_atoms 1 box']
bcc16 = ['lattice bcc 3.30', 'region b block 0 2 0 2 0 2', 'create_box 1 b', 'create_atoms 1 box']
for T in (800, 1200):
    for a in md_snaps(hcp16, T, 4, 1500): put(40, 'hotsolid', a)
for T in (1400, 1800):
    for a in md_snaps(bcc16, T, 4, 1500): put(41, 'hotsolid', a)
# 50 surfaces (vacuum 12 A, periodic in-plane)
for maker, kw in ((hcp0001, dict(a=2.93, c=4.64, size=(2, 2, 5))), (hcp10m10, dict(a=2.93, c=4.64, size=(2, 2, 4))), (bcc110, dict(a=3.27, size=(2, 2, 5)))):
    s0 = maker('Ti', vacuum=6.0, **kw)
    for s in (0.0, 0.05, 0.08): put(50, 'surface', rattled(s0, s))
# 60 bcc -> omega collapse path (3-atom cell): z shift of the two non-corner atoms
for d in (0.0, 0.0417, 0.0833, 0.125, 0.1667):
    a = 3.25 * 2**0.5; c = 3.25 * 3**0.5 / 2
    put(15, 'path', Atoms('Ti3', scaled_positions=[(0, 0, 0), (1/3, 2/3, 1/3 + d), (2/3, 1/3, 2/3 - d)],
                          cell=[[a, 0, 0], [-a/2, a*3**.5/2, 0], [0, 0, c]], pbc=True))
for ca in (1.50, 1.54, 1.62, 1.66, 1.70): put(15, 'path', bulk('Ti', 'hcp', a=2.93, c=2.93 * ca))
for a0 in (3.05, 3.15, 3.35, 3.45): put(15, 'path', bulk('Ti', 'bcc', a=a0))
files = sorted(os.listdir('queue'))
from collections import Counter
print(len(files), 'configurations:', dict(Counter(f.split('_')[1] for f in files)))
