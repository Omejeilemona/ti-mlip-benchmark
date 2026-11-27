"""Build the configurations of the finite-electron-temperature DFT campaign and fill queue/.
Each configuration is computed at k_B*Te = 0.1, 0.5, 1.0, 1.5 and 2.0 eV (Fermi-Dirac width).
Reused configurations are read from ../dft_v3/done (ground-state campaign); new ones are built here.
    python scripts/gen_te_structures.py            (run once, from dft_te/)"""
import os, glob, numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.cluster import Icosahedron
rng = np.random.default_rng(7)
TE = (0.1, 0.5, 1.0, 1.5, 2.0)
SRC = '../dft_v3/done'
cfg = []                                    # (name, atoms)
def add(name, at):
    at = Atoms(at.get_chemical_symbols(), positions=at.get_positions(), cell=at.get_cell(), pbc=at.pbc)
    cfg.append((name, at))
def reuse(stem, new):
    add(new, read(f'{SRC}/{stem}.xyz'))
def rattle(at, s):
    b = at.copy(); b.positions += rng.normal(0, s, b.positions.shape); return b
def omega(a=4.57, c=2.83):
    return Atoms('Ti3', cell=[[a, 0, 0], [-a/2, a*np.sqrt(3)/2, 0], [0, 0, c]], pbc=True,
                 scaled_positions=[[0, 0, 0], [1/3, 2/3, 1/2], [2/3, 1/3, 1/2]])
# A. equation of state of hcp, bcc, fcc, omega (3 volumes each)
base = {'hcp': bulk('Ti', 'hcp', a=2.93, c=4.64), 'bcc': bulk('Ti', 'bcc', a=3.25),
        'fcc': bulk('Ti', 'fcc', a=4.11), 'omega': omega()}
for ph, at in base.items():
    for v in (0.94, 1.00, 1.06):
        b = at.copy(); b.set_cell(at.cell * v**(1/3), scale_atoms=True); add(f'A_eos_{ph}_v{v:.2f}', b)
# B. hcp -> bcc Burgers path (orthorhombic 4-atom cell: lattice strain + shuffle of alternate basal planes)
ah, ch, ab = 2.93, 4.64, 3.25
cell_h = np.diag([ah, np.sqrt(3)*ah, ch]); cell_b = np.diag([ab, np.sqrt(2)*ab, np.sqrt(2)*ab])
fr_h = np.array([[0, 0, 0], [.5, .5, 0], [0, 1/3, .5], [.5, 5/6, .5]])
fr_b = np.array([[0, 0, 0], [.5, .5, 0], [0, .5, .5], [.5, 1.0, .5]])
for i, lam in enumerate(np.linspace(0, 1, 10)):
    cell = (1-lam)*cell_h + lam*cell_b; fr = (1-lam)*fr_h + lam*fr_b
    add(f'B_burgers_{i:02d}', Atoms('Ti4', cell=cell, scaled_positions=fr % 1.0, pbc=True))
# C. hcp -> omega path (reused, every second point), D. rattled omega (reused)
for k, s in enumerate(range(80, 94, 2)): reuse(f'15_path_{s:03d}', f'C_omegapath_{k:02d}')
reuse('60_path_051', 'C_omegapath_07')
for k, s in enumerate((10, 12, 14, 16)): reuse(f'11_omega_{s:03d}', f'D_omega_rattled_{k}')
# E. hot solid (16 atoms), F. liquid (24 atoms), G. surfaces (reused)
for k, s in enumerate((55, 58, 61, 63, 66, 69)): reuse(f'4{0 if s < 63 else 1}_hotsolid_{s:03d}', f'E_hotsolid_{k}')
for k, s in enumerate((35, 40, 45, 50)): reuse(f'30_liquid_{s:03d}', f'F_liquid_{k}')
for k, s in enumerate((71, 75, 79)): reuse(f'50_surface_{s:03d}', f'G_surface_{k}')
# H. small clusters in vacuum (for nanoparticle fragmentation), 12 A box margin
def boxed(at, vac=6.0):
    at.center(vacuum=vac); at.pbc = True; return at
ico = Icosahedron('Ti', noshells=2, latticeconstant=4.11)            # Ti13
add('H_cluster_Ti13_ico', boxed(ico.copy()))
add('H_cluster_Ti13_rattled', boxed(rattle(ico, 0.15)))
add('H_cluster_Ti13_expanded', boxed(Atoms(ico.get_chemical_symbols(), positions=ico.positions*1.12)))
hcp = bulk('Ti', 'hcp', a=2.93, c=4.64).repeat((5, 5, 4)); hcp.center()
c0 = hcp.positions.mean(0); d = np.linalg.norm(hcp.positions - c0, axis=1)
add('H_cluster_Ti19_hcp', boxed(hcp[np.argsort(d)[:19]]))
# I. expanded hot solid (low density, as during expansion/fragmentation)
for k, s in enumerate((55, 66)):
    for rho in (0.70, 0.85):
        b = read(f'{SRC}/4{0 if s < 63 else 1}_hotsolid_{s:03d}.xyz'); b.set_cell(b.cell * rho**(-1/3), scale_atoms=True)
        add(f'I_expanded_{k}_rho{rho:.2f}', b)
# write queue: one file per (configuration, Te); small jobs first within each Te
os.makedirs('queue', exist_ok=True)
for te in TE:
    for name, at in sorted(cfg, key=lambda x: len(x[1])):
        write(f'queue/Te{te:.1f}__{len(at):03d}__{name}.xyz', at, format='extxyz')
print(len(cfg), 'configurations x', len(TE), 'Te =', len(cfg)*len(TE), 'jobs')
