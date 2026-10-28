"""Align GPAW omega-Ti energies to the training-data energy zero and write fit/XYZ/Omega.xyz.
Shift = E_hcp(reference model) - E_hcp(GPAW), per atom, from data_omega/align.xyz.
The 7 omega configurations are repeated n times (default 12) to weight them in the POD fit."""
import sys
from ase.io import read
n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
E_HCP_REF = -7.8319          # eV/atom, hcp energy of POD v1 (fitted to the VASP-PBE training data)
hcp = [a for a in read('data_omega/align.xyz', ':') if a.info.get('config_type') == 'align_hcp'][0]
shift = E_HCP_REF - hcp.get_potential_energy() / len(hcp)
lines = []
for a in read('data_omega/omega.xyz', ':'):
    e = a.get_potential_energy() + shift * len(a); f = a.get_forces()
    lat = ' '.join(f'{x:.8f}' for x in a.cell.array.ravel())
    lines.append(f'{len(a)}\nLattice="{lat}" Properties=species:S:1:pos:R:3:forces:R:3 energy={e:.8f} pbc="T T T"\n')
    lines += [f'Ti {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {g[0]:.6f} {g[1]:.6f} {g[2]:.6f}\n' for p, g in zip(a.positions, f)]
open('fit/XYZ/Omega.xyz', 'w').write(''.join(lines) * n)
print(f'shift {shift:.4f} eV/atom; wrote {n} copies')
