"""Basal generalised stacking fault (hcp -> fcc stacking change) for the finite-Te DFT queue.
Bulk 8-layer hcp cell (8 atoms) whose c vector is tilted by s*b_p, b_p = a/sqrt(3) <1-100> (Shockley partial):
one fault per cell, unrelaxed; s = 0 ... 1 (s = 1: intrinsic fault I2). k_BTe = 0.1, 0.5, 1.0 eV.
Run from dft_te/:  python3 ../scripts/gen_gsf.py    then start the workers (bash run_te.sh 3)."""
import glob, numpy as np
from ase import Atoms
from ase.io import read, write
h = read(glob.glob('done/Te0.1__*__A_eos_hcp_v1.00.xyz')[0])
a = float(np.mean(h.cell.lengths()[:2])); c = float(h.cell.lengths()[2])      # PBE lattice of the EOS cell
print(f'hcp a = {a:.4f} A, c = {c:.4f} A (from the EOS cell)')
NL = 4                                                                          # 4 hcp cells = 8 layers
a1, a2, cz = np.array([a, 0, 0]), np.array([-a / 2, a * np.sqrt(3) / 2, 0]), np.array([0, 0, c])
pos = [i * cz + f for i in range(NL) for f in (0 * a1, a1 / 3 + 2 * a2 / 3 + cz / 2)]
bp = 2 * a1 / 3 + a2 / 3          # A -> C site, |bp| = a/sqrt(3); s = 1 gives ...ABAB|CACA... (I2)
n = 0
for T in (0.1, 0.5, 1.0):
    for i, s in enumerate(np.linspace(0, 1, 11)):
        at = Atoms('Ti8', positions=pos, cell=[a1, a2, NL * cz + s * bp], pbc=True)
        write(f'queue/Te{T:.1f}__008__K_gsf_basal_{i:02d}.xyz', at, format='extxyz'); n += 1
print(f'added {n} jobs to queue/ (s = 0, 0.1, ..., 1 at 0.1, 0.5 and 1.0 eV)')
