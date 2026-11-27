"""Electron heat capacity of hcp Ti from finite-temperature DFT (same settings as the Te campaign):
Mermin free energy F(Te) of the 2-atom hcp cell at fixed positions, kTe = 0.05 ... 2.5 eV in steps of 0.05 eV.
Run from the repository root:  python3 scripts/ce_dft.py   (restartable: finished temperatures are skipped)
Writes ce_hcp_F_vs_Te.txt; then run scripts/ce_analysis.py."""
import os, time, numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
OUT = 'ce_hcp_F_vs_Te.txt'
done = set(np.round(np.loadtxt(OUT, ndmin=2)[:, 0], 3)) if os.path.exists(OUT) else set()
if not done: open(OUT, 'w').write('# kTe(eV) F(eV/atom)\n')
for kT in np.round(np.arange(0.05, 2.51, 0.05), 3):
    if kT in done: continue
    t0 = time.time(); a = bulk('Ti', 'hcp', a=2.93, c=4.64)
    a.calc = GPAW(mode=PW(450), xc='PBE', kpts={'density': 3.5, 'gamma': True},
                  occupations=FermiDirac(kT), nbands='300%', txt=None)
    F = a.get_potential_energy(force_consistent=True) / len(a)
    open(OUT, 'a').write(f'{kT:.3f} {F:.6f}\n')
    print(f'kTe = {kT:.2f} eV: F = {F:.6f} eV/atom ({time.time() - t0:.0f} s)', flush=True)
