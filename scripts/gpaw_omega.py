"""Pilot GPAW-PBE set: (1) alignment structures hcp/bcc/fcc, (2) omega-Ti at 5 volumes + 4 rattled.
Writes extxyz with energy, forces, stress (GPAW, same PW500/PBE/FD 0.1 eV settings)."""
import numpy as np, time
from ase import Atoms
from ase.build import bulk
from ase.io import write
from gpaw import GPAW, PW, FermiDirac
def calc(): return GPAW(mode=PW(500), xc='PBE', kpts={'density': 6.0, 'gamma': True}, occupations=FermiDirac(0.1), txt=None)
def omega(a=4.57, c=2.83):
    return Atoms('Ti3', scaled_positions=[(0,0,0),(1/3,2/3,.5),(2/3,1/3,.5)], cell=[[a,0,0],[-a/2,a*3**.5/2,0],[0,0,c]], pbc=True)
jobs = [('align_hcp', bulk('Ti','hcp',a=2.93,c=4.64)), ('align_bcc', bulk('Ti','bcc',a=3.25)), ('align_fcc', bulk('Ti','fcc',a=4.11))]
for s in (0.94, 0.97, 1.0, 1.03, 1.06):
    jobs.append((f'omega_v{s:.2f}', omega(4.57*s, 2.83*s)))
rng = np.random.default_rng(3)
for k in range(4):
    at = omega().repeat((1,1,2)); at.rattle(0.08, seed=int(rng.integers(1e6))); jobs.append((f'omega_rattle{k}', at))
for name, at in jobs:
    t = time.time(); at.calc = calc()
    e = at.get_potential_energy(); f = at.get_forces(); s = at.get_stress(voigt=False)
    at.info['config_type'] = name
    write(f'{name.split("_")[0]}.xyz', at, format='extxyz', append=True)
    print(name, len(at), f'E/atom {e/len(at):.4f}  maxF {np.abs(f).max():.3f}  {time.time()-t:.0f}s', flush=True)
