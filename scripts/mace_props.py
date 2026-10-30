"""Static properties of Ti with a MACE model, using ASE (same strains as props.py; smaller cells,
48-atom hcp and 144-atom vacancy cell, because MACE on a CPU is slow).
usage: python scripts/mace_props.py <model>
Outputs JSON: hcp a, c/a, elastic constants, bcc/fcc/omega - hcp energies, vacancy formation energy."""
import sys, json, numpy as np
from ase import Atoms
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.units import GPa
from mace.calculators import MACECalculator
calc = MACECalculator(model_paths=sys.argv[1], device='cpu', default_dtype='float64')
def relax(at, cell=True, fmax=1e-3, steps=400):
    at.calc = calc
    BFGS(FrechetCellFilter(at, hydrostatic_strain=False) if cell else at, logfile=None).run(fmax=fmax, steps=steps)
    return at
def omega(a=4.6, c=2.85):
    return Atoms('Ti3', scaled_positions=[(0, 0, 0), (1/3, 2/3, .5), (2/3, 1/3, .5)],
                 cell=[[a, 0, 0], [-a/2, a*3**.5/2, 0], [0, 0, c]], pbc=True)
out = {}
h = relax(bulk('Ti', 'hcp', a=2.93, c=4.66, orthorhombic=True).repeat((3, 2, 2)))      # 48 atoms
N = len(h); e_hcp = h.get_potential_energy() / N
a = h.cell[0, 0] / 3; c = h.cell[2, 2] / 2
out.update(hcp_a=a, hcp_c_over_a=c / a, E_hcp=e_hcp, V_hcp=h.get_volume() / N)
# elastic constants: +-0.5% strain, internal relaxation at fixed cell, central differences
d = 0.005; C = np.zeros((6, 6)); cell0 = h.cell.array.copy()
pairs = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (1, 2), 4: (0, 2), 5: (0, 1)}
for j in range(6):
    s = []
    for sgn in (1, -1):
        eps = np.zeros((3, 3)); p, q = pairs[j]
        if p == q: eps[p, p] = sgn * d
        else: eps[p, q] = eps[q, p] = sgn * d / 2           # engineering shear strain = sgn*d
        b = h.copy(); b.set_cell(cell0 @ (np.eye(3) + eps), scale_atoms=True); relax(b, cell=False, fmax=1e-4)
        s.append(b.get_stress(voigt=True) / GPa)               # ASE: tensile positive; Voigt xx yy zz yz xz xy
    C[:, j] = (s[0] - s[1]) / (2 * d)
C = 0.5 * (C + C.T)
out.update({'C11': C[0, 0], 'C12': C[0, 1], 'C13': C[0, 2], 'C33': C[2, 2], 'C44': C[3, 3], 'C66': C[5, 5]})
for ph, at in (('bcc', bulk('Ti', 'bcc', a=3.25, cubic=True)), ('fcc', bulk('Ti', 'fcc', a=4.1, cubic=True)), ('omega', omega())):
    at = relax(at); e = at.get_potential_energy() / len(at)
    out[f'E_{ph}'] = e; out[f'dE_{ph}-hcp_meV'] = (e - e_hcp) * 1000
# vacancy: 4x3x3 orthorhombic hcp (144 atoms), relaxed at constant volume
v = relax(bulk('Ti', 'hcp', a=a, c=c, orthorhombic=True).repeat((4, 3, 3)), cell=False); n = len(v); e0 = v.get_potential_energy() / n
del v[0]; relax(v, cell=False)
out['E_vac_eV'] = v.get_potential_energy() - (n - 1) * e0
print(json.dumps({k: round(float(x), 4) for k, x in out.items()}, indent=1))
