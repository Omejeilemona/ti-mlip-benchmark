"""Static properties of Ti for a given interatomic potential (LAMMPS python).
usage: python props.py <snap|meam> <potdir>
Outputs JSON: hcp a, c/a, cohesive energies (hcp,bcc,fcc,omega), elastic constants (hcp), vacancy formation energy.
"""
import sys, json, numpy as np
from lammps import lammps
kind, potdir = sys.argv[1], sys.argv[2]

def pot_cmds():
    if kind == 'snap':
        return ['pair_style hybrid/overlay zbl 1.8 2.3 snap', 'pair_coeff * * zbl 22 22',
                f'pair_coeff * * snap {potdir}/Ti_snap.snapcoeff {potdir}/Ti_snap.snapparam Ti']
    if kind == 'meam':
        return ['pair_style meam/spline', f'pair_coeff * * {potdir}/Ti.meam.spline Ti']
    if kind == 'pod':
        return ['pair_style pod', f'pair_coeff * * {potdir}/Ti_param.pod {potdir}/Ti_coefficients.pod Ti']
    raise SystemExit('unknown')

def new():
    L = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p'])
    return L

def lattice_cmds(phase, a, ca=1.6, rep=(4, 4, 4)):
    nx, ny, nz = rep
    if phase == 'hcp':
        return [f'lattice custom {a} a1 1 0 0 a2 0 {np.sqrt(3)} 0 a3 0 0 {ca} '
                'basis 0 0 0 basis 0.5 0.5 0 basis 0.5 0.8333333333 0.5 basis 0 0.3333333333 0.5',
                f'region b prism 0 {nx} 0 {ny} 0 {nz} 0 0 0', 'create_box 1 b', 'create_atoms 1 box']
    if phase in ('bcc', 'fcc'):
        return [f'lattice {phase} {a}', f'region b prism 0 {nx} 0 {ny} 0 {nz} 0 0 0', 'create_box 1 b', 'create_atoms 1 box']
    if phase == 'omega':  # hexagonal omega, 3 atoms: (0,0,0),(1/3,2/3,1/2),(2/3,1/3,1/2), c/a~0.62
        return [f'lattice custom {a} a1 1 0 0 a2 0 {np.sqrt(3)} 0 a3 0 0 {ca} '
                'basis 0 0 0 basis 0.5 0.5 0 basis 0.5 0.1666666667 0.5 basis 0 0.6666666667 0.5 '
                'basis 0.5 0.8333333333 0.5 basis 0 0.3333333333 0.5',
                f'region b prism 0 {nx} 0 {ny} 0 {nz} 0 0 0', 'create_box 1 b', 'create_atoms 1 box']

def relaxed(phase, a0, ca=1.6, rep=(4, 4, 4), aniso=True):
    L = new(); L.commands_list(lattice_cmds(phase, a0, ca, rep)); L.command('mass 1 47.867'); L.commands_list(pot_cmds())
    L.commands_list(['thermo_style custom step pe press lx ly lz', 'min_modify line quadratic',
                     f'fix r all box/relax {"aniso" if aniso else "iso"} 0.0 vmax 0.001',
                     'minimize 1e-12 1e-10 10000 100000', 'unfix r', 'minimize 1e-12 1e-10 10000 100000'])
    N = L.get_natoms(); pe = L.get_thermo('pe') / N
    box = [L.get_thermo(k) for k in ('lx', 'ly', 'lz')]
    return L, pe, box, N

out = {}
# hcp
L, e_hcp, box, N = relaxed('hcp', 2.93, 1.59, (6, 4, 4))
a = box[0] / 6; c = box[2] / 4 / 1.0
out['hcp_a'] = a; out['hcp_c_over_a'] = (box[2] / 4) / a; out['E_hcp'] = e_hcp
out['V_hcp'] = box[0] * box[1] * box[2] / N
# elastic constants by finite strain on relaxed hcp (x along a, z along c)
def stress(L):
    L.command('run 0')
    return -np.array([L.get_thermo(k) for k in ('pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy')]) * 1e-4  # GPa, tensile +
L.command('thermo_style custom step pe pxx pyy pzz pyz pxz pxy')
L.command('change_box all triclinic')
d = 0.005; C = np.zeros((6, 6))
lx, ly, lz = box
for j in range(6):
    s = []
    for sgn in (+1, -1):
        L.command('write_restart /tmp/_ref.rst') if (j == 0 and sgn == 1) else None
        L2 = new(); L2.command('read_restart /tmp/_ref.rst'); L2.command('mass 1 47.867'); L2.commands_list(pot_cmds())
        L2.command('thermo_style custom step pe pxx pyy pzz pyz pxz pxy')
        e = sgn * d
        cmd = {0: f'change_box all x scale {1+e} remap', 1: f'change_box all y scale {1+e} remap', 2: f'change_box all z scale {1+e} remap',
               3: f'change_box all yz final {e*lz} remap', 4: f'change_box all xz final {e*lz} remap', 5: f'change_box all xy final {e*ly} remap'}[j]
        L2.command(cmd); L2.command('minimize 1e-12 1e-10 10000 100000')  # relax internal coords
        s.append(stress(L2)); L2.close()
    C[:, j] = (s[0] - s[1]) / (2 * d)
C = 0.5 * (C + C.T)
out.update({'C11': C[0, 0], 'C12': C[0, 1], 'C13': C[0, 2], 'C33': C[2, 2], 'C44': C[3, 3], 'C66': C[5, 5]})
L.close()
# other phases
for ph, a0, ca, rep in (('bcc', 3.25, 1, (5, 5, 5)), ('fcc', 4.1, 1, (4, 4, 4)), ('omega', 4.6, 0.62, (3, 2, 4))):
    Lp, e, b, n = relaxed(ph, a0, ca, rep, aniso=(ph == 'omega'))
    out[f'E_{ph}'] = e; out[f'dE_{ph}-hcp_meV'] = (e - e_hcp) * 1000; Lp.close()
# vacancy formation energy (relaxed, constant volume) in 6x4x4 hcp (384 atoms)
L, e0, box, N = relaxed('hcp', out['hcp_a'], out['hcp_c_over_a'], (6, 4, 4))
L.commands_list(['group del id 1', 'delete_atoms group del', 'minimize 1e-12 1e-10 10000 100000'])
ev = L.get_thermo('pe'); out['E_vac_eV'] = ev - (N - 1) * e0; L.close()
print(json.dumps({k: round(float(v), 4) for k, v in out.items()}, indent=1))
