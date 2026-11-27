"""Constant-Te MD: hcp Ti with hot electrons and a cold lattice (Langevin at T). Run from the repository root:
    python3 scripts/md_constant_te.py --kte 0.5 --mode npt --time 20
npt: zero pressure, fully flexible cell; nvt: fixed cell. Output md_te/kte{kTe}_{mode}.csv"""
import argparse, os, numpy as np
from lammps import lammps
ap = argparse.ArgumentParser()
ap.add_argument('--kte', type=float, required=True); ap.add_argument('--mode', default='nvt', choices=('nvt', 'npt'))
ap.add_argument('--time', type=float, default=20.0, help='ps with hot electrons'); ap.add_argument('--T', type=float, default=300.0)
ap.add_argument('--rep', default='6 4 4', help='replication of the 4-atom orthohexagonal cell (384 atoms)')
ap.add_argument('--seed', type=int, default=4321); ap.add_argument('--root', default='.'); ap.add_argument('--corr', default='potentials/te_snap')
a = ap.parse_args(); r = a.root
b = int(min(max(round(a.kte * 10), 0), 20))
os.makedirs('md_te', exist_ok=True); tag = f'md_te/kte{a.kte:.1f}_{a.mode}' + (f'_T{a.T:g}' if a.T != 300 else '') + ('_res' if 'res' in a.corr else '')
L = lammps(cmdargs=['-log', f'{tag}.log', '-screen', 'none'])
L.commands_list([
    'units metal', 'atom_style atomic', 'boundary p p p',
    'lattice custom 1.0 a1 2.94 0 0 a2 0 5.0922 0 a3 0 0 4.68 basis 0 0 0 basis 0.5 0.5 0 '
    'basis 0.5 0.8333333 0.5 basis 0 0.3333333 0.5',
    f'region box block 0 {a.rep.split()[0]} 0 {a.rep.split()[1]} 0 {a.rep.split()[2]}', 'create_box 1 box',
    'create_atoms 1 box', 'change_box all triclinic', 'mass 1 47.867',
    'pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401 snap',
    f'pair_coeff * * table {r}/potentials/pod_v4/wall.table WALL 1.6',
    f'pair_coeff * * pod {r}/potentials/pod_v4/Ti_param.pod {r}/potentials/pod_v4/Ti_coefficients.pod Ti',
    'pair_coeff * * zbl 22 22',
    f'pair_coeff * * snap {r}/{a.corr}/Ti_te.snapcoeff {r}/{a.corr}/Ti_te.snapparam Te00',
    'neighbor 1.0 bin', 'timestep 0.001', 'compute ack all ackland/atom',
    f'velocity all create {2 * a.T} {a.seed} mom yes rot no dist gaussian',
    f'fix eq all npt temp {a.T} {a.T} 0.1 aniso 0 0 1.0', 'thermo 1000', 'thermo_modify flush yes', 'run 5000', 'unfix eq'])
L.command(f'pair_coeff * * snap {r}/{a.corr}/Ti_te.snapcoeff {r}/{a.corr}/Ti_te.snapparam Te{b:02d}')
L.commands_list([f'fix th all langevin {a.T} {a.T} 0.1 {a.seed + 1} zero yes',
                 'fix ens all nve' if a.mode == 'nvt' else 'fix ens all nph tri 0 0 1.0'])
N = L.get_natoms()
with open(f'{tag}.csv', 'w') as f:
    f.write('t_ps,frac_hcp,frac_bcc,frac_fcc,frac_ico,frac_other,press_GPa,vol_A3_per_atom,pe_eV_per_atom,temp_K\n')
    for k in range(int(round(a.time / 0.1)) + 1):
        if k: L.command('run 100 pre no post no')
        else: L.command('run 0')
        ack = L.numpy.extract_compute('ack', 1, 1)[:N]
        fr = [np.mean(ack == i) for i in (3, 1, 2, 4, 0)]     # ackland: 0 unknown, 1 bcc, 2 fcc, 3 hcp, 4 ico
        f.write(f'{0.1 * k:.1f},' + ','.join(f'{q:.4f}' for q in fr) +
                f',{L.get_thermo("press") / 1e4:.3f},{L.get_thermo("vol") / N:.4f},{L.get_thermo("pe") / N:.5f},{L.get_thermo("temp"):.1f}\n'); f.flush()
L.command(f'write_dump all custom {tag}_final.dump id x y z c_ack modify sort id')
d = np.loadtxt(f'{tag}.csv', delimiter=',', skiprows=1)
print(f'kTe = {a.kte} eV ({a.mode}, T = {a.T:g} K): after {a.time} ps  hcp {d[-1,1]:.2f}  bcc {d[-1,2]:.2f}  fcc {d[-1,3]:.2f}  '
      f'other {d[-1,5]:.2f}  P = {d[-1,6]:.1f} GPa  V = {d[-1,7]:.2f} A3/atom')
