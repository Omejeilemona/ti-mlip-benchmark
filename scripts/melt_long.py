"""Longer bcc-liquid coexistence run for a firm melting point.
usage: python melt_long.py <potdir> <outprefix> <T_start> [nx ny nz] [nph_ps]
Defaults: 6 x 6 x 24 bcc cells (1,728 atoms), 60 ps NPH coexistence (about 3x longer and 50% larger than the run in the README). Uses pair_style pod."""
import sys, numpy as np
from lammps import lammps
potdir, pre, Tg = sys.argv[1], sys.argv[2], float(sys.argv[3])
nx, ny, nz = (int(x) for x in (sys.argv[4:7] if len(sys.argv) > 6 else (6, 6, 24)))
nph_ps = float(sys.argv[7]) if len(sys.argv) > 7 else 60.0
dt = 0.002; nph = int(nph_ps / dt)
L = lammps(cmdargs=['-log', f'{pre}.log', '-screen', 'none'])
L.commands_list([
    'units metal', 'atom_style atomic', 'boundary p p p', 'lattice bcc 3.31',
    f'region b block 0 {nx} 0 {ny} 0 {nz}', 'create_box 1 b', 'create_atoms 1 box', 'mass 1 47.867',
    'pair_style pod', f'pair_coeff * * {potdir}/Ti_param.pod {potdir}/Ti_coefficients.pod Ti',
    f'timestep {dt}', f'velocity all create {Tg} {int(Tg)} mom yes rot yes',
    'thermo 500', 'thermo_style custom step temp pe press lx ly lz',
    f'fix 1 all npt temp {Tg} {Tg} 0.2 aniso 0 0 2.0', 'run 5000', 'unfix 1',          # 10 ps bcc equilibration
    'variable zmid equal lz/2+zlo', 'region top block INF INF INF INF ${zmid} INF',
    'group liq region top', 'group sol subtract all liq',
    'velocity sol set 0 0 0', 'fix f sol setforce 0 0 0',
    'fix 2 liq nvt temp 4000 4000 0.2', 'run 5000',                                       # melt top half, 10 ps
    f'fix 2 liq nvt temp {Tg} {Tg} 0.2', 'run 2500', 'unfix 2', 'unfix f',              # cool liquid to T_start
    f'velocity sol create {Tg} 777 mom yes rot yes',
    'fix 3 all nph z 0 0 2.0', 'variable T equal temp',
    f'fix av all ave/time 50 1 50 v_T file {pre}_T.dat',
    f'dump d all custom {nph//10} {pre}_traj.dump id x y z',
    f'run {nph}'])
T = np.loadtxt(f'{pre}_T.dat', comments='#')[:, 1]
h = T[len(T) // 2:]
print(f'{pre}: start {Tg:.0f} K -> mean T (last half) = {h.mean():.0f} K, std {h.std():.0f} K')
L.command(f'write_dump all custom {pre}_final.dump id x y z')
