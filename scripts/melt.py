"""Melting point of Ti by solid(bcc)-liquid coexistence (NPH), LAMMPS python.
usage: python melt.py <snap|meam> <potdir> <outprefix> [nx ny nz] [T_guess]"""
import sys, numpy as np
from lammps import lammps
kind, potdir, pre = sys.argv[1:4]
nx, ny, nz = (int(x) for x in (sys.argv[4:7] if len(sys.argv) > 6 else (6, 6, 16)))
Tg = float(sys.argv[7]) if len(sys.argv) > 7 else 1950.0
pot = (['pair_style hybrid/overlay zbl 1.8 2.3 snap', 'pair_coeff * * zbl 22 22',
        f'pair_coeff * * snap {potdir}/Ti_snap.snapcoeff {potdir}/Ti_snap.snapparam Ti'] if kind == 'snap'
       else ['pair_style pod', f'pair_coeff * * {potdir}/Ti_param.pod {potdir}/Ti_coefficients.pod Ti'] if kind == 'pod'
       else ['pair_style meam/spline', f'pair_coeff * * {potdir}/Ti.meam.spline Ti'])
L = lammps(cmdargs=['-log', f'{pre}.log', '-screen', 'none'])
L.commands_list(['units metal', 'atom_style atomic', 'boundary p p p', f'lattice bcc 3.31',
                 f'region b block 0 {nx} 0 {ny} 0 {nz}', 'create_box 1 b', 'create_atoms 1 box', 'mass 1 47.867'] + pot +
                ['timestep 0.002', f'velocity all create {Tg} 4928 mom yes rot yes',
                 'thermo 100', 'thermo_style custom step temp pe press lx ly lz',
                 # 1) equilibrate bcc near Tm
                 f'fix 1 all npt temp {Tg} {Tg} 0.2 aniso 0 0 2.0', 'run 2500', 'unfix 1',
                 # 2) melt top half with bottom half frozen
                 f'variable zmid equal lz/2+zlo', 'region top block INF INF INF INF ${zmid} INF',
                 'group liq region top', 'group sol subtract all liq',
                 'velocity sol set 0 0 0', 'fix f sol setforce 0 0 0',
                 f'fix 2 liq nvt temp 4000 4000 0.2', 'run 2500',
                 f'fix 2 liq nvt temp {Tg} {Tg} 0.2', 'run 1500', 'unfix 2', 'unfix f',
                 f'velocity sol create {Tg} 777 mom yes rot yes',
                 # 3) coexistence: NPH, z free so interface normal can adjust
                 'fix 3 all nph z 0 0 2.0', 'compute pe all pe/atom',
                 'variable T equal temp', f'fix av all ave/time 50 1 50 v_T file {pre}_T.dat', 'run 15000'])
T = np.loadtxt(f'{pre}_T.dat', comments='#')[:, 1]
print(f'Tm estimate (last half mean) = {T[len(T)//2:].mean():.0f} K +/- {T[len(T)//2:].std()/np.sqrt(len(T)/20):.0f}')
L.command(f'write_dump all custom {pre}_final.dump id x y z c_pe')
