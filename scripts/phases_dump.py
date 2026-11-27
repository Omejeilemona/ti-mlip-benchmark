"""Crystal structure of the final coexistence snapshots after a quench (energy minimisation with the potential),
which removes thermal noise; then Ackland-Jones. usage: python3 ../scripts/phases_dump.py <potdir> run_*_final.dump"""
import sys, numpy as np
from lammps import lammps
pot = sys.argv[1]
for f in sys.argv[2:]:
    L = open(f).read().split('\n'); n = int(L[3]); bx = [list(map(float, L[i].split()[:2])) for i in (5, 6, 7)]
    d = np.loadtxt(f, skiprows=9)[:, 1:4]
    lmp = lammps(cmdargs=['-log', 'none', '-screen', 'none'])
    lmp.commands_list(['units metal', 'atom_style atomic', 'boundary p p p',
        f'region b block {bx[0][0]} {bx[0][1]} {bx[1][0]} {bx[1][1]} {bx[2][0]} {bx[2][1]}', 'create_box 1 b', 'mass 1 47.867'])
    lmp.create_atoms(n, None, [1] * n, d.ravel().tolist())
    lmp.commands_list(['pair_style pod', f'pair_coeff * * {pot}/Ti_param.pod {pot}/Ti_coefficients.pod Ti',
                       'compute a all ackland/atom', 'min_style fire', 'minimize 0 1e-3 2000 20000', 'run 0'])
    a = lmp.numpy.extract_compute('a', 1, 1)[:n].astype(int); lmp.close()
    fr = {k: np.mean(a == v) for k, v in (('bcc', 1), ('fcc', 2), ('hcp', 3), ('other', 0))}
    print(f'{f} (quenched): ' + '  '.join(f'{k} {v:.2f}' for k, v in fr.items()), flush=True)
