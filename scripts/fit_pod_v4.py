import os, glob, shutil, json, sys
import numpy as np
from ase.io import read
E_HCP_REF = -7.8319
COPIES = {'H': 3, 'I': 2, 'B': 4, 'A': 1}
NAMES = {'H': 'Cluster13', 'I': 'Expanded', 'B': 'Burgers', 'A': 'Eos'}
src = sorted(glob.glob('dft_te/done/Te0.1__*.xyz'))
ref = read('dft_te/done/Te0.1__002__A_eos_hcp_v1.00.xyz')
shift = E_HCP_REF - ref.info['E0'] / len(ref)
print(f'energy shift: {shift:.4f} eV/atom')
shutil.rmtree('fit/XYZ_v4', ignore_errors=True); shutil.copytree('fit/XYZ_v3', 'fit/XYZ_v4')
counts = {}
for g, n in COPIES.items():
    sel = [read(f) for f in src if os.path.basename(f).split('__')[2].startswith(g + '_')]
    counts[NAMES[g]] = len(sel); txt = []
    for a in sel:
        e = a.info['E0'] + shift * len(a); fr = a.get_forces()
        lat = ' '.join(f'{x:.8f}' for x in a.cell.array.ravel())
        txt.append(f'{len(a)}\nLattice="{lat}" Properties=species:S:1:pos:R:3:forces:R:3 energy={e:.8f} pbc="T T T"\n')
        txt += [f'Ti {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {q[0]:.6f} {q[1]:.6f} {q[2]:.6f}\n' for p, q in zip(a.positions, fr)]
    open(f'fit/XYZ_v4/{NAMES[g]}.xyz', 'w').write(''.join(txt) * n)
print('added:', counts, 'copies:', COPIES)
out = 'potentials/pod_v4'; os.makedirs(out, exist_ok=True)
shutil.copy('potentials/pod_v3/Ti_param.pod', out)
open(f'{out}/Ti_data.pod', 'w').write('\n'.join([
    'file_format extxyz', 'file_extension xyz', 'path_to_training_data_set "../../fit/XYZ_v4"',
    'path_to_test_data_set "../../fit/test"', 'fitting_weight_energy 1000.0', 'fitting_weight_force 1.0',
    'fitting_regularization_parameter 1e-10', 'error_analysis_for_training_data_set 1',
    'error_analysis_for_test_data_set 1', 'basename_for_output_files Ti', 'precision_for_pod_coefficients 10', '']))
from lammps import lammps
cwd = os.getcwd(); os.chdir(out)
L = lammps(cmdargs=['-log', 'fitpod.log', '-screen', 'none']); L.command('fitpod Ti_param.pod Ti_data.pod'); L.close()
os.chdir(cwd)
print(open(f'{out}/Ti_test_errors.pod').read())
os.makedirs('results', exist_ok=True)
json.dump({'copies': COPIES, 'counts': counts, 'shift': shift}, open('results/fit_v4_settings.json', 'w'), indent=1)
K, rw = 1.0e4, 1.55
r = np.linspace(0.2, 1.6, 1401)
with open(f'{out}/wall.table', 'w') as fw:
    fw.write(f'# repulsive wall {K:g}*(1.55-r)^3 eV below 1.55 A\n\nWALL\nN {len(r)} R {r[0]} {r[-1]}\n\n')
    for i, x in enumerate(r):
        fw.write(f'{i+1} {x:.6f} {K*max(rw-x,0)**3:.8e} {3*K*max(rw-x,0)**2:.8e}\n')
print('wrote potentials/pod_v4 (with wall.table)')
