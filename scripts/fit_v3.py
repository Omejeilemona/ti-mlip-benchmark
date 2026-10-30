"""Build the v3 training set from finished DFT jobs (done/) and fit POD v3.
Run from the repository root after get_data.sh, convert.py and make_testset.py.
usage: python scripts/fit_v3.py [category=copies ...]   e.g.  omega=6 vacancy=4
Steps: align GPAW energies to the reference data, check consistency, write fit/XYZ_v3, run fitpod,
then compute properties and held-out errors for potentials/pod_v3."""
import sys, os, glob, shutil, json, subprocess
import numpy as np
from ase.io import read
E_HCP_REF = -7.8319                      # eV/atom, hcp energy of POD v1 (fitted to the VASP-PBE data)
REF_DIFF = {'bcc': 105.3, 'fcc': 52.4}    # meV/atom, POD v1 phase energies (Init data well covered)
copies = dict(omega=4, defect=3, liquid=2, hotsolid=2, surface=3, path=3)
for a in sys.argv[1:]:
    k, v = a.split('='); copies[k] = int(v)
done = sorted(glob.glob(os.environ.get('DONE_DIR', 'dft_v3/done') + '/*.xyz'))
frames = {f: read(f) for f in done}
hcp = [a for f, a in frames.items() if os.path.basename(f) == '00_align_001.xyz']
if not hcp: sys.exit('the hcp alignment job (00_align_001) has not finished yet')
e_hcp = hcp[0].get_potential_energy() / len(hcp[0]); shift = E_HCP_REF - e_hcp
print(f'energy shift GPAW -> reference: {shift:.4f} eV/atom')
for f, a in frames.items():
    nm = os.path.basename(f)
    for ph, i in (('bcc', '002'), ('fcc', '003'), ('omega', '004')):
        if nm == f'00_align_{i}.xyz':
            d = (a.get_potential_energy() / len(a) - e_hcp) * 1000
            ref = REF_DIFF.get(ph)
            print(f'  {ph:5s} - hcp: GPAW {d:7.1f} meV/atom' + (f'   reference fit {ref:.1f}' if ref else '   (no reference: not in original data)'))
os.makedirs('fit/XYZ_v3', exist_ok=True)
for g in ('Init', 'Bulk'): shutil.copy(f'fit/XYZ/{g}.xyz', f'fit/XYZ_v3/{g}.xyz')
counts = {}
for cat, n in copies.items():
    sel = [a for f, a in frames.items() if f'_{cat}_' in os.path.basename(f)]
    if not sel: continue
    counts[cat] = len(sel); txt = []
    for a in sel:
        e = a.get_potential_energy() + shift * len(a); fr = a.get_forces()
        lat = ' '.join(f'{x:.8f}' for x in a.cell.array.ravel())
        txt.append(f'{len(a)}\nLattice="{lat}" Properties=species:S:1:pos:R:3:forces:R:3 energy={e:.8f} pbc="T T T"\n')
        txt += [f'Ti {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {g[0]:.6f} {g[1]:.6f} {g[2]:.6f}\n' for p, g in zip(a.positions, fr)]
    open(f'fit/XYZ_v3/{cat.capitalize()}.xyz', 'w').write(''.join(txt) * n)
print('added configurations:', counts, ' copies:', {k: copies[k] for k in counts})
out = 'potentials/pod_v3'; os.makedirs(out, exist_ok=True)
shutil.copy('potentials/pod_v1/Ti_param.pod', out)
open(f'{out}/Ti_data.pod', 'w').write('\n'.join([
    'file_format extxyz', 'file_extension xyz', 'path_to_training_data_set "../../fit/XYZ_v3"',
    'path_to_test_data_set "../../fit/test"', 'fitting_weight_energy 1000.0', 'fitting_weight_force 1.0',
    'fitting_regularization_parameter 1e-10', 'error_analysis_for_training_data_set 1',
    'error_analysis_for_test_data_set 1', 'basename_for_output_files Ti', 'precision_for_pod_coefficients 10', '']))
from lammps import lammps
cwd = os.getcwd(); os.chdir(out)
L = lammps(cmdargs=['-log', 'fitpod.log', '-screen', 'none']); L.command('fitpod Ti_param.pod Ti_data.pod'); L.close()
os.chdir(cwd)
print(open(f'{out}/Ti_test_errors.pod').read())
props = subprocess.run([sys.executable, 'scripts/props.py', 'pod', out], capture_output=True, text=True).stdout
open('results/props_pod_v3.json', 'w').write(props); print(props)
print(subprocess.run([sys.executable, 'scripts/heldout.py', 'pod', out, '150', 'results/heldout_pod_v3.npz'], capture_output=True, text=True).stdout)
json.dump({'copies': copies, 'counts': counts, 'shift_eV_per_atom': shift}, open('results/fit_v3_settings.json', 'w'), indent=1)
