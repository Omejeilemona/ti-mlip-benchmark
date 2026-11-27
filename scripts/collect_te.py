"""Collect the finite-electron-temperature DFT campaign into a release dataset. Run from the repository root."""
import glob, os, collections
from ase.io import read, write
files = sorted(glob.glob('dft_te/done/Te*__*.xyz'))
data = collections.defaultdict(dict)
for f in files:
    b = os.path.basename(f)[:-4]; kT = float(b.split('__')[0][2:]); data[b.split('__', 2)[2]][kT] = read(f)
temps = sorted({k for d in data.values() for k in d})
missing = [n for n, d in data.items() if 0.1 not in d]
print(f'{len(files)} calculations, {len(data)} configurations, kTe = {temps} eV')
print('configurations without a 0.1 eV reference:', missing or 'none')
os.makedirs('data_te', exist_ok=True)
with open('data_te/summary.csv', 'w') as s:
    s.write('config,kTe_eV,natoms,F_eV,E0_eV,F_minus_E0ref_eV_per_atom\n')
    for kT in temps:
        frames = []
        for n in sorted(data):
            if kT not in data[n]: continue
            a = data[n][kT]; frames.append(a)
            ref = data[n].get(0.1)
            d = (a.info['F'] - ref.info['E0']) / len(a) if ref is not None else float('nan')
            s.write(f'{n},{kT},{len(a)},{a.info["F"]:.6f},{a.info["E0"]:.6f},{d:.6f}\n')
        write(f'data_te/Ti_Te{kT:.1f}.extxyz', frames, format='extxyz')
        print(f'  kTe = {kT:.1f} eV: {len(frames)} configurations -> data_te/Ti_Te{kT:.1f}.extxyz')
