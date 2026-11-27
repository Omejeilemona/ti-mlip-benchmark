"""Static analysis of the finite-Te DFT data (DFT only). Run from the repository root."""
import glob, os, json
import numpy as np
from ase.io import read
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'stix'})
GPA = 160.21766
D = {}
for f in glob.glob('dft_te/done/Te*__*.xyz'):
    b = os.path.basename(f)[:-4]; kT = float(b.split('__')[0][2:]); D.setdefault(b.split('__', 2)[2], {})[kT] = read(f)
def F(name, kT):
    if kT == 0.0: a = D[name][0.1]; return a.info['E0'] / len(a)
    a = D[name][kT]; return a.info['F'] / len(a)
temps = [0.0] + sorted({k for n in D for k in D[n] if n.startswith('A_eos_hcp')})
os.makedirs('results', exist_ok=True); os.makedirs('figures', exist_ok=True)
out = {'kTe_eV': temps, 'fixed_volume': {}, 'pressure_GPa': {}, 'burgers': {}}
print('F(phase) - F(hcp) (meV/atom) at 0.94/1.00/1.06 of the reference volume; pressure at V0 (GPa)')
print(' kTe     bcc-hcp                  fcc-hcp                  omega-hcp               P_hcp  P_bcc')
for kT in temps:
    row = {ph: [1000 * (F(f'A_eos_{ph}_v{v:.2f}', kT) - F(f'A_eos_hcp_v{v:.2f}', kT)) for v in (0.94, 1.00, 1.06)]
           for ph in ('bcc', 'fcc', 'omega')}
    out['fixed_volume'][kT] = row
    ps = '   -      -'
    if kT > 0:
        P = {p: -np.trace(D[f'A_eos_{p}_v1.00'][kT].get_stress(voigt=False)) / 3 * GPA for p in ('hcp', 'bcc')}
        out['pressure_GPa'][kT] = P; ps = f'{P["hcp"]:6.1f} {P["bcc"]:6.1f}'
    print(f'{kT:4.1f}   ' + '   '.join('/'.join(f'{x:6.1f}' for x in row[p]) for p in ('bcc', 'fcc', 'omega')) + '   ' + ps)
print('\nBurgers path hcp -> bcc, F - F(hcp) (meV/atom)')
fig, ax = plt.subplots(figsize=(3.4, 2.6)); cm = plt.cm.plasma
for j, kT in enumerate(temps):
    y = np.array([F(f'B_burgers_{i:02d}', kT) for i in range(10)]); y = 1000 * (y - y[0]); out['burgers'][kT] = y.tolist()
    print(f'{kT:4.1f}  ' + ' '.join(f'{v:6.1f}' for v in y))
    ax.plot(np.linspace(0, 1, 10), y, 'o-', ms=2.5, lw=1.1, color=cm(j / (len(temps) - 1) * 0.9),
            label='ground state' if kT == 0 else f'{kT:g} eV')
ax.axhline(0, color='0.75', lw=0.6); ax.set_xlabel('Burgers path coordinate (hcp → bcc)'); ax.set_ylabel(r'$F - F_{\mathrm{hcp}}$ (meV/atom)')
ax.legend(frameon=False, fontsize=6.5, title=r'$k_BT_e$', title_fontsize=6.5); fig.tight_layout(); fig.savefig('figures/fig_te_burgers.png', dpi=200)
json.dump(out, open('results/te_phases.json', 'w'), indent=1)
print('\nwrote results/te_phases.json and figures/fig_te_burgers.png')
