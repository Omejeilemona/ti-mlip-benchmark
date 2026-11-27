"""DOS of hcp, bcc, fcc and omega Ti near the Fermi level (results/dos_*.dat) -> figures/fig_dos.pdf/.png"""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'stix'})
C = {'hcp': 'k', 'bcc': '#c0392b', 'fcc': '#2e8b57', 'omega': '#7d5ba6'}; L = {'omega': r'$\omega$'}
fig, ax = plt.subplots(figsize=(3.4, 2.6))
for ph in ('hcp', 'bcc', 'fcc', 'omega'):
    E, g = np.loadtxt(f'results/dos_{ph}.dat', unpack=True)
    ax.plot(E, g, color=C[ph], lw=1.1, label=L.get(ph, ph))
for t, a in ((0.3, 0.18), (1.0, 0.08)):
    ax.axvspan(-2 * t, 2 * t, color='0.5', alpha=a, lw=0)
ax.text(0.05, 0.93, r'$\pm 2k_BT_e$: 0.3 and 1 eV', transform=ax.transAxes, fontsize=6.8, color='0.35')
ax.axvline(0, color='0.3', lw=0.6, ls=':')
ax.set_xlim(-6, 4); ax.set_ylim(0, None); ax.set_xlabel(r'$E - E_F$ (eV)'); ax.set_ylabel('DOS (states eV$^{-1}$ atom$^{-1}$)')
ax.legend(frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(0.0, 0.9))
fig.savefig('figures/fig_dos.pdf', bbox_inches='tight'); fig.savefig('figures/fig_dos.png', dpi=200, bbox_inches='tight')
print('wrote figures/fig_dos.pdf/.png')
