"""Electron heat capacity of hcp Ti from F(Te) (ce_dft.py): S_e = -dF/dT, E_e = F + T S_e, c_e = dE_e/dT.
Compares with Lin & Zhigilei (nanoparticle/Ce_Ti_LinZhigilei.dat). Run from the repository root:
    python3 scripts/ce_analysis.py  -> figures/fig_ce.pdf/.png and the numbers for the paper"""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'stix'})
kB, NA, eV = 8.617333e-5, 6.02214076e23, 1.602176634e-19
n0 = 4.506e3 / (47.867 * 1.66054e-27)                        # atoms per m^3
kT, F = np.loadtxt('ce_hcp_F_vs_Te.txt', unpack=True)
T = kT / kB
S = -np.gradient(F, T); E = F + T * S; ce = np.gradient(E, T)  # eV/K/atom
kT, T, ce = kT[2:-2], T[2:-2], ce[2:-2]                        # drop end points (one-sided differences)
print(f'{len(T)} points used, kTe = {kT.min():.2f}-{kT.max():.2f} eV')
lT, lc = np.loadtxt('nanoparticle/Ce_Ti_LinZhigilei.dat', skiprows=1, unpack=True)  # K, eV/K/A^3
lc_atom = lc / (n0 * 1e-30)
m = (T >= lT.min()) & (T <= lT.max())
dev = ce[m] / np.interp(T[m], lT, lc_atom) - 1
print(f'deviation from Lin & Zhigilei over {T[m].min():.0f}-{T[m].max():.0f} K: '
      f'max |dev| = {100 * np.abs(dev).max():.1f}%, mean = {100 * dev.mean():+.1f}%')
for t in (5000, 10000, 15000, 20000, 25000):
    if T.min() <= t <= T.max():
        print(f'  Te = {t:5d} K: c_e = {np.interp(t, T, ce) * eV * n0 / 1e5:5.2f} x 1e5 J/m3/K   '
              f'(Lin & Zhigilei {np.interp(t, lT, lc_atom) * eV * n0 / 1e5:5.2f})')
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.plot(lT / 1e3, lc_atom * eV * n0 / 1e5, '-', color='0.55', lw=1.5, label='Lin & Zhigilei')
ax.plot(T / 1e3, ce * eV * n0 / 1e5, 'o', color='#c0392b', ms=3, label='This work (PBE, hcp)')
ax.set_xlabel('$T_e$ (10$^3$ K)'); ax.set_ylabel('$C_e$ (10$^5$ J m$^{-3}$ K$^{-1}$)')
ax.set_xlim(0, T.max() / 1e3 * 1.02); ax.set_ylim(0, None); ax.legend(frameon=False, fontsize=7)
fig.savefig('figures/fig_ce.pdf', bbox_inches='tight'); fig.savefig('figures/fig_ce.png', dpi=200, bbox_inches='tight')
print('wrote figures/fig_ce.pdf/.png')
