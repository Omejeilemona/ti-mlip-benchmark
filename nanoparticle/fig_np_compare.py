"""Nanoparticle figure: (a) snapshots EAM vs EAM+dF(Te); (b) largest fragment vs absorbed energy; (c) Te, Tl.
usage: python3 fig_np_compare.py [E_abs for snapshots, default 3.5]"""
import sys, os, numpy as np
def rd(f): return np.genfromtxt(f, delimiter=',', names=True)
def smooth(y, n=5): k = np.ones(n); return np.convolve(y, k, mode='same') / np.convolve(np.ones_like(y), k, mode='same')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from atomrender import render
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'stix'})
ESNAP = float(sys.argv[1]) if len(sys.argv) > 1 else 3.5
POTS = [('eam', 'EAM (Ti1)'), ('eamte', r'EAM + $\Delta F(T_e)$')]
COL = {'eam': '#3b6ea8', 'eamte': '#c0392b'}
ENERGIES = [1.5, 2, 2.5, 3, 3.5, 4, 5]
TIMES = [0, 2, 4, 8, 15]
C_BIG, C_FRAG, C_MONO = (0.40, 0.52, 0.68), (0.93, 0.55, 0.18), (0.80, 0.15, 0.15)
run = lambda p, e: f'{p}_d2_e{e:g}'
def colours(x):
    pr = cKDTree(x).query_pairs(3.6, output_type='ndarray'); n = len(x)
    nc, lab = connected_components(coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n)), directed=False)
    sz = np.bincount(lab); big = sz.argmax(); s = sz[lab]
    c = np.tile(C_FRAG, (n, 1)); c[lab == big] = C_BIG; c[s == 1] = C_MONO
    return c, x[lab == big].mean(0)
fig = plt.figure(figsize=(7.0, 5.0))
gs = fig.add_gridspec(4, len(TIMES), height_ratios=[1, 1, 0.15, 1.3], hspace=0.08, wspace=0.03)
for r, (p, lab) in enumerate(POTS):
    for c, t in enumerate(TIMES):
        ax = fig.add_subplot(gs[r, c]); ax.set_xticks([]); ax.set_yticks([])
        f = f'{run(p, ESNAP)}/snap_{t:06.2f}ps.dump'
        if os.path.exists(f):
            x = np.loadtxt(f, skiprows=9)[:, 1:4]; col, com = colours(x)
            ax.imshow(render(x - com, col, fov=40, size=420))
        for s in ax.spines.values(): s.set_color('0.8')
        if r == 0: ax.set_title(f'{t} ps', fontsize=8)
        if c == 0: ax.set_ylabel(lab, fontsize=7.5)
fig.text(0.07, 0.955, f'(a)  {ESNAP:g} eV/atom absorbed', fontsize=8.5, fontweight='bold')
ax0 = fig.axes[len(TIMES) - 1]
ax0.plot([0.70, 0.95], [0.06, 0.06], color='k', lw=1.2, transform=ax0.transAxes)
ax0.text(0.825, 0.10, '2 nm', transform=ax0.transAxes, ha='center', fontsize=7)
fig.legend([Line2D([], [], marker='o', ls='', color=cc, ms=6) for cc in (C_BIG, C_FRAG, C_MONO)],
           ['largest fragment', 'fragments', 'monomers'], loc='upper center', ncol=3, frameon=False,
           bbox_to_anchor=(0.62, 1.0), fontsize=7.5)
sub = gs[3, :].subgridspec(1, 2, wspace=0.3)
axb, axc = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
for p, lab in POTS:
    E, frac = [], []
    for e in ENERGIES:
        f = f'{run(p, e)}/ttm.csv'
        if os.path.exists(f):
            d = rd(f); E.append(e); frac.append(d['largest'][-1] / d['largest'][0])
    axb.plot(E, frac, 'o-', ms=3.5, color=COL[p], label=lab)
    f = f'{run(p, ESNAP)}/ttm.csv'
    if os.path.exists(f):
        d = rd(f)
        axc.plot(d['t_ps'], d['Te_K'] / 1e3, '-', color=COL[p], lw=1.2)
        axc.plot(d['t_ps'], smooth(d['Tl_K']) / 1e3, '--', color=COL[p], lw=1.0)
axb.set_xlabel('Absorbed energy (eV/atom)'); axb.set_title('after 30 ps', fontsize=7.5); axb.set_ylabel('Largest fragment / $N$'); axb.set_ylim(0, 1.05)
axb.legend(frameon=False, fontsize=7); axb.text(-0.24, 1.06, '(b)', transform=axb.transAxes, fontsize=9, fontweight='bold')
axc.plot([], [], 'k-', lw=1.2, label='$T_e$'); axc.plot([], [], 'k--', lw=1.0, label='$T_l$'); axc.legend(frameon=False, fontsize=7)
axc.set_xlim(0, 30); axc.set_xlabel('Time (ps)'); axc.set_ylabel('Temperature (10$^3$ K)')
axc.text(-0.20, 1.06, '(c)', transform=axc.transAxes, fontsize=9, fontweight='bold')
fig.savefig('fig_np_compare.pdf', bbox_inches='tight', dpi=300); fig.savefig('fig_np_compare.png', bbox_inches='tight', dpi=200)
print('wrote fig_np_compare.pdf / .png')
