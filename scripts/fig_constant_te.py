"""Constant-Te MD figure: (a) phase fractions vs time; (b) final structures coloured by Ackland-Jones analysis."""
import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, 'nanoparticle'); from atomrender import render
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'stix'})
RUNS = [('kte0.5_npt_res', '0.5 eV, flexible cell'), ('kte0.5_nvt_res', '0.5 eV, fixed cell'),
        ('kte1.0_nvt_res', '1.0 eV, fixed cell'), ('kte0.5_npt', '0.5 eV, flexible cell,\nfirst correction')]
CL = ['#222222', '#e08a1e', '#7d5ba6', '0.65']
COLS = {3: (0.85, 0.25, 0.25), 2: (0.25, 0.65, 0.35), 1: (0.25, 0.45, 0.80), 0: (0.75, 0.75, 0.75), 4: (0.9, 0.8, 0.2)}
fig = plt.figure(figsize=(7.0, 4.9))
gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.15], hspace=0.32, wspace=0.08)
axa = fig.add_subplot(gs[0, :3])
for (tag, lab), c in zip(RUNS, CL):
    f = f'md_te/{tag}.csv'
    if not os.path.exists(f): continue
    d = np.loadtxt(f, delimiter=',', skiprows=1)
    axa.plot(d[:, 0], d[:, 1], '-', color=c, lw=1.2); axa.plot(d[:, 0], d[:, 3], '--', color=c, lw=1.0)
    axa.plot(d[:, 0], d[:, 2], ':', color=c, lw=1.3)
axa.set_xlabel('Time after switching on $T_e$ (ps)'); axa.set_ylabel('Fraction of atoms'); axa.set_xlim(0, 20); axa.set_ylim(-0.02, 1.02)
axa.text(-0.08, 1.05, '(a)', transform=axa.transAxes, fontsize=9, fontweight='bold')
axl = fig.add_subplot(gs[0, 3]); axl.axis('off')
h1 = [Line2D([], [], color=c, lw=1.6) for c in CL]
h2 = [Line2D([], [], color='k', ls=s, lw=1.2) for s in ('-', '--', ':')]
lg = axl.legend(h1, [l for _, l in RUNS], frameon=False, fontsize=7, loc='upper left', bbox_to_anchor=(0.02, 1.02))
axl.add_artist(lg)
axl.legend(h2, ['hcp', 'fcc', 'bcc'], frameon=False, fontsize=7, loc='lower left', bbox_to_anchor=(0.02, 0.0), title='phase', title_fontsize=7)
for i, (tag, lab) in enumerate(RUNS):
    ax = fig.add_subplot(gs[1, i]); ax.set_xticks([]); ax.set_yticks([])
    f = f'md_te/{tag}_final.dump'
    if os.path.exists(f):
        d = np.loadtxt(f, skiprows=9); x = d[:, 1:4]; ack = d[:, 4].astype(int)
        col = np.array([COLS.get(a, COLS[0]) for a in ack])
        ax.imshow(render(x - x.mean(0), col, radius=1.2, fov=13, size=420, view=(-1.0, -0.12, -0.18)))
    ax.set_title(lab, fontsize=7, color=CL[i] if i < 3 else '0.4')
    for s in ax.spines.values(): s.set_color('0.8')
fig.legend([Line2D([], [], marker='o', ls='', color=COLS[k], ms=6) for k in (3, 2, 1, 0)], ['hcp', 'fcc', 'bcc', 'other'],
           loc='lower center', ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.07))
fig.text(0.075, 0.49, '(b)', fontsize=9, fontweight='bold')
fig.savefig('figures/fig_constant_te.pdf', bbox_inches='tight'); fig.savefig('figures/fig_constant_te.png', dpi=200, bbox_inches='tight')
print('wrote figures/fig_constant_te.pdf/.png')
