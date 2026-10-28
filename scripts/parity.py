import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
models=[('meam','MEAM-spline (Hennig 2008)','#888888'),('snap','SNAP (linear, 2J=10)','#1f77b4'),('pod','POD (2-, 3-, 4-body)','#d62728')]
fig,ax=plt.subplots(2,3,figsize=(10,6.4))
for j,(m,lab,c) in enumerate(models):
    for g,mk in (('Init','o'),('Bulk','s')):
        d=np.load(f'heldout_{m}_{g}.npz')
        ax[0,j].scatter(d['ed'],d['ep'],s=8,marker=mk,c=c,alpha=.6,label=g)
        idx=np.random.default_rng(0).choice(d['fd'].size,6000,replace=False)
        ax[1,j].scatter(d['fd'].ravel()[idx],d['fp'].ravel()[idx],s=2,c=c,alpha=.3)
    for a,lim in ((ax[0,j],(-7.9,-6.5)),(ax[1,j],(-6,6))):
        a.plot(lim,lim,'k-',lw=.7); a.set_xlim(lim); a.set_ylim(lim); a.set_aspect('equal')
    ax[0,j].set_title(lab); ax[0,j].set_xlabel('DFT energy (eV/atom)'); ax[0,j].set_ylabel('Model energy (eV/atom)')
    ax[1,j].set_xlabel('DFT force (eV/Å)'); ax[1,j].set_ylabel('Model force (eV/Å)')
    ax[0,j].legend(fontsize=7,loc='upper left')
fig.suptitle('Ti: held-out DFT test set (300 configurations, 76,800 atoms, never used in fitting)\nMEAM energies shifted by a constant (different energy zero)',fontsize=10)
fig.tight_layout(); fig.savefig('parity.png',dpi=160); print('saved')
