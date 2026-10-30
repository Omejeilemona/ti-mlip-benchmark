"""Held-out errors of a MACE model on the same 300 configurations used for SNAP and POD.
usage: python scripts/mace_eval.py <model> [out_prefix]   (run from the repository root)"""
import sys, numpy as np
from ase.io import read
from mace.calculators import MACECalculator
model = sys.argv[1]; pre = sys.argv[2] if len(sys.argv) > 2 else 'results/heldout_mace'
calc = MACECalculator(model_paths=model, device='cpu', default_dtype='float64')
for g in ('Init', 'Bulk'):
    ed, ep, fd, fp = [], [], [], []
    for at in read(f'mace/test_{g}.xyz', ':'):
        n = len(at); ed.append(at.info['REF_energy'] / n); fd.append(at.arrays['REF_forces'].copy())
        at.calc = calc; ep.append(at.get_potential_energy() / n); fp.append(at.get_forces())
    ed, ep = np.array(ed), np.array(ep); fd, fp = np.concatenate(fd), np.concatenate(fp)
    r = (fp - fd).ravel()
    print(g, dict(n=len(ed), E_rmse_meV=round(1000 * np.sqrt(np.mean((ep - ed) ** 2)), 2),
                  F_rmse=round(np.sqrt(np.mean(r ** 2)), 4), F_mae=round(np.mean(np.abs(r)), 4),
                  F_std_dft=round(fd.std(), 4)), flush=True)
    np.savez(f'{pre}_{g}.npz', ed=ed, ep=ep, fd=fd, fp=fp)
