"""Convert Röcken & Zavadlav (2024) Ti DFT data (DeePMD-style npy) to FitSNAP extended-XYZ.
Units in: eV, Angstrom, eV/A, virial in eV (LAMMPS sign: positive = compressive).
FitSNAP xyz 'stress' = LAMMPS pressure tensor in bar = virial/V * 1.602176634e6.
"""
import numpy as np, sys, os
src, out = sys.argv[1], sys.argv[2]
rng = np.random.default_rng(0)
L = lambda n: np.load(os.path.join(src, n + '.npy'))
E, F, C, B, W, T = L('energy'), L('force'), L('coord'), L('box'), L('virial'), L('types')
B = B.reshape(-1, 3, 3); C = C.reshape(len(E), -1, 3); F = F.reshape(len(E), -1, 3); W = W.reshape(-1, 3, 3)
EV_A3_TO_BAR = 1.602176634e6
idx = {'Init': np.where(T == 0)[0], 'Bulk': np.where(T == 1)[0]}
nfit = {'Init': int(sys.argv[3]), 'Bulk': int(sys.argv[4])}
holdout = {}
for g, ids in idx.items():
    ids = rng.permutation(ids)
    use, rest = ids[:nfit[g]], ids[nfit[g]:]
    holdout[g] = rest
    with open(os.path.join(out, f'{g}.xyz'), 'w') as f:
        for i in use:
            V = abs(np.linalg.det(B[i])); P = W[i] / V * EV_A3_TO_BAR
            lat = ' '.join(f'{x:.8f}' for x in B[i].ravel())
            st = ' '.join(f'{x:.6f}' for x in P.ravel())
            f.write(f'{C.shape[1]}\nLattice="{lat}" Properties=species:S:1:pos:R:3:forces:R:3 energy={E[i]:.8f} stress="{st}" pbc="T T T"\n')
            for r, fr in zip(C[i], F[i]):
                f.write(f'Ti {r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {fr[0]:.6f} {fr[1]:.6f} {fr[2]:.6f}\n')
np.savez(os.path.join(out, '..', 'holdout_ids.npz'), **holdout)
print({g: (nfit[g], len(holdout[g])) for g in idx})
