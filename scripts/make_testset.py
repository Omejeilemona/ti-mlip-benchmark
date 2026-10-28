"""Write the 300 held-out test configurations (150 Init + 150 Bulk) to extended XYZ for POD testing.
Uses the same random selection (seed 1) as heldout.py."""
import numpy as np, os
D = 'data/dft/261022_AllInitAndBulk_256atoms_with_types_curatedData/'
H = np.load('fit/holdout_ids.npz'); rng = np.random.default_rng(1); os.makedirs('fit/test', exist_ok=True)
E = np.load(D+'energy.npy'); F = np.load(D+'force.npy').reshape(len(E), -1, 3)
C = np.load(D+'coord.npy').reshape(len(E), -1, 3); B = np.load(D+'box.npy').reshape(-1, 3, 3)
for g in ('Init', 'Bulk'):
    ids = rng.choice(H[g], size=150, replace=False)
    with open(f'fit/test/{g}.xyz', 'w') as f:
        for i in ids:
            lat = ' '.join(f'{x:.8f}' for x in B[i].ravel())
            f.write(f'{C.shape[1]}\nLattice="{lat}" Properties=species:S:1:pos:R:3:forces:R:3 energy={E[i]:.8f} pbc="T T T"\n')
            for r, fr in zip(C[i], F[i]):
                f.write(f'Ti {r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {fr[0]:.6f} {fr[1]:.6f} {fr[2]:.6f}\n')
