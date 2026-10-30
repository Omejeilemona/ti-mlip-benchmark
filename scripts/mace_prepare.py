"""Write MACE training, validation and test files from the same split used for SNAP and POD.
usage: python scripts/mace_prepare.py [n_init n_bulk]   (default 150 250)   -> mace/
Keys are renamed to REF_energy / REF_forces so ASE never mistakes them for a model prediction."""
import sys, os, re, random
ni, nb = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (150, 250)
def frames(f):
    L = open(f).read().splitlines(); i = 0; out = []
    while i < len(L):
        n = int(L[i]); h = L[i + 1].replace(' energy=', ' REF_energy=').replace('forces:R:3', 'REF_forces:R:3')
        h = re.sub(r' stress="[^"]*"', '', h)
        out.append('\n'.join([L[i], h] + L[i + 2:i + 2 + n])); i += 2 + n
    return out
I, B = frames('fit/XYZ/Init.xyz'), frames('fit/XYZ/Bulk.xyz')   # the POD/SNAP training pool only
random.seed(0); random.shuffle(I); random.shuffle(B)
os.makedirs('mace', exist_ok=True)
open('mace/train.xyz', 'w').write('\n'.join(I[:ni] + B[:nb]) + '\n')
open('mace/valid.xyz', 'w').write('\n'.join(I[ni:ni + 20] + B[nb:nb + 30]) + '\n')
for g in ('Init', 'Bulk'):                                      # the 300 held-out configurations
    open(f'mace/test_{g}.xyz', 'w').write('\n'.join(frames(f'fit/test/{g}.xyz')) + '\n')
print(f'train {ni + nb}, valid 50, test 300 -> mace/')
