"""Write MACE training/validation/test files from the full v3 training pool (all 8 categories) and the
300 held-out test configurations used for SNAP and POD, and pack them into mace_v3.zip for Colab.
Run from the repository root:  python3 scripts/mace_prepare_v3.py"""
import os, re, random, zipfile
def frames(f):
    L = open(f).read().splitlines(); i = 0; out = []
    while i < len(L):
        n = int(L[i]); h = L[i + 1].replace(' energy=', ' REF_energy=').replace('forces:R:3', 'REF_forces:R:3')
        h = re.sub(r' stress="[^"]*"', '', h)
        out.append('\n'.join([L[i], h] + L[i + 2:i + 2 + n])); i += 2 + n
    return out
CATS = ('Init', 'Bulk', 'Omega', 'Path', 'Defect', 'Surface', 'Hotsolid', 'Liquid')
random.seed(0); train, valid = [], []
for c in CATS:
    fr = frames(f'fit/XYZ_v3/{c}.xyz'); random.shuffle(fr)
    nv = max(2, round(0.05 * len(fr)))                       # 5% of every category for validation
    valid += fr[:nv]; train += fr[nv:]
    print(f'{c:9s} {len(fr) - nv:4d} train  {nv:3d} valid')
os.makedirs('mace_v3', exist_ok=True)
open('mace_v3/train.xyz', 'w').write('\n'.join(train) + '\n')
open('mace_v3/valid.xyz', 'w').write('\n'.join(valid) + '\n')
for g in ('Init', 'Bulk'):
    open(f'mace_v3/test_{g}.xyz', 'w').write('\n'.join(frames(f'fit/test/{g}.xyz')) + '\n')
with zipfile.ZipFile('mace_v3.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for f in ('train.xyz', 'valid.xyz', 'test_Init.xyz', 'test_Bulk.xyz'): z.write(f'mace_v3/{f}', f)
print(f'train {len(train)}, valid {len(valid)}, test 300 -> mace_v3.zip ({os.path.getsize("mace_v3.zip") / 1e6:.1f} MB)')
