"""Add the EOS cells and the Burgers path at k_BTe = 0.2, 0.3 and 0.4 eV to the queue (hcp-bcc crossover).
Run from dft_te/:  python3 ../scripts/gen_extra_te.py"""
import glob, os
from ase.io import read, write
n = 0
for te in (0.2, 0.3, 0.4):
    for f in sorted(glob.glob('done/Te0.1__*__A_eos_*.xyz') + glob.glob('done/Te0.1__*__B_burgers_*.xyz')):
        dst = 'queue/' + os.path.basename(f).replace('Te0.1__', f'Te{te:.1f}__')
        if not os.path.exists(dst) and not os.path.exists(dst.replace('queue/', 'done/')):
            a = read(f); a.calc = None
            for k in list(a.info): del a.info[k]
            write(dst, a, format='extxyz'); n += 1
print(f'added {n} jobs to queue/')
