"""One DFT worker: repeatedly claims the next configuration in queue/, runs GPAW, writes done/<name>.
Safe to run several in parallel and to stop at any time: an interrupted job is moved back with
`python scripts/dft_worker.py --requeue`. Settings are identical for every configuration."""
import os, sys, time, glob, socket, traceback
from ase.io import read, write
SETTINGS = dict(pw=450, kdensity=5.0, smearing=0.1)   # eV, points per 1/A, eV (Fermi-Dirac)
for d in ('queue', 'running', 'done', 'failed', 'logs'): os.makedirs(d, exist_ok=True)
if '--requeue' in sys.argv:
    for f in glob.glob('running/*'):
        os.rename(f, 'queue/' + os.path.basename(f).split('__', 1)[1]); print('requeued', f)
    sys.exit()
from gpaw import GPAW, PW, FermiDirac
tag = f'{socket.gethostname()}-{os.getpid()}'
while True:
    claimed = None
    for f in sorted(glob.glob('queue/*.xyz')):
        dst = f'running/{tag}__{os.path.basename(f)}'
        try:
            os.rename(f, dst); claimed = dst; break          # atomic: only one worker wins
        except FileNotFoundError:
            continue
    if claimed is None:
        print(tag, 'queue empty, stopping', flush=True); break
    name = claimed.split('__', 1)[1]; t0 = time.time()
    try:
        at = read(claimed)
        at.calc = GPAW(mode=PW(SETTINGS['pw']), xc='PBE', kpts={'density': SETTINGS['kdensity'], 'gamma': True},
                       occupations=FermiDirac(SETTINGS['smearing']), txt=f'logs/{name}.txt')
        e = at.get_potential_energy(); f = at.get_forces(); s = at.get_stress(voigt=False)
        at.info.update(config=name, seconds=round(time.time() - t0), **{k: v for k, v in SETTINGS.items()})
        write(f'done/{name}', at, format='extxyz')
        os.remove(claimed)
        print(f'{time.strftime("%H:%M")} {name}: {len(at)} atoms, E/atom {e/len(at):.4f} eV, {time.time()-t0:.0f} s', flush=True)
    except Exception:
        traceback.print_exc(); os.rename(claimed, f'failed/{name}')
