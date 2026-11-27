"""Finite-electron-temperature DFT worker (GPAW, PBE).
Claims the next job in queue/ (file name Te<kTe>__<natoms>__<config>.xyz), runs it, writes done/<job>.xyz.
Stored per job: free energy F = E - Te*S (energy='F', the quantity the forces derive from), the
zero-width extrapolated energy E0, forces, stress, kTe, number of bands and occupation of the highest band.
Settings for every Te and every configuration (so that dF(Te) = F(Te) - F(0.1 eV) is consistent):
    plane waves 450 eV, k-point density 3.5 per 1/A (Gamma-centred; 1 k-point for clusters in vacuum),
    Fermi-Dirac occupations with width kTe, bands 200% (kTe <= 1 eV) or 300% (kTe >= 1.5 eV);
    clusters in vacuum: 1.2 x occupied + 16/60/140 bands at kTe = 0.1/0.5/1 eV (vacuum states), and only
    kTe <= 1 eV (at higher Te a box-size-dependent fraction of electrons is thermally emitted into the vacuum).
Large jobs (>= 16 atoms at kTe >= 1 eV) take a lock so that only one runs at a time (memory).
Interrupted jobs are put back with:  python3 scripts/dft_worker_te.py --requeue"""
import os, sys, time, glob, socket, traceback
import numpy as np
from ase.io import read, write
for d in ('queue', 'running', 'done', 'failed', 'logs'): os.makedirs(d, exist_ok=True)
LOCK = 'running/BIG.lock'
if '--requeue' in sys.argv:
    for f in glob.glob('running/*__*.xyz'):
        os.rename(f, 'queue/' + os.path.basename(f).split('__', 1)[1]); print('requeued', f)
    if os.path.exists(LOCK): os.remove(LOCK)
    sys.exit()
from gpaw import GPAW, PW, FermiDirac
tag = f'{socket.gethostname()}-{os.getpid()}'
def is_big(name):
    kT, n = float(name.split('__')[0][2:]), int(name.split('__')[1]); return False
def claim():
    for f in sorted(glob.glob('queue/*.xyz')):
        name = os.path.basename(f); big = is_big(name)
        if big:
            try: os.close(os.open(LOCK, os.O_CREAT | os.O_EXCL))   # only one big job at a time
            except FileExistsError: continue
        try:
            os.rename(f, f'running/{tag}__{name}'); return name, big
        except FileNotFoundError:
            if big: os.remove(LOCK)
    return None, False
while True:
    name, big = claim()
    if name is None:
        if glob.glob('queue/*.xyz'):            # only big jobs left and one is running: wait
            time.sleep(120); continue
        print(tag, 'queue empty, stopping', flush=True); break
    path = f'running/{tag}__{name}'; t0 = time.time(); kT = float(name.split('__')[0][2:])
    try:
        at = read(path)
        kp = (1, 1, 1) if min(at.cell.lengths()) > 12 else {'density': 3.5, 'gamma': True}
        if kp == (1, 1, 1):      # cluster in vacuum: free-electron-like vacuum states need explicit extra bands
            nb = int(1.2 * 6 * len(at)) + int(round({0.1: 16, 0.5: 60, 1.0: 140}[kT] * max(1.0, at.get_volume() / 1953.0)))
        else:
            nb = '200%' if kT <= 1.0 else '300%'
        # clusters in vacuum: many nearly empty smeared bands converge slowly; 1e-6 eV^2 per electron is ample
        conv = {'eigenstates': 1e-6} if kp == (1, 1, 1) else {}
        calc = GPAW(mode=PW(450), xc='PBE', kpts=kp, occupations=FermiDirac(kT), nbands=nb, convergence=conv,
                    txt=f'logs/{name[:-4]}.txt')
        at.calc = calc
        F = at.get_potential_energy(force_consistent=True); E0 = at.get_potential_energy()
        f = at.get_forces(); s = at.get_stress()
        occ = max(calc.get_occupation_numbers(k, sp)[-1] for k in range(len(calc.get_ibz_k_points()))
                  for sp in range(calc.get_number_of_spins())) / (2 / calc.get_number_of_spins())
        out = at.copy()
        from ase.calculators.singlepoint import SinglePointCalculator
        out.calc = SinglePointCalculator(out, energy=F, free_energy=F, forces=f, stress=s)
        out.info.update(config=name.split('__', 2)[2][:-4], kTe=kT, F=F, E0=E0, nbands=calc.get_number_of_bands(),
                        top_band_occ=float(occ), seconds=round(time.time() - t0), pw=450, kdensity=3.5)
        write(f'done/{name}', out, format='extxyz')
        os.remove(path)
        flag = '  WARNING: increase bands' if occ > 1e-4 else ''
        print(f'{time.strftime("%H:%M")} {name}: F/atom {F/len(at):.4f} eV, top band occ {occ:.1e}, {time.time()-t0:.0f} s{flag}', flush=True)
    except Exception:
        traceback.print_exc(); os.rename(path, f'failed/{name}')
    finally:
        if big and os.path.exists(LOCK): os.remove(LOCK)
