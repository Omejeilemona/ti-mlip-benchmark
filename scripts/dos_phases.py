"""Electronic density of states of hcp, bcc, fcc and omega Ti (same cells as the finite-Te EOS set, V0) and the
electronic specific-heat coefficient gamma = (pi^2/3) kB^2 g(EF). Sommerfeld estimate of the change of
F(phase) - F(hcp) with Te: dF ~ -(pi^2/6) (kB Te)^2 [g_phase(EF) - g_hcp(EF)].
Run from the repository root:  python3 scripts/dos_phases.py [k-point density, default 8]
Writes results/dos_{phase}.dat, results/dos_summary.txt; then scripts/fig_dos.py makes the figure."""
import sys, glob, os, time, numpy as np
from ase.io import read
from gpaw import GPAW, PW, FermiDirac
KD = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
SIG, E = 0.10, np.linspace(-8, 6, 1401)            # Gaussian broadening (eV), energy grid relative to EF
os.makedirs('results', exist_ok=True); summ = []
for ph in ('hcp', 'bcc', 'fcc', 'omega'):
    t0 = time.time(); a = read(glob.glob(f'dft_te/done/Te0.1__*__A_eos_{ph}_v1.00.xyz')[0]); a.calc = None
    a.calc = GPAW(mode=PW(450), xc='PBE', kpts={'density': KD, 'gamma': True}, occupations=FermiDirac(0.02),
                  nbands='150%', txt=f'results/dos_{ph}.gpaw.txt')
    a.get_potential_energy(); c = a.calc; ef = c.get_fermi_level(); w = c.get_k_point_weights()
    g = np.zeros_like(E)
    for k, wk in enumerate(w):
        eps = c.get_eigenvalues(kpt=k) - ef
        g += 2 * wk * np.exp(-0.5 * ((E[:, None] - eps[None, :]) / SIG) ** 2).sum(1) / (SIG * np.sqrt(2 * np.pi))
    g /= len(a)                                     # states / eV / atom, both spins
    np.savetxt(f'results/dos_{ph}.dat', np.c_[E, g], header='E-EF(eV) g(states/eV/atom)')
    gEF = np.interp(0, E, g); gam = np.pi ** 2 / 3 * (8.617333e-5) ** 2 * gEF * 96485.33 * 1e3
    summ.append((ph, gEF, gam)); print(f'{ph:5s}: g(EF) = {gEF:.3f} states/eV/atom, gamma = {gam:.2f} mJ/mol/K^2 '
                                       f'({len(w)} k-points, {time.time() - t0:.0f} s)', flush=True)
gh = summ[0][1]
lines = [f'{p:5s} g(EF) = {g:.3f} states/eV/atom   gamma = {gm:.2f} mJ/mol/K^2' for p, g, gm in summ]
lines.append('Sommerfeld estimate of the change of F(phase) - F(hcp) from kTe = 0 (meV/atom):')
for p, g, _ in summ[1:]:
    lines.append(f'  {p:5s} ' + '  '.join(f'{t:.1f} eV: {-1000 * np.pi ** 2 / 6 * t ** 2 * (g - gh):7.1f}' for t in (0.1, 0.2, 0.3)))
open('results/dos_summary.txt', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
