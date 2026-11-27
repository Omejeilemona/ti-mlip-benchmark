"""0D TTM-MD of a free Ti nanoparticle (uniform Te). Pots: eam, eamte (EAM + DFT Te correction). Run from nanoparticle/."""
import argparse, os, json, time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from lammps import lammps

KB = 8.617333262e-5           # eV/K
MASS = 47.867                 # g/mol
MVV2E = 1.0364269e-4          # (g/mol)(A/ps)^2 -> eV  (LAMMPS metal units)
POTS = {'eam': ['pair_style eam/fs', 'pair_coeff * * {dir}/Ti1_Mendelev2016.eam.fs Ti'],
        'eamte': ['pair_style hybrid/overlay eam/fs snap', 'pair_coeff * * eam/fs {dir}/Ti1_Mendelev2016.eam.fs Ti Ti',
                  'pair_coeff * * snap {dir}/te_snap/Ti_te.snapcoeff {dir}/te_snap/Ti_te.snapparam Te00 Te00']}
NTYPES = {'eamte': 2}
MINFRAG = 13                  # atoms in fragments smaller than this get no electron-temperature correction
RFLOW, RCLUS = 6.0, 3.6       # A: local-flow averaging radius; bond cutoff for fragment analysis

ap = argparse.ArgumentParser()
ap.add_argument('--pot', default='eam', choices=['eam', 'eamte'])
ap.add_argument('--potdir', default='../potentials')
ap.add_argument('--d', type=float, default=2.0, help='particle diameter (nm)')
ap.add_argument('--eabs', type=float, default=2.0, help='absorbed energy (eV/atom)')
ap.add_argument('--fwhm', type=float, default=0.5, help='pulse FWHM (ps)')
ap.add_argument('--tpeak', type=float, default=1.0, help='pulse maximum (ps)')
ap.add_argument('--tend', type=float, default=30.0, help='simulated time after equilibration (ps)')
ap.add_argument('--T0', type=float, default=300.0)
ap.add_argument('--G', type=float, default=None, help='constant electron-phonon coupling (W m^-3 K^-1)')
ap.add_argument('--gtable', default='G_Ti_LinZhigilei.dat', help='count line, then Te(K) G(eV/ps/A^3/K)')
ap.add_argument('--ce', default='Ce_Ti_LinZhigilei.dat', help='Lin-Zhigilei Ce table: count line, then Te(K) Ce(eV/K/A^3)')
ap.add_argument('--dt', type=float, default=0.001, help='MD time step (ps)')
ap.add_argument('--ncouple', type=int, default=10, help='MD steps between TTM updates')
ap.add_argument('--seed', type=int, default=12345)
ap.add_argument('--out', default='run')
ap.add_argument('--dumpat', default='', help='comma-separated times (ps) for snapshots, e.g. 0,2,5,10,20,30')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
json.dump(vars(a), open(f'{a.out}/settings.json', 'w'), indent=1)

n0 = 4.506e3 / (MASS * 1.66054e-27)            # atoms per m^3 (rho = 4.506 g/cm^3)
ce_T, ce_v = np.loadtxt(a.ce, skiprows=1, unpack=True); ce_v = ce_v / (n0 * 1e-30)   # eV/K/atom
def c_e(T): return np.interp(T, ce_T, ce_v)
trapz = getattr(np, 'trapezoid', None) or np.trapz
def E_e(T):                                   # electronic internal energy per atom relative to T0
    Ts = np.linspace(a.T0, T, 200) if T != a.T0 else np.array([a.T0, a.T0])
    return trapz(c_e(Ts), Ts)
if a.G is not None:
    gT, gv = np.array([1.0, 1e6]), np.full(2, a.G / n0 / 1.602176634e-19 * 1e-12)   # eV/(ps K atom)
else:
    gT, gv = np.loadtxt(a.gtable, skiprows=1, unpack=True); gv = gv / (n0 * 1e-30)        # per atom
def g(T): return np.interp(T, gT, gv)
sig = a.fwhm / 2.3548
def s(t): return a.eabs * np.exp(-0.5 * ((t - a.tpeak) / sig)**2) / (sig * np.sqrt(2*np.pi))   # eV/ps/atom

# ---------------- build and equilibrate the particle ----------------
L = lammps(cmdargs=['-log', f'{a.out}/log.lammps', '-screen', 'none'])
R = 10 * a.d / 2
L.commands_list(['units metal', 'atom_style atomic', 'boundary s s s', 'atom_modify map array sort 0 0',
    f'lattice custom 1.0 a1 2.95 0 0 a2 -1.475 2.5548 0 a3 0 0 4.68 '
    'basis 0.3333333 0.6666667 0.25 basis 0.6666667 0.3333333 0.75',
    f'region box block -{R+3} {R+3} -{R+3} {R+3} -{R+3} {R+3} units box',
    f'create_box {NTYPES.get(a.pot, 1)} box', f'region ball sphere 0 0 0 {R} units box',
    'create_atoms 1 region ball'] + [f'mass {i} {MASS}' for i in range(1, NTYPES.get(a.pot, 1) + 1)] +
    [c.format(dir=a.potdir) for c in POTS[a.pot]] +
    ['neighbor 1.0 nsq', 'neigh_modify every 1 delay 0 check yes',
     'compute pe all pe/atom', 'compute ack all ackland/atom',
     f'velocity all create {2*a.T0} {a.seed} mom yes rot yes dist gaussian',
     f'timestep {a.dt}', 'thermo 1000',
     'minimize 1e-8 1e-10 1000 10000',
     f'fix eq all langevin {a.T0} {a.T0} 0.2 {a.seed}', 'fix nve all nve', 'run 10000',
     'unfix eq', 'velocity all zero linear', 'velocity all zero angular', 'run 5000'])
N = L.get_natoms()

def state():
    x = L.numpy.extract_atom('x')[:N]; v = L.numpy.extract_atom('v')[:N]
    tree = cKDTree(x); nb = tree.query_ball_point(x, RFLOW)
    vf = np.array([v[j].mean(0) for j in nb]); nn = np.array([len(j) for j in nb])
    vT = v - vf
    corr = np.where(nn > 1, nn / np.maximum(nn - 1, 1), 0.0)       # neighbour-average dof correction
    KT = 0.5 * MASS * MVV2E * (vT**2).sum(1)
    coupled = nn > 1
    Tl = 2 * (KT * corr)[coupled].sum() / (3 * KB * max(coupled.sum(), 1))
    return x, v, vf, vT, KT, coupled, Tl

def fragments(x, labels=False):
    pairs = cKDTree(x).query_pairs(RCLUS, output_type='ndarray')
    m = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(N, N))
    nc, lab = connected_components(m, directed=False); sizes = np.bincount(lab)
    if labels: return sizes[lab]
    return nc, sizes.max(), (sizes == 1).sum()

cur_bin, Esw = 0, 0.0
ESP = np.zeros(3)                                   # E_sw parts: monomers, clusters 2-12, Te bin (eV)
def update_te_potential(x, Te):
    """eamte: set the Te bin and per-atom types in three steps, recording the energy change of each."""
    global cur_bin
    b = int(min(max(round(KB * Te / 0.1), 0), 20))
    size = fragments(x, labels=True)
    newtype = np.where(size >= MINFRAG, 1, 2).astype(np.int32)
    typ = L.numpy.extract_atom('type')
    if b == cur_bin and np.array_equal(typ[:N], newtype): return 0.0
    def pe():
        L.command('run 0 pre yes post no'); return L.get_thermo('pe')
    e0 = pe(); de = np.zeros(3)
    for k, sel in enumerate((size == 1, size > 1)):    # 0: monomers, 1: everything else (clusters, rejoining)
        ch = sel & (typ[:N] != newtype)
        if ch.any():
            typ = L.numpy.extract_atom('type'); typ[np.where(ch)[0]] = newtype[ch]; e1 = pe(); de[k] = e1 - e0; e0 = e1
    if b != cur_bin:
        L.command(f'pair_coeff * * snap {a.potdir}/te_snap/Ti_te.snapcoeff {a.potdir}/te_snap/Ti_te.snapparam Te{b:02d} Te00')
        cur_bin = b; e1 = pe(); de[2] = e1 - e0
    ESP[:] += de
    return de.sum()

# ---------------- TTM-MD ----------------
Te, t, Eabs, Elat = a.T0, 0.0, 0.0, 0.0
dumps = sorted(float(x) for x in a.dumpat.split(',') if x.strip())
nsub = 50
L.command('run 0'); E0_atoms = L.get_thermo('etotal')
fs = open(f'{a.out}/esw_split.csv', 'w'); fs.write('t_ps,Esw_monomers,Esw_clusters,Esw_Tebin,Esw_total\n')
f = open(f'{a.out}/ttm.csv', 'w')
f.write('t_ps,Te_K,Tl_K,Eabs_eV_per_atom,Ee_eV_per_atom,Eatoms_eV_per_atom,Econs_eV_per_atom,'
        'n_fragments,largest,monomers,frac_hcp,frac_bcc,frac_fcc,frac_ico,frac_other,Rg_A,Te_bin,Eswitch_eV_per_atom\n')
nchunk = int(round(a.tend / (a.dt * a.ncouple))); t0 = time.time()
for k in range(nchunk + 1):
    x, v, vf, vT, KT, coupled, Tl = state()
    if dumps and abs(t - dumps[0]) < 0.5 * a.dt * a.ncouple:
        L.command('run 0 pre no post no')
        L.command(f'write_dump all custom {a.out}/snap_{dumps.pop(0):06.2f}ps.dump id x y z vx vy vz c_ack c_pe modify sort id')
    if k % 10 == 0:                                   # output every 10 chunks (0.1 ps)
        L.command('run 0 pre no post no')
        Eat = (L.get_thermo('etotal') - E0_atoms) / N
        ack = L.numpy.extract_compute('ack', 1, 1)[:N]
        fr = [np.mean(ack == i) for i in (3, 2, 1, 4, 0)]   # ackland: 0 unk, 1 bcc, 2 fcc, 3 hcp, 4 ico
        nc, big, mono = fragments(x)
        rg = np.sqrt(((x - x.mean(0))**2).sum(1).mean())
        Ee = E_e(Te)
        f.write(f'{t:.3f},{Te:.1f},{Tl:.1f},{Eabs:.4f},{Ee:.4f},{Eat:.4f},{Eat + Ee - Eabs - Esw / N:.4f},'
                f'{nc},{big},{mono},' + ','.join(f'{q:.3f}' for q in fr) + f',{rg:.2f},{cur_bin},{Esw / N:.4f}\n'); f.flush()
        fs.write(f'{t:.3f},' + ','.join(f'{q / N:.4f}' for q in ESP) + f',{Esw / N:.4f}\n'); fs.flush()
    if k == nchunk: break
    dtc = a.dt * a.ncouple; h = dtc / nsub; dE = 0.0     # electron ODE with sub-steps, Tl frozen in chunk
    fc = coupled.sum() / N
    for i in range(nsub):
        ts = t + (i + 0.5) * h
        q = g(Te) * (Te - Tl) * h; dE += q
        Te += (s(ts) * h - q * fc) / c_e(Te); Eabs += s(ts) * h      # electrons lose only what coupled atoms receive
    # give dE per coupled atom to the lattice by rescaling thermal velocities (exact energy)
    dEtot = dE * coupled.sum(); A = KT[coupled].sum()
    B = MASS * MVV2E * (vf[coupled] * vT[coupled]).sum()
    lam = (-B + np.sqrt(max(B*B + 4*A*(A + B + dEtot), 0.0))) / (2*A) if A > 0 else 1.0
    vnew = v.copy(); vnew[coupled] = vf[coupled] + lam * vT[coupled]
    L.numpy.extract_atom('v')[:N] = vnew
    if a.pot == 'eamte': Esw += update_te_potential(x, Te)
    L.command(f'run {a.ncouple} pre no post no')
    t += dtc
f.close()
L.command(f'write_dump all custom {a.out}/final.dump id x y z vx vy vz c_ack modify sort id')
print(f'{a.out}: {N} atoms, {a.tend} ps in {time.time()-t0:.0f} s')
