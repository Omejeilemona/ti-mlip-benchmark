"""Create ttm0d_split.py: same as ttm0d.py, but the switching energy E_sw is split into
(1) retyping of monomers (force-free offset), (2) retyping of clusters of 2-12 atoms (and atoms rejoining),
(3) change of the Te bin for the atoms in the particle. Written to {out}/esw_split.csv every 0.1 ps.
Run from nanoparticle/:  python3 make_split.py"""
s = open('ttm0d.py').read()
a = s.index('def update_te_potential'); b = s.index('# ---------------- TTM-MD')
new = '''ESP = np.zeros(3)                                   # E_sw parts: monomers, clusters 2-12, Te bin (eV)
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

'''
s = s[:a] + new + s[b:]
s = s.replace("f = open(f'{a.out}/ttm.csv', 'w')",
              "fs = open(f'{a.out}/esw_split.csv', 'w'); fs.write('t_ps,Esw_monomers,Esw_clusters,Esw_Tebin,Esw_total\\n')\n"
              "f = open(f'{a.out}/ttm.csv', 'w')")
s = s.replace("f',{rg:.2f},{cur_bin},{Esw / N:.4f}\\n'); f.flush()",
              "f',{rg:.2f},{cur_bin},{Esw / N:.4f}\\n'); f.flush()\n"
              "        fs.write(f'{t:.3f},' + ','.join(f'{q / N:.4f}' for q in ESP) + f',{Esw / N:.4f}\\n'); fs.flush()")
assert s.count('esw_split.csv') == 1 and s.count('fs.write(f') == 1
open('ttm0d_split.py', 'w').write(s); print('wrote ttm0d_split.py')
