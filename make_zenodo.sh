#!/bin/bash
# Build the Zenodo release folder and zip from the ti-mlip-benchmark repository.
# Run from the repository root AFTER the v1.0 commit and tag:   bash make_zenodo.sh
# Output: zenodo_release/ and ti_mlip_zenodo_v1.0.zip (upload the zip contents, or the zip itself, to Zenodo).
set -u
R=zenodo_release
rm -rf "$R" ti_mlip_zenodo_v1.0.zip
missing=0
copy() {  # copy <source...> <destination folder>; warns about missing sources
    local dst="${@: -1}"; mkdir -p "$dst"
    for src in "${@:1:$#-1}"; do
        if compgen -G "$src" > /dev/null; then cp -r $src "$dst"/; else echo "MISSING: $src"; missing=$((missing+1)); fi
    done
}

# Description and license
mkdir -p "$R"; cp zenodo_README.md "$R/README.md"
cat > "$R/LICENSE.txt" <<'EOF'
Data and potentials in this record: Creative Commons Attribution 4.0 International (CC BY 4.0),
https://creativecommons.org/licenses/by/4.0/
Exceptions: mace_training_sets/ (Apache-2.0, derived from Roecken & Zavadlav 2024 / Wen et al. 2021)
and code_v1.0.zip (MIT, see the LICENSE file inside).
EOF

# Potentials (the third-party EAM file is not included)
for p in snap pod_v1 pod_v2_omega pod_v3 pod_v4 mace mace_v3 te_snap te_snap_res; do
    copy "potentials/$p" "$R/potentials"
done
rm -f "$R"/potentials/*/Ti_training_analysis.pod "$R"/potentials/*/Ti_test_analysis.pod
compgen -G "$R/potentials/mace_v3/*.model" > /dev/null || { echo "MISSING: potentials/mace_v3/*.model (download from Google Drive)"; missing=$((missing+1)); }

# Ground-state DFT
copy "dft_v3/done/*.xyz" "$R/dft_groundstate_v3"
copy "dft_v3/logs/*.txt" "$R/dft_groundstate_v3/logs"
copy "data_omega/*.xyz" "$R/dft_omega"

# Finite-Te DFT
copy "dft_te/done/*.xyz" "$R/dft_finite_te/done"
copy "data_te/*" "$R/dft_finite_te/by_temperature"
copy "dft_te/logs/Te*.txt" "$R/dft_finite_te/logs"
n=$(ls "$R"/dft_finite_te/done/*.xyz 2>/dev/null | wc -l); echo "finite-Te calculations: $n (expected 414)"

# MACE training sets
copy "data_mace/*.xyz" "$R/mace_training_sets/compact"
copy "mace_v3/*.xyz" "$R/mace_training_sets/v3"

# Results
copy "md_te/*.csv" "$R/results/constant_te_md"
copy md_te_summary.txt md_te_res_summary.txt "$R/results/constant_te_md"
for d in nanoparticle/eam_d2_* nanoparticle/eamte_d2_* nanoparticle/split_e*; do
    mkdir -p "$R/results/$d"; cp "$d"/*.csv "$d"/settings.json "$R/results/$d"/ 2>/dev/null
done
copy "melt_v3/*_T.dat" "$R/results/melting_pod_v3"
copy "results/mace_v3/*" "$R/results/mace_v3"
copy "results/*.json" "results/*.txt" "results/*.dat" "results/*.pod" ce_hcp_F_vs_Te.txt "$R/results/properties"

# TTM input tables
copy nanoparticle/Ce_Ti_LinZhigilei.dat nanoparticle/G_Ti_LinZhigilei.dat "$R/ttm_inputs"
[ -f nanoparticle/ce_LinZhigilei_per_atom.txt ] && cp nanoparticle/ce_LinZhigilei_per_atom.txt "$R/ttm_inputs/"

# Code snapshot of the tagged release
if git rev-parse -q --verify v1.0 > /dev/null; then
    git archive --format=zip --prefix=ti-mlip-benchmark-v1.0/ -o "$R/code_v1.0.zip" v1.0
else
    echo "MISSING: git tag v1.0 (create it first: git tag -a v1.0 -m 'Version submitted to Phys. Rev. Materials')"; missing=$((missing+1))
fi

# Checksums and zip
( cd "$R" && find . -type f ! -name MD5SUMS.txt -print0 | sort -z | xargs -0 md5sum > MD5SUMS.txt )
zip -qr ti_mlip_zenodo_v1.0.zip "$R"
echo "files: $(find "$R" -type f | wc -l), size: $(du -sh "$R" | cut -f1), zip: $(du -sh ti_mlip_zenodo_v1.0.zip | cut -f1)"
[ $missing -eq 0 ] && echo "OK: nothing missing" || echo "WARNING: $missing item(s) missing, see the MISSING lines above"
