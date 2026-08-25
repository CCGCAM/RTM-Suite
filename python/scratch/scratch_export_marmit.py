"""One-time export of MARMIT reference data (R side -> Python package data).

Reads ToolsRTM/inst/extdata/marmit/ directly (no R needed -- these are plain
CSV/tab-separated text files) and writes consolidated CSVs into
python/toolsrtm/src/toolsrtm/data/marmit/. Only bundles the driest spectrum
per soil ID (17 for Bablet_2016) since that's the only file
get.marmit.rsoil()/get_marmit_rsoil() ever actually reads -- MARMIT computes
wet reflectance FROM that one dry reference, it doesn't need the other ~90
wetter-condition spectra in the database. Mirrors the R side's own choice to
bundle only Bablet_2016 (of 8 available databases) to stay lightweight.

Run once from repo root: python python/scratch/scratch_export_marmit.py
"""
import csv
import os

REPO_ROOT = os.path.dirname(os.path.abspath(__file__)) + "/.."
R_MARMIT = os.path.join(REPO_ROOT, "ToolsRTM", "inst", "extdata", "marmit")
OUT_DIR = os.path.join(REPO_ROOT, "python", "toolsrtm", "src", "toolsrtm", "data", "marmit")
os.makedirs(OUT_DIR, exist_ok=True)


def read_tab_csv(path):
    with open(path, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = [row for row in reader if row]
    return header, rows


# --- 1. Water optics: n_segelstein.csv + alpha_buikouwie.csv -> water_optics.csv ---
_, n_rows = read_tab_csv(os.path.join(R_MARMIT, "parameters", "n_segelstein.csv"))
_, alpha_rows = read_tab_csv(os.path.join(R_MARMIT, "parameters", "alpha_buikouwie.csv"))
n_by_wl = {float(r[0]): float(r[1]) for r in n_rows}
alpha_by_wl = {float(r[0]): float(r[1]) for r in alpha_rows}
wls = sorted(set(n_by_wl) & set(alpha_by_wl))

with open(os.path.join(OUT_DIR, "water_optics.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Wvl", "n", "alpha"])
    for wl in wls:
        w.writerow([wl, n_by_wl[wl], alpha_by_wl[wl]])

# --- 2. Bablet_2016: driest spectrum per ID + K/a/psi metadata ---
db_dir = os.path.join(R_MARMIT, "databases", "Bablet_2016")
# The index file is comma-separated (unlike the tab-separated spectra files).
with open(os.path.join(db_dir, "Bablet_2016.csv"), "r", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    index_rows = list(reader)

by_id = {}
for row in index_rows:
    sid = int(row["ID"])
    smcg = float(row["SMCg"])
    if sid not in by_id or smcg < by_id[sid]["SMCg"]:
        by_id[sid] = {
            "SMCg": smcg, "Refl_file": row["Refl_file"], "Name": row["Name"],
            "K": row["K"], "a": row["a"], "psi": row["psi"],
        }

with open(os.path.join(OUT_DIR, "bablet_2016_index.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ID", "Name", "K", "a", "psi"])
    for sid in sorted(by_id):
        rec = by_id[sid]
        w.writerow([sid, rec["Name"], rec["K"], rec["a"], rec["psi"]])

with open(os.path.join(OUT_DIR, "bablet_2016_spectra.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["ID", "Wvl", "R"])
    for sid in sorted(by_id):
        spec_path = os.path.join(db_dir, "spectra", by_id[sid]["Refl_file"])
        _, spec_rows = read_tab_csv(spec_path)
        for r in spec_rows:
            w.writerow([sid, r[0], r[1]])

print(f"Wrote {len(wls)} water-optics rows, {len(by_id)} soil IDs to {OUT_DIR}")
