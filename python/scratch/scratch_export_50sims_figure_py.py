"""One-off: generate the REAL Python output figure for the "50 simulations"
code example in docs/index.html. Run from the repo root:
    python python/scratch/scratch_export_50sims_figure_py.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from toolsrtm import inform

root = Path(__file__).resolve().parents[2]
outdir = root / "assets" / "examples"
outdir.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(42)
n = 50
wl = np.arange(400, 2501)
rsoil = np.full(2101, 0.15)

rows = []
for _ in range(n):
    rows.append(dict(
        N=rng.uniform(1.2, 2.5), Cab=rng.uniform(10, 70), Car=rng.uniform(5, 15),
        Anth=rng.uniform(0, 2), Cbrown=rng.uniform(0, 0.3), EWT=rng.uniform(0.005, 0.03),
        LMA=rng.uniform(0.003, 0.015), alpha=40, Prot=0, CBC=0,
        LAI=rng.uniform(0.5, 7.0), hspot=0.01, LIDFa=rng.uniform(-0.5, 0.5),
        LIDFb=rng.uniform(-0.3, 0.3), TypeLidf=1, tts=30, tto=0, psi=0,
        LAIu=rng.uniform(0.5, 3.0), sd=rng.uniform(200, 800), cd=rng.uniform(3, 10),
        h=rng.uniform(8, 20), skyl=0.1,
    ))

lai_values = [r["LAI"] for r in rows]
order = np.argsort(lai_values)
colors = plt.cm.viridis(np.linspace(0.05, 0.95, n))

fig, ax = plt.subplots(figsize=(9, 5.6), dpi=130)
for rank, idx in enumerate(order):
    refl = inform(rows[idx], rsoil, leaf_model="PROSPECT-PRO")
    ax.plot(wl, refl, color=colors[rank], lw=1)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Reflectance")
ax.set_ylim(0, 0.5)
ax.grid(color="#dddddd")
fig.tight_layout()
fig.savefig(outdir / "py_50sims.png")
plt.close(fig)
print("Saved py_50sims.png")
