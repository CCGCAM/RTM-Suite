"""Regression evaluation metrics shared by the classical (3_inversion_ml.py)
and deep-learning (4_inversion_dl.py) inversion scripts.

RMSE/NRMSE/NMB/FGE match the cost functions used in the course material
this pipeline mirrors (RMSE = sqrt(mean((sim-obs)^2)); NRMSE = RMSE scaled
by the observed data range; NMB = normalized mean bias; FGE = fractional
gross error).
"""
from __future__ import annotations

import numpy as np


def rmse(obs: np.ndarray, sim: np.ndarray) -> float:
    obs, sim = np.asarray(obs, dtype=float), np.asarray(sim, dtype=float)
    return float(np.sqrt(np.mean((sim - obs) ** 2)))


def nrmse(obs: np.ndarray, sim: np.ndarray) -> float:
    """RMSE normalized by the observed value range."""
    obs = np.asarray(obs, dtype=float)
    obs_range = obs.max() - obs.min()
    return float(rmse(obs, sim) / obs_range) if obs_range > 0 else float("nan")


def nmb(obs: np.ndarray, sim: np.ndarray) -> float:
    """Normalized mean bias: (mean(sim) - mean(obs)) / mean(obs)."""
    obs, sim = np.asarray(obs, dtype=float), np.asarray(sim, dtype=float)
    return float((sim.mean() - obs.mean()) / obs.mean()) if obs.mean() != 0 else float("nan")


def fge(obs: np.ndarray, sim: np.ndarray) -> float:
    """Fractional gross error: mean(2*|sim-obs|/(|sim|+|obs|))."""
    obs, sim = np.asarray(obs, dtype=float), np.asarray(sim, dtype=float)
    denom = np.abs(sim) + np.abs(obs)
    valid = denom > 0
    return float(np.mean(2 * np.abs(sim[valid] - obs[valid]) / denom[valid]))


def r2(obs: np.ndarray, sim: np.ndarray) -> float:
    obs, sim = np.asarray(obs, dtype=float), np.asarray(sim, dtype=float)
    ss_res = np.sum((obs - sim) ** 2)
    ss_tot = np.sum((obs - obs.mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def all_metrics(obs: np.ndarray, sim: np.ndarray) -> dict:
    return {"R2": r2(obs, sim), "RMSE": rmse(obs, sim), "NRMSE": nrmse(obs, sim), "NMB": nmb(obs, sim), "FGE": fge(obs, sim)}
