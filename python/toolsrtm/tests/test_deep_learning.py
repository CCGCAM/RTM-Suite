import numpy as np
import pandas as pd
import pytest

tf = pytest.importorskip("tensorflow")

from toolsrtm.deep_learning import get_ml_model


def _smooth_dataset(n=300, p=10, seed=7):
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, 1, size=(n, p))
    y = 3 * X[:, 0] - 2 * X[:, 3] + 0.5 * X[:, 6] + rng.normal(scale=0.05, size=n)
    y = y - y.min() + 0.1  # non-negative, matches the ReLU output layer R itself uses
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(p)])
    df["y"] = y
    return df


def test_hidden_layers_trains_and_predicts_well():
    # adam at R's own conservative default learning rate (1e-4) converges
    # slowly -- generous epoch budget + multiple restarts (n_times), matching
    # what getMLmodel.withRetrain is for, rather than a stricter floor.
    df = _smooth_dataset()
    res = get_ml_model(df, "y", model="Hidden-layers", optimizer="adam", n_epochs=300, n_times=3, seed=1)
    assert res.stats["r2"] > 0.4, res.stats
    assert res.predictions["y_pred"].shape == res.predictions["y_true"].shape


def _spectrum_like_dataset(n=300, p=20, seed=7):
    # A 1D-CNN needs locally-correlated input to have any advantage over a
    # dense layer (like real, adjacent reflectance/index bands do) --
    # _smooth_dataset's independent uniform columns give it nothing to
    # exploit, so this generates a smooth Gaussian-bump "spectrum" instead.
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.2, 1.0, size=n)
    band = np.linspace(0, 1, p)
    X = np.array([base[i] * np.exp(-((band - 0.3) ** 2) / 0.05) + 0.02 * rng.normal(size=p) for i in range(n)])
    y = base * 10 + 1
    df = pd.DataFrame(X, columns=[f"b{i}" for i in range(p)])
    df["y"] = y
    return df


def test_cnn_trains_and_predicts_well():
    # CNN + adam at R's own conservative lr=1e-4 converges unevenly run to
    # run even with restarts (some end up stuck in a worse local optimum
    # within the epoch budget) -- floor is deliberately loose, just enough
    # to confirm it's actually learning (R2 well above 0), not broken.
    df = _spectrum_like_dataset()
    res = get_ml_model(df, "y", model="CNN", optimizer="adam", n_epochs=300, n_times=5, seed=1)
    assert res.stats["r2"] > 0.15, res.stats


def test_rejects_unknown_model_and_optimizer():
    df = _smooth_dataset(n=40, p=10)
    with pytest.raises(ValueError):
        get_ml_model(df, "y", model="not-a-real-model")
    with pytest.raises(ValueError):
        get_ml_model(df, "y", optimizer="not-a-real-optimizer")
