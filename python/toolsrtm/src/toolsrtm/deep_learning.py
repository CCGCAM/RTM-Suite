"""Deep-learning trait inversion: dense ("Hidden-layers") and 1D-CNN Keras
models with a configurable optimizer, matching R's ``getMLmodel`` /
``getMLmodel.withRetrain``.

Needs the optional ``dl`` extra (``pip install toolsrtm[dl]``: tensorflow).
Like :mod:`toolsrtm.inversion`, nothing here is imported by
``toolsrtm/__init__.py`` and TensorFlow is imported lazily inside
:func:`get_ml_model`, so a plain ``import toolsrtm`` never requires it.

Unlike R's own (non-reproducible, GPU/BLAS-order-dependent) Keras training,
this is not verified to floating-point precision against R -- what's
verified is that both architectures train to convergence and produce sane
held-out R^2/RMSE on synthetic data (see ``tests/test_deep_learning.py``),
the same standard already used for
``Scripts/Python/*/3_inversion_dl.py``/``4_inversion_dl.py``, which this
module formalizes into an installable, tested package function.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np


_OPTIMIZERS = {
    "adam": lambda tf: tf.keras.optimizers.Adam(learning_rate=1e-4, beta_1=0.9, beta_2=0.999),
    "adadelta": lambda tf: tf.keras.optimizers.Adadelta(learning_rate=1.0, rho=0.95),
    "adagrad": lambda tf: tf.keras.optimizers.Adagrad(learning_rate=0.01),
    "adamax": lambda tf: tf.keras.optimizers.Adamax(learning_rate=0.002, beta_1=0.9, beta_2=0.999),
    "nadam": lambda tf: tf.keras.optimizers.Nadam(learning_rate=0.002, beta_1=0.9, beta_2=0.999),
    "rmsprop": lambda tf: tf.keras.optimizers.RMSprop(learning_rate=0.001, rho=0.9),
    "sgd": lambda tf: tf.keras.optimizers.SGD(learning_rate=0.01, momentum=0.0, nesterov=False),
}


@dataclass
class MLModelResult:
    model: object  #: the fitted `keras.Model`
    history: dict  #: per-epoch training history (`keras.callbacks.History.history`)
    stats: dict  #: {"r2":.., "rmse":..} on the held-out validation split
    predictions: dict  #: {"y_true": np.ndarray, "y_pred": np.ndarray} on the held-out validation split
    x_scaler: object  #: fitted `sklearn.preprocessing.StandardScaler` for the predictors


def get_ml_model(
    dataset,
    dep_var: str,
    model: Literal["Hidden-layers", "CNN"] = "Hidden-layers",
    optimizer: str = "adam",
    batch_size: int = 125,
    n_epochs: int = 100,
    prop_split: tuple[float, float] = (0.8, 0.2),
    n_layers: int = 3,
    n_neurons: int = 64,
    n_times: int = 1,
    seed: int = 123,
    verbose: int = 0,
) -> MLModelResult:
    """Train a dense or 1D-CNN Keras regression model to predict ``dep_var``
    from every other column of ``dataset``.

    Python port of ``getMLmodel``/``getMLmodel.withRetrain`` (R). Predictors
    are standardized (`sklearn.preprocessing.StandardScaler`) before
    training, matching R's own ``data.trans='preProcess'`` default; the
    response is left on its original scale (matching R's own
    ``depVar.trans=FALSE`` default).

    :param dataset: `pandas.DataFrame` containing ``dep_var`` and predictor columns.
    :param dep_var: name of the column to predict.
    :param model: ``"Hidden-layers"`` (dense MLP: ``n_layers`` hidden layers
        of ``n_neurons`` units, ReLU, dropout 0.1 after the first hidden
        layer, matching R's 3-layer 64/32(dropout)/16 default when
        ``n_layers=3, n_neurons=64``) or ``"CNN"`` (1D convolution over the
        predictor vector: conv(64,k=4) -> pool -> conv(32,k=2) -> pool ->
        dense(16) -> dropout(0.1) -> output).
    :param optimizer: one of ``"adam"``, ``"adadelta"``, ``"adagrad"``,
        ``"adamax"``, ``"nadam"``, ``"rmsprop"``, ``"sgd"`` (same learning
        rates/momenta as the R defaults for each).
    :param batch_size: training batch size.
    :param n_epochs: maximum training epochs (early stopping on
        ``val_loss``, patience 5, restores best weights -- matches R).
    :param prop_split: ``(train_fraction, val_fraction)``.
    :param n_layers: number of hidden layers for ``"Hidden-layers"`` (ignored for ``"CNN"``).
    :param n_neurons: units in the first hidden layer for ``"Hidden-layers"``
        (subsequent layers halve down to a floor of 8; ignored for ``"CNN"``).
    :param n_times: fit this many times with different random initializations
        and keep the run with the lowest validation loss (matches
        ``getMLmodel.withRetrain``'s ``n.times``).
    :param seed: random seed for the train/val split and Keras initialization.
    :param verbose: Keras ``fit()`` verbosity (0, 1, or 2).
    :return: :class:`MLModelResult`.
    """
    import tensorflow as tf
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler

    if optimizer not in _OPTIMIZERS:
        raise ValueError(f"Unknown optimizer {optimizer!r}. Choose from {list(_OPTIMIZERS)}.")
    if model not in ("Hidden-layers", "CNN"):
        raise ValueError("model must be 'Hidden-layers' or 'CNN'.")

    inputs = [c for c in dataset.columns if c != dep_var]
    X = dataset[inputs].to_numpy(dtype=float)
    y = dataset[dep_var].to_numpy(dtype=float)

    train_frac, val_frac = prop_split
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=val_frac / (train_frac + val_frac), random_state=seed
    )

    x_scaler = StandardScaler().fit(X_train)
    X_train_s = x_scaler.transform(X_train)
    X_val_s = x_scaler.transform(X_val)

    def _build_model(attempt_seed):
        tf.random.set_seed(attempt_seed)
        if model == "Hidden-layers":
            layers = [tf.keras.layers.Input(shape=(X_train_s.shape[1],))]
            units = n_neurons
            for i in range(n_layers):
                layers.append(tf.keras.layers.Dense(max(units, 8), activation="relu"))
                if i == 1:  # dropout after the 2nd hidden layer, matching R's Hidden-layers architecture
                    layers.append(tf.keras.layers.Dropout(0.1))
                units = max(units // 2, 8)
            # linear, not relu: a relu output can get permanently stuck at 0
            # (zero gradient too) if its pre-activation goes negative early
            # in training -- fixed on the R side in the same way, see
            # getMLmodel.R's comment on this exact issue.
            layers.append(tf.keras.layers.Dense(1, activation="linear"))
            net = tf.keras.Sequential(layers)
        else:  # CNN
            net = tf.keras.Sequential([
                tf.keras.layers.Input(shape=(X_train_s.shape[1], 1)),
                tf.keras.layers.Conv1D(filters=64, kernel_size=4, activation="relu"),
                tf.keras.layers.MaxPooling1D(pool_size=2),
                tf.keras.layers.Conv1D(filters=32, kernel_size=2, activation="relu"),
                tf.keras.layers.MaxPooling1D(pool_size=2),
                tf.keras.layers.Flatten(),
                tf.keras.layers.Dense(16, activation="relu"),
                tf.keras.layers.Dropout(0.1),
                tf.keras.layers.Dense(1, activation="linear"),  # see the Hidden-layers branch above
            ])
        net.compile(loss="mse", optimizer=_OPTIMIZERS[optimizer](tf), metrics=["mean_absolute_error"])
        return net

    if model == "CNN":
        X_train_fit = X_train_s.reshape(*X_train_s.shape, 1)
        X_val_fit = X_val_s.reshape(*X_val_s.shape, 1)
    else:
        X_train_fit, X_val_fit = X_train_s, X_val_s

    best_val_loss = np.inf
    best_net, best_history = None, None
    for attempt in range(max(1, n_times)):
        net = _build_model(seed + attempt)
        callback = tf.keras.callbacks.EarlyStopping(monitor="val_loss", mode="min", patience=5,
                                                      restore_best_weights=True)
        history = net.fit(X_train_fit, y_train, epochs=n_epochs, batch_size=batch_size, verbose=verbose,
                           shuffle=False, callbacks=[callback], validation_split=0.2)
        val_loss = min(history.history["val_loss"])
        if val_loss < best_val_loss:
            best_val_loss, best_net, best_history = val_loss, net, history

    y_pred = np.ravel(best_net.predict(X_val_fit, verbose=0))
    ss_res = np.sum((y_val - y_pred) ** 2)
    ss_tot = np.sum((y_val - y_val.mean()) ** 2)
    stats = {"r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
              "rmse": float(np.sqrt(np.mean((y_val - y_pred) ** 2)))}

    return MLModelResult(model=best_net, history=best_history.history, stats=stats,
                          predictions={"y_true": y_val, "y_pred": y_pred}, x_scaler=x_scaler)
