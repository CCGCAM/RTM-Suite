toolsrtm.inversion
====================

Trait-inversion tools: CARS-PLS and VIF predictor selection, LUT nearest-
neighbour ("merit function") matching, and a 12-algorithm ML dispatcher
(PLSR/SVM/RF/GB/NN/Bayesian/AdaBag/BRNN/xGB/RVM/qLASSO/Ensemble) built on
scikit-learn/xgboost. Direct port of ``ToolsRTM::carspls``/``get.cars.pls``,
``getVIF``, ``get.inversionOpt``, ``get.inversion``, ``hybrid_inversion``/
``hybrid_inversionE``.

.. note::
   Needs the optional ``ml`` extra: ``pip install "toolsrtm[ml]"``. Nothing
   in this module is imported by ``toolsrtm/__init__.py``'s own import chain
   -- a plain ``import toolsrtm`` never requires scikit-learn/xgboost.

.. note::
   R's ``get.inversion``/``hybrid_inversion`` dispatch to specific ``caret``
   methods (``bartMachine``, ``rqlasso``, ``rvmLinear``, ``AdaBag``, ``brnn``,
   ...). See :data:`~toolsrtm.inversion.ALGORITHMS` for exactly which
   scikit-learn/xgboost estimator each algorithm name maps to, and, where
   there's no direct equivalent, what was substituted and why.

.. automodule:: toolsrtm.inversion
   :members:
   :undoc-members:
   :show-inheritance:
