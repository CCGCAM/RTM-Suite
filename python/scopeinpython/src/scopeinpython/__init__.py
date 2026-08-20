"""scopeinpython (SCOPEinPython): Python port of the R SCOPEinR package.

Ported: the BSM soil reflectance model, the RTMo optical (top-of-canopy
BRDF) radiative transfer pipeline (leaf+soil+geometry ->
rdd/rsd/rdo/rso/refl), driven by the ported ``toolsrtm`` leaf/canopy
models, and leaf-level biochemistry (Farquhar/Collatz photosynthesis +
van der Tol et al. 2014 fluorescence yield, ``get_biochemical``) given an
assumed leaf micro-environment. The full thermal energy-balance (ebal
iteration itself) and the fluorescence (RTMf)/Zeaxanthin (RTMz) canopy
radiative-transfer modules are NOT ported yet. See ``python/README.md`` at
the repo root for exact scope and numerical verification against the
original R package.
"""
from .biochemical import (
    BiochemResult,
    LeafBio,
    MeteoLeaf,
    get_ball_berry,
    get_biochemical,
    get_ci_next,
    get_compute_a,
    get_fluorescence_model,
    get_gs_fun,
    get_high_temp_inhibtion_c3,
    get_temperature_function_c3,
    sel_root,
)
from .fluspect import FluspectCxScopeResult, get_fluspect_cx_scope
from .fluspect_mscope import FluspectMScopeResult, MultiLayerLeafBio, fluspect_mscope
from .rtmf import RTMfResult, rtmf
from .rtmz import RTMzResult, rtmz
from .rtmt_sb import RTMtSbResult, rtmt_sb
from .ebal import EbalCanopyParams, EbalMeteo, EbalResult, EbalSoilParams, aggregator_ebal, ebal
from .soil import SoilParams, WettingParams, get_bsm, soilwat
from .thermal import (
    HeatFluxes,
    ResistanceParams,
    Resistances,
    get_heatfluxes,
    get_resistances,
    monin_obukhov,
    phstar,
    psih,
    psim,
    stefan_boltzmann,
)
from .spectral import SpectralConfig, get_spectra_scope
from .rtmo import CanopyStructure, NetRadiationLite, RTMoResult, net_radiation_lite, run_rtmo
from .scope import ScopeOptions, ScopeResult, TDP_DEFAULT, get_scope, get_zo_and_d

__version__ = "0.1.0"

__all__ = [
    "SoilParams",
    "WettingParams",
    "get_bsm",
    "soilwat",
    "SpectralConfig",
    "get_spectra_scope",
    "CanopyStructure",
    "RTMoResult",
    "run_rtmo",
    "NetRadiationLite",
    "net_radiation_lite",
    "BiochemResult",
    "LeafBio",
    "MeteoLeaf",
    "get_ball_berry",
    "get_biochemical",
    "get_ci_next",
    "get_compute_a",
    "get_fluorescence_model",
    "get_gs_fun",
    "get_high_temp_inhibtion_c3",
    "get_temperature_function_c3",
    "sel_root",
    "FluspectCxScopeResult",
    "get_fluspect_cx_scope",
    "FluspectMScopeResult",
    "MultiLayerLeafBio",
    "fluspect_mscope",
    "RTMfResult",
    "rtmf",
    "RTMzResult",
    "rtmz",
    "monin_obukhov",
    "psim",
    "psih",
    "phstar",
    "ResistanceParams",
    "Resistances",
    "get_resistances",
    "HeatFluxes",
    "get_heatfluxes",
    "stefan_boltzmann",
    "RTMtSbResult",
    "rtmt_sb",
    "EbalCanopyParams",
    "EbalMeteo",
    "EbalResult",
    "EbalSoilParams",
    "aggregator_ebal",
    "ebal",
    "ScopeOptions",
    "ScopeResult",
    "TDP_DEFAULT",
    "get_scope",
    "get_zo_and_d",
]
