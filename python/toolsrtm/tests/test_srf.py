import numpy as np
import pytest

from toolsrtm.srf import (
    fwhm_prisma,
    spectral_convolution_srf,
    srf_prisma,
    srf_sentinel2a,
    srf_sentinel2b,
    enmap_characteristics,
    sensor_characteristics,
    spectral_convolution_gaussian,
)


def _toy_spectrum():
    wave = np.arange(400, 2501)
    values = 0.05 + 0.3 * wave / 2500
    return wave, values


def test_srf_prisma_matches_r_reference():
    # Cross-checked against ToolsRTM::get.spectral.convolution.srf() (R) run
    # with the same toy spectrum -- first 3 bands, band-center wavelength
    # (SRF-weighted mean) and convolved reflectance.
    wave, values = _toy_spectrum()
    prisma = srf_prisma()
    fwhm = fwhm_prisma()
    res = spectral_convolution_srf(wave, values, prisma, fwhm=fwhm)

    assert len(res.band_names) == 234
    assert not np.isnan(res.wl).any()
    assert not np.isnan(res.rfl).any()
    np.testing.assert_allclose(res.wl[:3], [410.53720265, 416.02688819, 423.78502123], atol=1e-4)
    np.testing.assert_allclose(res.rfl[:3], [0.09926446, 0.09992323, 0.10085420], atol=1e-6)
    # PRISMA's bundled FWHM must be used verbatim (not the coarser
    # SRF-derived estimate) when supplied.
    np.testing.assert_allclose(res.fwhm[:3], [11.352248, 10.377187, 9.750846], atol=1e-4)


def test_srf_sentinel2a_matches_r_reference():
    wave, values = _toy_spectrum()
    s2a = srf_sentinel2a()
    res = spectral_convolution_srf(wave, values, s2a)

    assert res.band_names == ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"]
    np.testing.assert_allclose(
        res.wl,
        [442.6950461, 492.71521346, 559.84905548, 664.62175292, 704.11493622, 740.49182088,
         782.75291733, 832.79041114, 864.71078924, 945.05447044, 1373.45554874, 1613.68050304, 2202.36780085],
        atol=1e-3,
    )
    np.testing.assert_allclose(res.rfl[0], 0.10312341, atol=1e-6)
    # No bundled FWHM table for Sentinel-2A -> empirical half-max-crossing
    # estimate, an integer number of nm.
    assert not np.isnan(res.fwhm).any()


def test_srf_sentinel2a_and_2b_are_real_distinct_satellites():
    # Sentinel-2A/B are twin satellites with genuinely different measured
    # SRFs, not two copies of the same curve -- confirmed against ToolsRTM
    # (max band-center difference ~17nm, in SWIR2/B12).
    wave, values = _toy_spectrum()
    res_a = spectral_convolution_srf(wave, values, srf_sentinel2a())
    res_b = spectral_convolution_srf(wave, values, srf_sentinel2b())

    assert res_a.band_names == res_b.band_names
    diffs = np.abs(res_a.wl - res_b.wl)
    assert diffs.max() > 10  # B12 differs by ~17nm
    assert diffs.min() > 0  # no band is bit-identical between A and B


def test_spectral_convolution_gaussian_enmap_matches_r_reference():
    # Cross-checked against ToolsRTM::get.spectral.convolution.gaussian()
    # (R), same toy spectrum, first 6 EnMAP bands.
    wave, values = _toy_spectrum()
    res = spectral_convolution_gaussian(wave, values, sensor="EnMAP")

    assert len(res.wl) == 242
    assert not np.isnan(res.rfl).any()
    np.testing.assert_allclose(res.wl[:6], [423.03, 428.80, 434.29, 439.58, 444.72, 449.75], atol=1e-2)
    np.testing.assert_allclose(
        res.rfl[:6], [0.1007636, 0.1014560, 0.1021148, 0.1027496, 0.1033664, 0.1039700], atol=1e-6
    )


@pytest.mark.parametrize(
    "sensor,expected_nbands",
    [("ALI", 9), ("Hyperion", 242), ("Landsat4", 6), ("Landsat8", 8), ("MODIS", 19),
     ("Quickbird", 4), ("RapidEye", 5), ("Sentinel2a", 13), ("WorldView2-8", 8)],
)
def test_spectral_convolution_gaussian_sensor_characteristics(sensor, expected_nbands):
    wave, values = _toy_spectrum()
    res = spectral_convolution_gaussian(wave, values, sensor=sensor)
    assert len(res.wl) == expected_nbands
    # Every band is Gaussian-truncated to its own published [lb, ub] edges
    # -- bands entirely outside the toy spectrum's 400-2500nm range (e.g.
    # Hyperion's tail bands) are the only expected NaNs.
    center, lb, ub = sensor_characteristics(sensor)
    in_range = (ub >= 400) & (lb <= 2500)
    assert not np.isnan(res.rfl[np.argsort(center)][in_range[np.argsort(center)]]).any()


def test_spectral_convolution_gaussian_unknown_sensor_raises():
    wave, values = _toy_spectrum()
    with pytest.raises(ValueError, match="Unknown sensor"):
        spectral_convolution_gaussian(wave, values, sensor="Sentinel3A")


def test_spectral_convolution_gaussian_own_sensor_explicit_fwhm():
    # The AEO-Course "3-camera 15-band synchronized rig" worked example --
    # cross-checked against ToolsRTM::get.spectral.convolution.gaussian().
    wave, values = _toy_spectrum()
    centers = [444, 475, 502, 531, 550, 560, 570, 650, 668, 678, 705, 717, 740, 754, 842]
    fwhm = [28, 32, 18, 14, 12, 27, 14, 16, 14, 14, 10, 12, 18, 10, 57]
    res = spectral_convolution_gaussian(wave, values, centers=centers, fwhm=fwhm)

    assert len(res.wl) == 15
    assert not np.isnan(res.rfl).any()
    np.testing.assert_allclose(res.wl, centers)
    np.testing.assert_allclose(res.fwhm, fwhm)
    np.testing.assert_allclose(res.rfl[0], 0.103281, atol=1e-6)
    np.testing.assert_allclose(res.rfl[-1], 0.15104, atol=1e-5)


def test_spectral_convolution_gaussian_own_sensor_centers_only():
    # Headwall-style ENVI header: only band centers known, FWHM derived
    # from spacing (the contiguous-pushbroom-spectrometer assumption).
    wave, values = _toy_spectrum()
    centers = [398.0166, 400.247, 402.4774, 404.7078, 406.9382, 409.1686, 411.399, 413.6294, 415.8598, 418.0902]
    res = spectral_convolution_gaussian(wave, values, centers=centers)

    assert len(res.wl) == 10
    np.testing.assert_allclose(res.fwhm, np.full(10, 2.2304), atol=1e-3)
    assert not np.isnan(res.rfl).any()
