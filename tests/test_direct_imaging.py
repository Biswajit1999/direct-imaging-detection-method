"""Executable checks on the PSF/aperture photometry and the Mawet et al.
(2014) small-sample correction, including a regression guard on the
separation-dependent independent-aperture count (the earlier version of
this repo used a fixed count of 12 regardless of separation)."""

import numpy as np
import direct_imaging_demo as di


def test_gaussian_psf_total_flux_matches_analytic_integral():
    psf = di.gaussian_psf(0, 0, 1.0)
    analytic_total = 1.0 * 2 * np.pi * di.PSF_SIGMA_PX**2
    assert np.isclose(psf.sum(), analytic_total, rtol=1e-6)


def test_aperture_photometry_matches_analytic_enclosed_fraction():
    # For a 2D Gaussian, the fraction of total flux inside radius r is
    # 1 - exp(-r^2 / (2*sigma^2)) -- an independent closed-form check on
    # the pixel-sum aperture photometry, not just an internal consistency
    # check against the module's own code.
    psf = di.gaussian_psf(0, 0, 1.0)
    total_flux = psf.sum()
    radius_px = 4.0
    ap_flux = di.aperture_photometry(psf, 0, 0, radius_px=radius_px)
    analytic_fraction = 1 - np.exp(-(radius_px**2) / (2 * di.PSF_SIGMA_PX**2))
    assert abs(ap_flux / total_flux - analytic_fraction) < 0.02


def test_independent_apertures_scale_with_separation():
    # Regression guard: independent-aperture count must grow with
    # separation (circumference / resolution element), not stay fixed.
    counts = [di.independent_apertures_at_separation(sep) for sep in (8, 16, 32, 50)]
    assert counts == sorted(counts)
    assert counts[0] < counts[-1]
    for sep, n in zip((8, 16, 32, 50), counts):
        expected = max(3, int(np.floor(2 * np.pi * sep / (2 * 4.0))))
        assert n == expected


def test_small_sample_factor_converges_to_gaussian_for_large_n():
    # As n_independent -> infinity, the Student-t small-sample correction
    # must converge to the plain Gaussian sigma_level itself.
    factor = di.small_sample_sigma_factor(5.0, 100000)
    assert abs(factor - 5.0) < 0.01


def test_small_sample_factor_penalizes_small_n():
    # For very few independent samples, the correction must be
    # substantially larger than the naive Gaussian threshold -- this is
    # the entire point of the Mawet et al. (2014) correction.
    factor_small_n = di.small_sample_sigma_factor(5.0, 3)
    factor_large_n = di.small_sample_sigma_factor(5.0, 1000)
    assert factor_small_n > 10 * factor_large_n


def test_annulus_noise_excludes_companion_position_angle():
    rng = np.random.default_rng(0)
    image = rng.normal(0, 1e-4, size=(di.GRID, di.GRID))
    companion_pa = 40.0
    # Inject a large fake signal exactly at the excluded position angle;
    # it must not leak into the noise estimate.
    yy, xx = np.mgrid[0 : di.GRID, 0 : di.GRID]
    x0 = di.CENTER + 32.0 * np.cos(np.radians(companion_pa))
    y0 = di.CENTER + 32.0 * np.sin(np.radians(companion_pa))
    image[(yy - y0) ** 2 + (xx - x0) ** 2 <= 16] += 10.0

    noise_excluding, _ = di.annulus_noise(image, 32.0, exclude_pa_deg=companion_pa)
    noise_including, _ = di.annulus_noise(image, 32.0, exclude_pa_deg=999.0)
    assert noise_excluding < noise_including
