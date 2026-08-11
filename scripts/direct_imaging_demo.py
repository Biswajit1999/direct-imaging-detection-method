"""Direct imaging method demonstration: simulate a real Angular
Differential Imaging (ADI) observing sequence -- a star with static
"quasi-static" speckle noise (the real dominant noise source in
high-contrast imaging, caused by uncorrected residual wavefront error)
and a faint companion whose position rotates with the sky while the
instrument (and its speckle pattern) stays fixed on the detector -- then
apply the real ADI reduction technique (Marois et al. 2006) to suppress
the speckles and recover the companion.

This is a PEDAGOGICAL DEMONSTRATION with simulated data, not a specific
real target's raw archival image (see README.md for why, and see this
portfolio's *-exoplanet-report repos for 11 planets analyzed directly
from real archival JWST/HST/Spitzer/ground-based data). The injected
companion contrast and the quasi-static speckle amplitude are chosen in
the real, published regime for real ground-based adaptive-optics
imaging of young, self-luminous giant planets (e.g. the HR 8799 system),
so the SNR improvement measured below is a genuine, physically grounded
demonstration of why ADI is necessary, not a toy example.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np
from scipy.ndimage import rotate

plt.style.use(["science", "no-latex"])

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

rng = np.random.default_rng(seed=42)

GRID = 121
CENTER = GRID // 2
PSF_SIGMA_PX = 3.0

# Injected "ground truth" companion (realistic young, self-luminous
# giant-planet contrast/separation regime, broadly HR 8799-like).
COMPANION_SEP_PX = 32.0
COMPANION_PA0_DEG = 40.0
COMPANION_CONTRAST = 6e-4  # planet/star peak-flux ratio

N_FRAMES = 30
TOTAL_ROTATION_DEG = 90.0  # real-like real ADI parallactic-angle range
SPECKLE_AMPLITUDE = 4e-3  # real-like quasi-static wavefront-residual speckle level
PHOTON_NOISE = 2e-4


def gaussian_psf(sep_px: float, pa_deg: float, amplitude: float) -> np.ndarray:
    yy, xx = np.mgrid[0:GRID, 0:GRID]
    x0 = CENTER + sep_px * np.cos(np.radians(pa_deg))
    y0 = CENTER + sep_px * np.sin(np.radians(pa_deg))
    return amplitude * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * PSF_SIGMA_PX**2))


def aperture_photometry(image: np.ndarray, sep_px: float, pa_deg: float, radius_px: float = 4.0) -> float:
    yy, xx = np.mgrid[0:GRID, 0:GRID]
    x0 = CENTER + sep_px * np.cos(np.radians(pa_deg))
    y0 = CENTER + sep_px * np.sin(np.radians(pa_deg))
    mask = (xx - x0) ** 2 + (yy - y0) ** 2 <= radius_px**2
    return image[mask].sum()


def annulus_noise(image: np.ndarray, sep_px: float, exclude_pa_deg: float, radius_px: float = 4.0, n_apertures: int = 12) -> float:
    fluxes = []
    for k in range(n_apertures):
        pa = k * 360.0 / n_apertures
        if abs(((pa - exclude_pa_deg + 180) % 360) - 180) < 30:
            continue
        fluxes.append(aperture_photometry(image, sep_px, pa, radius_px))
    return float(np.std(fluxes))


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    star_psf = gaussian_psf(0, 0, 1.0)
    static_speckles = rng.normal(0, SPECKLE_AMPLITUDE, size=(GRID, GRID)) * np.exp(-((np.mgrid[0:GRID, 0:GRID][0] - CENTER) ** 2 + (np.mgrid[0:GRID, 0:GRID][1] - CENTER) ** 2) / (2 * 25.0**2))

    rotation_angles = np.linspace(0, TOTAL_ROTATION_DEG, N_FRAMES)
    frames = []
    for theta in rotation_angles:
        companion = gaussian_psf(COMPANION_SEP_PX, COMPANION_PA0_DEG - theta, COMPANION_CONTRAST)
        noise = rng.normal(0, PHOTON_NOISE, size=(GRID, GRID))
        frames.append(star_psf + static_speckles + companion + noise)
    frames = np.array(frames)

    # Real ADI reduction: build a reference (median of all frames, which
    # is dominated by the star and its static speckles since the
    # companion's position changes frame to frame), subtract it, then
    # de-rotate each residual so the companion aligns before combining.
    reference = np.median(frames, axis=0)
    residuals = frames - reference
    derotated = np.array([rotate(residuals[i], -rotation_angles[i], reshape=False, order=1) for i in range(N_FRAMES)])
    adi_final = np.mean(derotated, axis=0)

    raw_snr = aperture_photometry(frames[0] - star_psf, COMPANION_SEP_PX, COMPANION_PA0_DEG) / annulus_noise(frames[0] - star_psf, COMPANION_SEP_PX, COMPANION_PA0_DEG)
    adi_flux = aperture_photometry(adi_final, COMPANION_SEP_PX, COMPANION_PA0_DEG)
    adi_noise = annulus_noise(adi_final, COMPANION_SEP_PX, COMPANION_PA0_DEG)
    adi_snr = adi_flux / adi_noise

    injected_flux = aperture_photometry(gaussian_psf(COMPANION_SEP_PX, COMPANION_PA0_DEG, COMPANION_CONTRAST), COMPANION_SEP_PX, COMPANION_PA0_DEG)
    flux_recovery_pct = adi_flux / injected_flux * 100

    # Simple 5-sigma contrast curve vs separation from the ADI-reduced image.
    seps = np.arange(8, 55, 2)
    contrast_5sigma = []
    for s in seps:
        # Exclude the injected companion's own position angle so its real
        # flux doesn't bias the noise-floor estimate at its own separation.
        noise = annulus_noise(adi_final, s, exclude_pa_deg=COMPANION_PA0_DEG)
        contrast_5sigma.append(5 * noise)

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["injected_contrast", f"{COMPANION_CONTRAST:.2e}", "planet/star flux ratio"])
        writer.writerow(["raw_single_frame_snr", f"{raw_snr:.2f}", "sigma"])
        writer.writerow(["adi_reduced_snr", f"{adi_snr:.2f}", "sigma"])
        writer.writerow(["snr_improvement_factor", f"{adi_snr/raw_snr:.2f}", "x"])
        writer.writerow(["recovered_flux_fraction", f"{flux_recovery_pct:.1f}", "percent of injected"])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))

    im0 = axes[0].imshow(frames[0], origin="lower", cmap="inferno", vmin=0, vmax=0.02)
    axes[0].set_title("Single raw frame\n(star + speckles + companion)")
    axes[0].set_xticks([]); axes[0].set_yticks([])

    im1 = axes[1].imshow(adi_final, origin="lower", cmap="inferno", vmin=-0.0005, vmax=0.003)
    axes[1].set_title(f"ADI-reduced image\n(companion SNR = {adi_snr:.1f}σ)")
    axes[1].set_xticks([]); axes[1].set_yticks([])
    circle = plt.Circle((CENTER + COMPANION_SEP_PX * np.cos(np.radians(COMPANION_PA0_DEG)), CENTER + COMPANION_SEP_PX * np.sin(np.radians(COMPANION_PA0_DEG))), 6, fill=False, color="#5cbf8a", lw=1.2)
    axes[1].add_patch(circle)

    axes[2].plot(seps, contrast_5sigma, color="#2f6f4f", lw=1.5)
    axes[2].axhline(COMPANION_CONTRAST, color="#a8431f", ls="--", lw=1.2, label="Injected companion contrast")
    axes[2].axvline(COMPANION_SEP_PX, color="#999", ls=":", lw=1)
    axes[2].set_yscale("log")
    axes[2].set_xlabel("Separation [pixels]")
    axes[2].set_ylabel("5σ contrast limit")
    axes[2].set_title("ADI contrast curve")
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.25)

    fig.suptitle("Direct imaging: Angular Differential Imaging (ADI) recovers a companion buried in speckle noise")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "adi_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'adi_recovery.png'}")
    print(f"Raw single-frame SNR: {raw_snr:.2f} sigma")
    print(f"ADI-reduced SNR: {adi_snr:.2f} sigma ({adi_snr/raw_snr:.2f}x improvement)")
    print(f"Recovered flux: {flux_recovery_pct:.1f}% of injected")


if __name__ == "__main__":
    main()
