"""Independent synthetic ADI robustness checks (not observed planet measurements)."""
import numpy as np
from scipy.ndimage import rotate

def adi_reduce(cube, angles, exclusion_deg=0):
    """Independent-frame median reference and sky derotation."""
    cube = np.asarray(cube, float)
    angles = np.asarray(angles, float)
    if cube.ndim != 3 or cube.shape[0] != len(angles) or cube.shape[1] != cube.shape[2]:
        raise ValueError("cube must be N x P x P matching angles")
    if not np.all(np.isfinite(cube)) or not np.all(np.isfinite(angles)) or exclusion_deg < 0:
        raise ValueError("finite inputs and nonnegative exclusion required")
    residuals = []
    for i, angle in enumerate(angles):
        eligible = np.abs(angles-angle) > exclusion_deg
        if not np.any(eligible):
            raise ValueError("no independent reference frames")
        residuals.append(rotate(cube[i]-np.median(cube[eligible], axis=0),
                                -angle, reshape=False, order=1))
    return np.mean(residuals, axis=0)

def matched_filter_amplitude(image, template, variance):
    """Independent-pixel inverse-variance PSF amplitude and formal uncertainty."""
    y, t, v = (np.asarray(z, float) for z in (image, template, variance))
    if y.shape != t.shape or y.shape != v.shape or not all(np.all(np.isfinite(z)) for z in (y,t,v)) or np.any(v <= 0):
        raise ValueError("matching finite arrays and positive variance required")
    information = np.sum(t*t/v)
    if information <= 0:
        raise ValueError("zero template information")
    return float(np.sum(y*t/v)/information), float(information**-0.5)
