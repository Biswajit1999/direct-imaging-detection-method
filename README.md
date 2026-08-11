# Direct Imaging — Exoplanet Detection Method

The only method that produces an actual picture of a planet: resolve
its light as a separate point source next to its host star. This repo
explains the physics and implements a real Angular Differential Imaging
(ADI) reduction pipeline in Python from scratch, validated by injecting
a known companion buried in realistic speckle noise and recovering it.

## The physics

Direct imaging faces two combined, brutal real challenges: **contrast**
(a Jupiter-mass planet is roughly $10^{-4}$ to $10^{-9}$ times fainter
than its star, depending on age and wavelength) and **angular
separation** (most known planets are far too close to their star to be
resolved from typical distances). Adaptive optics correct atmospheric
turbulence in real time, but residual, slowly-evolving "quasi-static"
speckles — imperfections in the optics themselves — remain the real
dominant noise source at small separations, often far exceeding
ordinary photon noise. Angular Differential Imaging (Marois et al.
2006) exploits the fact that on an alt-az telescope with the field
derotator off, the sky (and any real companion) appears to rotate
around the star across an observing sequence while the instrument's own
speckle pattern stays fixed on the detector. Subtracting a
reference frame built from the sequence removes most of the
star and its speckles while a real companion — whose position moves
between frames — survives, and can be recovered by de-rotating each
residual frame back to sky orientation before combining.

## Why this method matters

Direct imaging is the only real technique that yields a planet's own
light directly: its spectrum, its temperature, sometimes even
resolved orbital motion over years, all without needing a fortunate
transit alignment or a large reflex velocity. It works best for young,
still-warm, self-luminous giant planets on wide orbits — real systems
like HR 8799 (four real imaged giant planets) and Beta Pictoris b
were found and characterized this way.

**Real limitation:** direct imaging is strongly biased toward young
(hot, still glowing from formation), massive, wide-separation planets
around nearby stars — it is currently essentially blind to older,
cooler, close-in planets like most of the real archival JWST/HST
targets covered elsewhere in this portfolio, which is exactly why
transit and radial-velocity spectroscopy remain necessary for
characterizing the bulk of the known exoplanet population.

## What this repo's code does

`scripts/direct_imaging_demo.py`:

1. Simulates a 30-frame ADI sequence: a fixed stellar point-spread
   function plus a **quasi-static speckle pattern** at a real-like
   amplitude, with a faint injected companion at a real young-giant-
   planet-like contrast ($6\times10^{-4}$) whose position angle rotates
   with the sky across 90 degrees of real-like parallactic-angle
   coverage while the speckles stay fixed on the detector.
2. Builds a reference PSF from the median of all frames, subtracts it
   from each frame, de-rotates each residual to align the sky, and
   combines them — the real ADI algorithm (Marois et al. 2006).
3. Measures the companion's detection significance via aperture
   photometry against an annulus noise estimate, in a single raw frame
   versus the final ADI-reduced image, and computes a 5-sigma contrast
   curve versus separation.

Run it yourself:

```bash
pip install -r requirements.txt
python scripts/direct_imaging_demo.py
```

## Result

| Quantity | Value |
|---|---|
| Injected contrast | 6.0×10⁻⁴ |
| Raw single-frame SNR | 0.73σ — **not detectable** |
| ADI-reduced SNR | 38.8σ — **clear detection** |
| SNR improvement | 53.4x |
| Recovered flux | 88.3% of injected |

In a single raw frame the companion is completely buried in speckle
noise (0.73σ — indistinguishable from a random fluctuation). After ADI
reduction it becomes an unambiguous, isolated point source at 38.8σ —
a genuine, quantified demonstration of why this technique, not
brute-force integration time, is what makes real direct imaging of
exoplanets possible.

## Honest limitation

The 5-sigma contrast curve shows a real, well-documented ADI artifact:
elevated noise right around the companion's own separation, caused by
**self-subtraction** — with only 90 degrees of real parallactic-angle
rotation, the companion's own signal partially contaminates the median
reference frame and leaves negative "side lobes" near its true
position (Milli et al. 2012). This is not a bug; it is a genuine,
physically real limitation of ADI with limited field rotation, and real
observing sequences are often planned specifically to maximize
parallactic-angle coverage to reduce it.

## Why this repo uses simulated (not raw archival) data

This repo demonstrates the *method itself* — how ADI turns an
undetectable signal into a clear one, and its own real limitations —
which is best shown with a known "ground truth" to validate recovery
against. This portfolio's companion `*-exoplanet-report` repositories
instead each analyze one real target's actual archival JWST/HST/
Spitzer/ground-based spectra directly, with zero simulated data. Both
approaches are stated plainly here rather than blurring the two.

## Repository structure

```text
scripts/direct_imaging_demo.py   ADI simulation + reduction pipeline + injection-recovery test
figures/                         generated plot + summary_statistics.csv
```

## References

1. Marois, C. et al., 2006. Angular Differential Imaging: A Powerful
   High-Contrast Imaging Technique. *The Astrophysical Journal*,
   641(1), pp.556-564.
2. Marois, C. et al., 2008. Direct Imaging of Multiple Planets Orbiting
   the Star HR 8799. *Science*, 322(5906), pp.1348-1352.
3. Chauvin, G. et al., 2004. A giant planet candidate near a young
   brown dwarf. *Astronomy & Astrophysics*, 425(2), L29-L32.
4. Milli, J. et al., 2012. Impact of angular differential imaging on
   circumstellar disk images. *Astronomy & Astrophysics*, 545, A111 —
   real self-subtraction bias in ADI.
5. Mawet, D. et al., 2014. Fundamental Limitations of High Contrast
   Imaging Set by Small Sample Statistics. *The Astrophysical
   Journal*, 792(2), 97.
6. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
