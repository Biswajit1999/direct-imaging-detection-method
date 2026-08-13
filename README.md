# Direct Imaging — Exoplanet Detection Method

The only method that produces an actual picture of a planet: resolve
its light as a separate point source next to its host star. This repo
works through the physics, implements an Angular Differential Imaging
(ADI) reduction pipeline in Python from scratch, and validates it by
injecting a companion into a simulated speckle-limited sequence and
recovering it — including a contrast curve calibrated in proper flux-
ratio units with a correction for how few independent noise samples
exist close to the star.

## The physics

### Two problems at once

Direct imaging faces two combined challenges: **contrast** (a Jupiter-
mass planet is roughly $10^{-4}$ to $10^{-9}$ times fainter than its
star, depending on age and wavelength) and **angular separation** (most
known planets sit far too close to their star to be resolved from
typical distances — the angular separation of a planet at 1 AU seen
from 10 parsecs is only 0.1 arcseconds, comparable to or smaller than
the diffraction limit of even a large telescope at visible wavelengths).
Both have to be overcome simultaneously: a wide-orbit planet is easier
to separate from its star angularly but usually fainter in absolute
terms too, since it intercepts less starlight and formed farther from
the heat of formation.

### Why speckles, not just diffraction, set the floor

A telescope's diffraction pattern (the Airy pattern for a circular
aperture) sets a hard geometric limit on angular resolution, but real
high-contrast imaging is usually limited by something else entirely:
adaptive optics correct atmospheric turbulence in real time, but leave
behind residual, slowly-evolving "quasi-static" speckles — small
imperfections in the optics and imperfect AO correction that produce a
mottled pattern of bright and dark spots across the field, each one
looking exactly like a faint point source would. These speckles evolve
on minutes-to-hours timescales (as temperature and mechanical flexure
shift the optics slightly) and typically dominate over ordinary photon
noise at the separations where planets are found, which is why simply
integrating longer doesn't help much — the noise doesn't average down
the way photon noise does, because it isn't random from exposure to
exposure.

### How ADI tells a planet from a speckle

Angular Differential Imaging (Marois et al. 2006) exploits geometry: on
an alt-az telescope with the image derotator switched off (or
deliberately disabled), the sky appears to rotate around the field
center over the course of a night while the telescope's own optics —
and their speckle pattern — stay fixed relative to the detector. A real
companion, fixed on the sky, traces an arc across the detector as the
field rotates; the speckles don't move at all. Building a reference
image from the sequence (the simplest version: the median of all
frames) captures the star and its speckles well, since they're present
in every frame in the same place, while the companion is smeared across
many different positions and contributes little to the median. Subtract
that reference from each frame, and what's left is mostly the moving
companion signal on a much fainter noise floor. De-rotating each
residual frame back to a common sky orientation before combining makes
the companion's signal add up coherently while the leftover speckle
residue, no longer at a fixed detector position, partially cancels.

## Why this method matters

Direct imaging is the only technique that yields a planet's own light
directly: its spectrum, its temperature, sometimes resolved orbital
motion over years, all without needing a fortunate transit alignment or
a large reflex velocity. It works best for young, still-warm, self-
luminous giant planets on wide orbits — systems like HR 8799 (four
imaged giant planets) and Beta Pictoris b were found and characterized
this way, and it remains the main way to study a planet's atmosphere
independent of transmission or emission spectroscopy during transit.

**Limitation:** direct imaging is strongly biased toward young (hot,
still glowing from formation), massive, wide-separation planets around
nearby stars — it's currently close to blind to older, cooler, close-in
planets like most of the archival JWST/HST targets covered elsewhere in
this portfolio, which is exactly why transit and radial-velocity
spectroscopy remain necessary for characterizing the bulk of the known
exoplanet population.

## What this repo's code does

`scripts/direct_imaging_demo.py`:

1. Simulates a 30-frame ADI sequence: a fixed stellar point-spread
   function plus a quasi-static speckle pattern at a realistic
   amplitude, with a faint injected companion at a young-giant-
   planet-like contrast ($6\times10^{-4}$) whose position angle rotates
   with the sky across 90 degrees of parallactic-angle coverage while
   the speckles stay fixed on the detector.
2. Builds a reference PSF from the median of all frames, subtracts it
   from each frame, de-rotates each residual to align the sky, and
   combines them — the ADI algorithm from Marois et al. (2006).
3. Measures the companion's detection significance via aperture
   photometry against an annulus noise estimate, comparing a single raw
   frame against the final ADI-reduced image.
4. Computes a 5-sigma contrast curve versus separation, normalized to
   the star's own aperture flux (so it's a genuine dimensionless
   contrast, comparable directly to the injected value) and corrected
   for the small number of independent noise samples available in a
   narrow annulus close to the star (Mawet et al. 2014) rather than
   assuming a flat Gaussian threshold everywhere.

Run it yourself:

```bash
pip install -r requirements.txt
python scripts/direct_imaging_demo.py
```

## Worked example with a real target

HR 8799 b, one of the first planets ever directly imaged (Marois et
al. 2008), separates the two challenges from "The physics" above
cleanly. Its real measured angular separation is 1.713 arcseconds at a
distance of 39.4 parsecs — a physical separation of about 67 AU. The
discovery used Keck (10 m primary mirror) and Gemini North (8 m) in
the near-infrared H band (1.6 microns); the diffraction limit for a
10 m telescope there is:

```
theta = 1.22 * lambda / D
      = 1.22 * 1.6e-6 m / 10 m
      = 1.95e-7 rad = 0.040 arcsec
```

The real separation (1.71 arcsec) is roughly 40 times larger than that
diffraction limit — for this particular planet, at this separation,
resolution was never the bottleneck. What made the detection hard
was contrast: HR 8799 b is roughly $10^{-5}$ times fainter than its
star in the near-infrared, well into the regime where quasi-static
speckle noise, not the diffraction limit, sets the real detection
floor, which is exactly the problem ADI (used in this repo's
simulation) was built to solve.

## Result

| Quantity | Value |
|---|---|
| Injected contrast | 6.0×10⁻⁴ |
| Raw single-frame SNR | 0.73σ — not detectable |
| ADI-reduced SNR | 38.8σ — clear detection |
| SNR improvement | 53.4x |
| Recovered flux | 88.3% of injected |
| 5σ contrast limit at the companion's separation | 1.9×10⁻⁴ (10 independent apertures) |

In a single raw frame the companion is buried in speckle noise (0.73σ —
indistinguishable from a random fluctuation). After ADI reduction it
becomes an isolated point source at 38.8σ, and the injected contrast
sits comfortably above the 5σ detection limit everywhere on the
contrast curve — a quantified demonstration of why this technique, not
longer integration time on its own, is what makes direct imaging of
exoplanets practical.

## Limitations

The contrast curve shows a documented ADI artifact: elevated noise
right around the companion's own separation, caused by
**self-subtraction** — with only 90 degrees of parallactic-angle
rotation, the companion's own signal partially contaminates the median
reference frame and leaves negative "side lobes" near its true
position (Milli et al. 2012). That's not a bug; it's a real limitation
of ADI with limited field rotation, and real observing sequences are
often planned specifically to maximize parallactic-angle coverage to
reduce it. Separately, the quasi-static speckle field here is a static
Gaussian-random texture rather than a correlated, slowly evolving PSF
residual, and algorithmic throughput (how much real companion flux ADI
itself removes through self-subtraction, beyond what this repo's
88.3%-recovery number already shows) isn't independently calibrated via
fake-planet injection at multiple separations, which a published
contrast curve would do.

## Extending this

A natural next step: inject fake companions at a grid of separations
and position angles, run the same pipeline on each, and use the
fraction of flux recovered at each point as an empirical throughput
correction for the contrast curve — this is standard practice in
real high-contrast imaging pipelines and would tighten the gap between
this repo's simplified curve and a publication-grade one. You could
also replace the median-combination reference PSF with a more capable
algorithm like KLIP (Karhunen-Loève Image Projection, Soummer et al.
2012) or LOCI (Lafrenière et al. 2007), both of which build a smarter
reference from a weighted combination of the other frames and typically
recover more companion flux at fixed self-subtraction — implemented in
real pipelines such as `pyKLIP`.

## Why this repo uses simulated (not raw archival) data

This repo demonstrates the *method itself* — how ADI turns an
undetectable signal into a clear one, and where its own approximations
break down — which is best shown with a known "ground truth" to
validate recovery against. This portfolio's companion `*-exoplanet-
report` repositories instead each analyze one real target's archival
JWST/HST/Spitzer/ground-based spectra directly, with no simulated data.
Both approaches are stated plainly here rather than blurring the two.

## Repository structure

```text
scripts/direct_imaging_demo.py   ADI simulation + reduction pipeline + calibrated contrast curve
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
   self-subtraction bias in ADI.
5. Mawet, D. et al., 2014. Fundamental Limitations of High Contrast
   Imaging Set by Small Sample Statistics. *The Astrophysical
   Journal*, 792(2), 97 — the small-sample correction applied above.
6. Soummer, R., Pueyo, L. and Larkin, J., 2012. Detection and
   Characterization of Exoplanets and Disks Using Projections on
   Karhunen-Loeve Eigenimages. *The Astrophysical Journal Letters*,
   755(2), L28.
7. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
