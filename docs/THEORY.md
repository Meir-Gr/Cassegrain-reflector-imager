# Physical and Mathematical Foundations of the Cassegrain Defocus Simulations

This document explains the theory behind every computation in this project
(cassegrain_defocus_sim.py, cassegrain_extended_source_sim.py,
cassegrain_parameter_study.py, cassegrain_vibration_sim.py,
cassegrain_arago_vs_source.py, cassegrain_gui.py and its MATLAB port
cassegrain_simulator.m), with emphasis on the
**assumptions**, the **approximations**, and the **regimes where the model
stops being valid** and would need additional physics.

---

## 1. The model chain

Every image produced by the suite is computed as

$$I(\mathbf{u}) = \left[\, |\mathcal{F}\{P\}|^2 \ast S_{source} \right](\mathbf{u}) \,\times\, M_{blur}(\mathbf{u})$$

with the complex pupil function

$$P(x,y) = A(x,y)\, e^{\,i\,2\pi\,[\,W_{020}\rho^2 \;(+\; W_{040}\rho^4)\,]},\qquad \rho = r/R .$$

- $A(x,y)$: binary aperture transmission — annulus (primary minus central
  obstruction) minus three spider-vane strips.
- $W_{020}$: defocus in **waves at the pupil edge**; $W_{040}$: primary
  spherical aberration (optional).
- $\mathcal{F}$: 2-D Fourier transform, realized as an FFT.
- $S_{source}$: geometric image of the source (uniform disk).
- $M_{blur}$: vibration and drift MTF factors, applied multiplicatively in
  the frequency domain.

In the interactive simulators the point-source PSF is additionally an
incoherent sum over a flat spectral band of fractional width
$\Delta\lambda/\lambda$ (§7.3); the batch study scripts are monochromatic.

Each stage rests on specific physical assumptions, detailed below.

---

## 2. Scalar Fourier optics

### 2.1 What is assumed

| Assumption | Statement | Where it breaks |
| --- | --- | --- |
| Scalar diffraction | Light treated as one complex scalar field; polarization ignored | Fast beams (roughly $F\# \lesssim 1.5$, NA $\gtrsim 0.3$), polarization-sensitive instruments, vector effects at edges |
| Paraxial / Fraunhofer imaging | The exit-pupil field maps to the focal plane by a Fourier transform | Very fast beams; systems with large pupil aberrations of high order; near-field imaging |
| Thin-screen (Kirchhoff) aperture | The aperture multiplies the field by 0/1 at a single plane | Thick baffles, edge scattering, conductive-edge effects at scales $\sim\lambda$ |
| Monochromatic light | Single wavelength $\lambda$ (batch scripts); interactive simulators add flat-band polychromatic light (§7.3) | weighted spectra, chromatic optics, atmospheric dispersion (§7.3) |
| Point source at infinity, on axis | Incident field = unit plane wave | Wide field angles (needs field-dependent aberrations), finite-distance sources (extra quadratic phase) |

Under these assumptions the focal-plane field of a point source is the
Fraunhofer relation

$$U(\mathbf{u}) \,\propto\, \iint P(x,y)\, e^{-i \frac{2\pi}{\lambda f} (xu + yv)}\, dx\, dy, \qquad \mathrm{PSF} = |U|^2 .$$

The Fourier-transforming property of a converging reference sphere is exact
within Fresnel diffraction; corrections are of order $(\theta_{max})^2 = 1/(2F\#)^2$
— about $0.25\%$ at $f/10$, which is why an $f/10$ Cassegrain is a
near-ideal subject for this method.

### 2.2 Defocus as a quadratic pupil phase

Displacing the observation plane by $\Delta z$ adds an optical path
difference across the converging cone $\mathrm{OPD}(\theta) = \Delta z\,(1-\cos\theta) \approx \Delta z\,\theta^2/2$.
With the marginal-ray angle $\theta_{max} = 1/(2F\#)$, the edge OPD in waves is

$$W_{020} = \frac{\Delta z}{8\,\lambda\,F\#^2} \qquad\Longleftrightarrow\qquad \Delta z = 8\,\lambda\,F\#^2\,W_{020} .$$

**Approximations involved:** the $\cos\theta$ expansion is truncated at
$\theta^2$ (error $\sim \theta^2/12$ relative, i.e. $<0.1\%$ at $f/10$), and
the mapping is linearized around focus. For $F\# \lesssim 4$ or defocus
comparable to the focal length, the quartic term of $(1-\cos\theta)$ and
pupil-scaling effects must be added; the quadratic form here silently
overestimates edge OPD for very fast cones.

### 2.3 Geometric (large-defocus) limit

The stationary-phase evaluation of the diffraction integral for
$W_{020} \gg 1$ reproduces ray optics: each pupil point maps to the image
point $\mathbf{u} \propto \nabla W$. Two used consequences:

- **Donut outer diameter** $\varnothing = 8\,W_{020}\,\lambda F\# = \Delta z/F\#$; the image converges to a scaled shadow of the pupil (verified to ~3 %, the residual being the diffraction-softened edge).
- **Sign of defocus:** for a real aperture, $\mathrm{PSF}(-W) (\mathbf{u}) = \mathrm{PSF}(+W)(-\mathbf{u})$ — an exact point reflection (this is why the 3-fold spider pattern appears twisted by 60° across focus). This exactness holds only for *pure defocus with a real* $A$; any odd aberration (coma) breaks it.

---

## 3. The aperture model

$A$ is built as: disk of radius $R = D/2$, minus disk of radius
$\varepsilon R$ (obstruction), minus three half-infinite strips of width $w$
from the center outward at the given angles.

**Assumptions and their consequences:**

1. **All obscurations lie in one plane (the pupil).** Physically the vanes
   and secondary sit at different axial positions than the entrance pupil.
   Treating them as coplanar ignores Fresnel propagation between those
   planes. The error scales with the Fresnel number of the vane width over
   the propagation distance to the pupil image; for amateur/laboratory
   geometries it slightly softens vane-edge diffraction detail but leaves
   widths, depths and all scalings intact. High-fidelity stray-light or
   apodization studies require multi-plane (angular-spectrum) propagation.
2. **Perfectly absorbing, infinitely thin vanes.** No reflection or
   scattering off vane faces, no edge birefringence. Reflective vane faces
   in real telescopes add faint glints not modeled here.
3. **Binary sampled mask (no anti-aliasing).** Edges are staircase at the
   grid pitch $D/N_{pupil}$. This injects weak high-frequency artifacts far
   in the PSF wings — irrelevant at the intensity levels analysed here
   ($\gtrsim 10^{-6}$ of peak), but a limitation for coronagraph-grade
   contrast studies ($<10^{-8}$), which need gray-pixel edges or analytic
   aperture transforms.
4. **No wavefront error map.** The optics are perfect except for the
   explicit $W_{020}, W_{040}$ terms. Real mirrors add zonal/figure errors;
   these are straightforward to include by adding a phase map
   $e^{i2\pi W_{err}(x,y)}$ but change the star-test appearance
   qualitatively when $W_{err} \gtrsim \lambda/8$.

---

## 4. Numerical realization (DFT) and its validity conditions

### 4.1 Sampling of the pupil phase

Between adjacent pupil samples the defocus phase changes by at most
$\delta\phi = 8\pi W_{020}/N_{pupil}$. The suite enforces
$N_{pupil} \geq 16\,|W_{020}|$ (i.e. $\delta\phi \leq \pi/2$), a factor-2
margin over the Nyquist limit $\delta\phi < \pi$. **If violated, the
quadratic phase aliases and spurious inner rings/ghost donuts appear.** This
is the single most important numerical constraint; it is why grid size grows
with defocus (adaptive tiers in the GUI, `n_pupil_for()` in the scripts).

### 4.2 Image-plane sampling

With zero-padding factor 2 the image pixel is $\lambda F\#/2$ — the Nyquist
spacing for the *coherent field* (cutoff $D/\lambda f$). The intensity has
twice that bandwidth, and with pad = 2 its spectrum extends exactly to the
Nyquist edge: the computed samples are exact, and frequency-domain
operations on the PSF (source/vibration MTFs) are alias-free *by margin
zero*. Sub-pixel interpolation of intensity (used in the fine radial
profiles) is therefore slightly band-edge-limited; the analysis code uses
cubic spline sampling at $\leq 0.25$ px steps to stay safe.

### 4.3 Periodic wrap-around

The DFT is periodic: light diffracted beyond the computational field of view
re-enters from the opposite side. The crop guards keep the geometric blur
below ~75 % of the half-field, so wrap-around only affects the faint spike
tails (below the $10^{-4}$ level in the displayed stretches).

### 4.4 Finite precision

Arrays are float32/complex64. Reliable intensity dynamic range is roughly
$10^{-7}$ of the peak; the $10^{-6}$ floor of the log displays is safe, but
contrast statements below $\sim10^{-7}$ would require float64.

---

## 5. Extended sources

### 5.1 Incoherent superposition

The image of an extended source is computed as

$$I = \mathrm{PSF} \ast S, \qquad \hat{I}(\mathbf{q}) = \mathrm{OTF}(\mathbf{q})\cdot \hat{S}(\mathbf{q}),$$

with, for a uniform disk of angular diameter $\theta_d$ ($q$ in
cycles/radian),

$$\hat{S}(q) = \frac{2\, J_1(\pi\, \theta_d\, q)}{\pi\, \theta_d\, q} .$$

This rests on **three assumptions**:

1. **Spatial incoherence** of the source (each source point radiates
   independently). True for thermal/astronomical sources by the van
   Cittert–Zernike theorem at telescope scales. **Fails for laser-illuminated
   pinholes or any partially coherent illumination** — then intensities do
   not add; one must propagate the mutual coherence function or sum coherent
   fields, and the convolution model (and all conclusions built on it) no
   longer applies.
2. **Isoplanatism**: the PSF is identical for every point of the source.
   Valid for sources of arcseconds against telescope field aberrations that
   change over arcminutes. For wide-field simulation, PSF must become
   field-dependent.
3. **Uniform, circular source.** Real planetary disks are limb-darkened
   (slightly weaker high-frequency sidelobes than $2J_1(x)/x$) and not
   necessarily circular; replace $\hat S$ accordingly if that matters.

### 5.2 Seeing

Long-exposure atmospheric seeing was treated as "an extra source diameter."
Strictly, the long-exposure seeing OTF is
$\exp[-3.44\,(\lambda q/r_0)^{5/3}]$ (Kolmogorov), which is close to, but not
exactly, a disk/Gaussian kernel. The qualitative feature-survival hierarchy
is unaffected; quantitative work in strong turbulence — and *any* statement
about short exposures (speckle) or scintillation — requires phase-screen
Monte Carlo, which is **not** modeled.

---

## 6. Vibration and motion blur

All time-dependent line-of-sight (LOS) motion is modeled in the
**fully-averaged limit**: exposure time $\gg$ motion time-scale, so the
image is the PSF convolved with the *time-marginal distribution* of the LOS
displacement. Each mechanism contributes a multiplicative MTF factor:

| Mechanism | Kernel (dwell distribution) | MTF factor |
| --- | --- | --- |
| Random jitter, RMS $\sigma$ (many uncorrelated cycles) | Gaussian | $\exp(-2\pi^2\sigma^2 q^2)$ |
| Sinusoidal resonance, amplitude $a$, along $\hat{e}$ | Arcsine $\frac{1}{\pi\sqrt{a^2-x^2}}$ | $J_0(2\pi a\, q_\parallel)$ |
| Linear drift, length $L$ per integration, along $\hat{e}$ | Top-hat of width $L$ | $\mathrm{sinc}(L\, q_\parallel)$ |

Derivation sketch for the sinusoid: $x(t) = a\sin\omega t$ spends
$dt \propto dx/\sqrt{a^2 - x^2}$ near displacement $x$; the Fourier
transform of that arcsine density is $J_0(2\pi a q)$ — oscillating in sign,
hence the ring-contrast nulls and reversed-contrast revivals verified in the
simulations.

**Assumptions and limits:**

1. **Rigid image motion (tilt only).** The vibration translates the PSF
   without deforming it. Mirror bending modes, defocus oscillations
   (despace vibration), or wavefront-error modulation need time-resolved
   wavefront modeling.
2. **Full averaging.** If the exposure covers only ~1–3 vibration periods,
   the effective kernel depends on the phase at shutter-open: images become
   asymmetric double-spots with random orientation jitter frame-to-frame.
   Model by Monte-Carlo averaging of partial-cycle kernels.
3. **Stationarity.** Amplitude/frequency assumed constant during the
   exposure.
4. **Multiplicativity of the MTF factors** (source × jitter × drift) is
   exact *given* linear shift-invariant incoherent imaging — it inherits
   the assumptions of §5.1, and was verified numerically to ~10 %
   (limited by the ring-contrast metric, not the physics).

---

## 7. Additional approximations worth knowing about

### 7.1 Normalization

Images are peak-normalized for display; no absolute photometry (source flux,
throughput, QE) is modeled. Energy bookkeeping (Parseval) is used only as an
internal sanity check.

### 7.2 The axial (Arago) spot

The bright on-axis feature inside the obstruction shadow emerges naturally
from the diffraction integral (it is not added by hand). Its intensity
oscillates rapidly with defocus, so single-wavelength, single-defocus values
of $I(0)$ quoted in the reports are exact but *fragile*: a few percent
bandwidth or a fraction of a wave of defocus shifts them substantially. Its
convolution with a resolved source (making a central image of the source) is
robust.

### 7.3 Polychromatic light (implemented in the interactive simulators)

The interactive simulators (cassegrain_gui.py, cassegrain_simulator.m) model
finite bandwidth; the batch study scripts remain monochromatic. The model is
a **flat (top-hat) spectrum** of fractional width $b = \Delta\lambda/\lambda_0$,
realized as an incoherent, equal-energy sum of $N_\lambda$ monochromatic PSFs:

$$I_{poly} = \frac{1}{N_\lambda}\sum_{k} I_{\lambda_k}, \qquad \lambda_k = \lambda_0\left(1 - \frac{b}{2} + b\,\frac{k - 1/2}{N_\lambda}\right),\quad k = 1..N_\lambda .$$

Implementation choices and their justification:

- **Common angular grid without interpolation.** For each sub-wavelength the
  *pupil* is redrawn on the fixed grid, scaled by $\lambda_0/\lambda_k$. This
  makes the angular image pixel exactly $\lambda$-independent, so sub-images
  add sample-by-sample with no resampling error.
- **Achromatic focus.** An all-mirror telescope imposes one physical $\Delta z$
  on all wavelengths, hence
  $W_{020}(\lambda_k) = W_{020}(\lambda_0)\,\lambda_0/\lambda_k$ — bluer light
  carries more waves of the same defocus.
- **Sampling invariance.** Under the pupil rescaling the phase step per
  sample, $8\pi W(\lambda)/N_{pupil}(\lambda)$, is $\lambda$-independent, so
  the anti-aliasing criterion of §4.1 is inherited unchanged for every
  sub-wavelength.
- **Adaptive $N_\lambda$ (3–15).** Chosen so adjacent sub-images shift by
  less than about one ring width at the donut rim; too few samples would
  produce spurious families of discrete rings instead of a smooth chromatic
  smear.

Physical consequences (visible in the simulators): fringes and rings smear
radially in proportion to $b$ times their distance from the achromatic
anchoring edges, while the donut diameter, the hole, and the vane-shadow
geometry are strictly achromatic ($\varnothing = \Delta z/F\#$ contains no
$\lambda$).

Remaining assumptions and limits of the polychromatic model:

- **Flat spectrum in energy** — no source color, filter transmission curve,
  or detector QE weighting (a weighted sum over $\lambda_k$ is a trivial
  extension).
- **All-reflective optics** — no chromatic focus shift. Catadioptric systems
  (corrector plates) add real longitudinal chromatism that must be modeled
  as a $\lambda$-dependent $\Delta z$.
- **No atmospheric dispersion** — the prism-like elongation of broadband
  images at low elevation is not modeled.
- **Blue-edge sampling margin** — for $b \gtrsim 0.2$ the bluest sub-pupils
  slightly exceed the pad-2 Nyquist margin of §4.2 (intensity spectrum
  extends beyond the band edge by up to $1/(1-b/2)$); the aliased content is
  far below the displayed contrast floor but would matter for precision
  far-wing photometry.

### 7.4 Where each headline result stops holding

| Result | Valid while… | Extension needed when… |
| --- | --- | --- |
| $\varnothing = 8 W_{020} \lambda F\#$ | $F\# \gtrsim 4$, $W_{020} \gtrsim 2$ | fast beams (exact OPD), near focus (use quadrature with diffraction floor) |
| 60° twist across focus | pure defocus, real aperture | any odd aberration present (coma) — the symmetry argument fails |
| Ring period ≈ 7 λ/D (interior) | this pupil geometry ($\varepsilon=0.35$) | different obstruction ratios shift the edge-beat frequency |
| Feature-survival hierarchy vs source size | incoherent, isoplanatic, uniform disk | coherent sources; limb-darkened/elongated sources (quantitative shifts) |
| Vibration MTF factors | fully averaged, rigid motion | short exposures, non-rigid (bending/despace) vibration |
| Fresnel floor $\sqrt{\lambda \Delta z}$ for vane shadows | shadow width ≪ donut, scalar edge diffraction | vanes far from pupil plane with large Fresnel numbers (multi-plane propagation) |

---

## 8. Explicit list of *unmodeled* physics

If your use case touches any of these, the present suite must be extended:

1. **Polarization / vector diffraction** (fast optics, coatings, polarimetry).
2. **Spectral weighting & chromatic optics** — flat-band integration *is*
   implemented in the interactive simulators (§7.3); still missing:
   source-color/QE-weighted spectra, the chromatic focus shift of
   catadioptric correctors, and atmospheric dispersion.
3. **Atmospheric turbulence beyond a blur kernel** — speckle, scintillation,
   anisoplanatism: needs Kolmogorov phase screens.
4. **Multi-plane obstruction diffraction** (vane/secondary axial positions,
   baffle edges): angular-spectrum propagation between planes.
5. **Stray light**: baffle leaks, ghost reflections, vane-face glints,
   scatter from dust/roughness — non-sequential ray tracing or dedicated
   coherent models (the central-spot *case study* modeled a ghost simply as
   an incoherent second defocused image, which captures geometry but not
   coating physics).
6. **Mirror figure errors and thermal gradients**: add measured/synthetic
   wavefront maps to the pupil phase.
7. **Detector physics**: pixel-aperture MTF, undersampling/aliasing on real
   pixels, noise (photon/read), nonlinearity, saturation, rolling shutter.
8. **Partially coherent illumination** (laboratory pinholes with lasers).
9. **Finite-conjugate imaging** (source not at infinity): extra quadratic
   phase and pupil-magnification bookkeeping.

---

## 9. Numerical cross-checks that anchor the implementation

The implementation was validated against closed-form theory wherever one
exists — these agreements are the evidence that the code computes what the
formulas above claim:

| Check | Theory | Simulation |
| --- | --- | --- |
| EE80 diameter vs defocus | $8 W \lambda F\#$ asymptote + diffraction floor | matches over 1–128 waves (log-log) |
| Hole/outer diameter ratio | $= \varepsilon$ exactly | $\varepsilon$ to a few % at all $\varepsilon, W$ tested |
| Ring-contrast decay vs source size | disk visibility $2J_1(x)/x$ at the *measured* ring frequency | 0.95/0.87/0.64 vs 0.98/0.91/0.67 |
| Sinusoidal-vibration null | first zero of $J_0$: $a = 2.405/(2\pi q_{ring})$ | null observed at the predicted amplitude; side-lobe revivals match $J_0$ lobes (0.36/0.31 vs 0.31/0.38) |
| Intra/extra symmetry | exact 180° point reflection | panel-by-panel identical (mirrored) images |
| MTF multiplicativity | product of factors | 0.57 measured vs 0.51 predicted (metric-limited) |

The residual discrepancies are dominated by measurement metrics (band-pass
of the detrending, broadband ring spectra), not by the diffraction engine.
