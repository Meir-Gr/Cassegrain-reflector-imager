"""
Out-of-focus point-source image (star test) for a Cassegrain reflector
whose secondary mirror is held by a 3-vane spider (triangular support,
vanes 120 degrees apart).

Physics (scalar Fourier optics)
-------------------------------
The image of a point source is the point spread function (PSF):

    PSF = | FFT{ P(x, y) } |^2
    P(x, y) = A(x, y) * exp[ i * 2*pi * W020 * (r/R)^2 ]

A(x, y)  : aperture transmission = annulus (primary minus secondary
           obstruction) minus the three spider vanes.
W020     : defocus coefficient in waves at the pupil edge.

A defocused star therefore shows the classic "donut" (shadow of the
secondary) crossed by the shadows of the three vanes, and the 3-vane
spider additionally produces 6 diffraction spikes (each vane creates a
pair of spikes perpendicular to itself; 3 vanes at 120 deg -> 6 spikes
every 60 deg).

Outputs
-------
cassegrain_pupil_and_psfs.png  : pupil + PSF gallery for small defocus values
cassegrain_star_test.png       : intra-focal vs extra-focal comparison
cassegrain_defocus_series.png  : through-focus morphology, intra- and extra-focal
cassegrain_spot_size.png       : spot diameter vs defocus (diffraction -> geometric)
cassegrain_star_test_spherical.png : through-focus asymmetry with spherical aberration
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Output folder for the figures of every study (others use optics.FIG_DIR).
FIG_DIR = Path(__file__).resolve().parent.parent / "docs" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------
# User parameters (edit freely)
# ----------------------------------------------------------------------
D_PRIMARY   = 0.200        # primary mirror diameter [m]
OBSTRUCTION = 0.35         # secondary obstruction ratio (d_sec / D)
VANE_WIDTH  = 0.004        # spider vane width [m]
VANE_ANGLES = np.deg2rad([90.0, 210.0, 330.0])  # 3 vanes, 120 deg apart
WAVELENGTH  = 550e-9       # observation wavelength [m]
F_NUMBER    = 10.0         # effective focal ratio (typical Cassegrain)

N_PUPIL     = 1024         # base grid samples across the pupil diameter
PAD_FACTOR  = 2            # FFT padding -> image pixel = lambda*F#/PAD (Nyquist)
DEFOCUS_SET = [0.0, 3.0, 8.0, 16.0]   # W020 [waves] for the small-defocus gallery
STAR_TEST_W = 10.0         # +/- defocus [waves] for intra/extra comparison
SPHERICAL_W040 = 0.5       # primary spherical aberration [waves] for figure 5
SA_TEST_SET = [2.0, 4.0, 8.0, 16.0]   # |defocus from best focus| [waves]

# Through-focus series: |W020| doubling steps, applied with both signs
# (negative = intra-focal, positive = extra-focal).
GALLERY_SET = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]
MAX_N_PUPIL = 4096         # memory guard for the adaptive pupil sampling


def n_pupil_for(w):
    """Pupil sampling so the defocus phase steps < ~pi/2 per sample."""
    n = N_PUPIL
    while n < 16.0 * abs(w):
        n *= 2
    if n > MAX_N_PUPIL:
        raise ValueError(f"Defocus {w} waves would need N_PUPIL={n}; "
                         f"raise MAX_N_PUPIL if you have the RAM.")
    return n


# ----------------------------------------------------------------------
# Pupil construction
# ----------------------------------------------------------------------
_MASK_CACHE = {}
_RHO2_CACHE = {}


def make_pupil(n_pupil=N_PUPIL, obstruction=None, vane_width=None):
    """Return (aperture mask, (r/R)^2 map) on the padded grid (cached).

    obstruction / vane_width override the module defaults (parameter studies).
    """
    eps = OBSTRUCTION if obstruction is None else obstruction
    wv = VANE_WIDTH if vane_width is None else vane_width
    key = (n_pupil, eps, wv)
    if key in _MASK_CACHE:
        return _MASK_CACHE[key], _RHO2_CACHE[n_pupil]
    n = n_pupil * PAD_FACTOR
    R = D_PRIMARY / 2.0
    dx = D_PRIMARY / n_pupil
    coords = ((np.arange(n) - n // 2) * dx).astype(np.float32)
    X, Y = np.meshgrid(coords, coords)
    rr = np.hypot(X, Y)

    mask = (rr <= R) & (rr >= eps * R)

    # Spider vanes: half-infinite strips from the centre outwards.
    for ang in VANE_ANGLES:
        ux, uy = np.cos(ang), np.sin(ang)
        along = X * ux + Y * uy          # coordinate along the vane
        across = -X * uy + Y * ux        # perpendicular distance
        mask &= ~((along >= 0.0) & (np.abs(across) <= wv / 2.0))

    if n_pupil not in _RHO2_CACHE:
        _RHO2_CACHE[n_pupil] = ((rr / R) ** 2).astype(np.float32)
    m = mask.astype(np.float32)
    _MASK_CACHE[key] = m
    return m, _RHO2_CACHE[n_pupil]


def compute_psf(mask, rho2, defocus_waves, spherical_waves=0.0):
    """Normalized PSF for defocus W020 and primary spherical W040 [waves]."""
    phase = np.float32(2.0 * np.pi * defocus_waves) * rho2
    if spherical_waves:
        phase = phase + np.float32(2.0 * np.pi * spherical_waves) * rho2 * rho2
    pupil = mask * np.exp(1j * phase)          # complex64
    field = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(pupil)))
    psf = np.abs(field).astype(np.float32)
    psf **= 2
    return psf / psf.max()


def crop_center(img, half):
    c = img.shape[0] // 2
    return img[c - half:c + half, c - half:c + half]


def ee_diameter_um(psf, frac=0.8):
    """Diameter [um] of the circle enclosing `frac` of the PSF energy."""
    n = psf.shape[0]
    ax = np.arange(n, dtype=np.float32) - n // 2
    r = np.hypot(ax[None, :], ax[:, None]).astype(np.int32)
    prof = np.bincount(r.ravel(), weights=psf.ravel().astype(np.float64))
    cum = np.cumsum(prof)
    r_frac = np.searchsorted(cum, frac * cum[-1])
    return 2.0 * float(r_frac) * PIX_UM


# ----------------------------------------------------------------------
# Physical scales
# ----------------------------------------------------------------------
PIX_UM = WAVELENGTH * F_NUMBER / PAD_FACTOR * 1e6   # focal-plane pixel [um]
PIX_ARCSEC = np.degrees(WAVELENGTH / (PAD_FACTOR * D_PRIMARY)) * 3600.0
DZ_PER_WAVE_MM = 8.0 * WAVELENGTH * F_NUMBER**2 * 1e3  # focus shift per wave of W020

def blur_half_px(w):
    """Geometric defocus blur radius in image pixels (~4*W*lambda/D)."""
    return 4.0 * abs(w) * PAD_FACTOR


def show(ax, img, half, title, mode):
    extent_um = half * PIX_UM
    if mode == "log":            # reveals faint rings + spikes
        data = np.log10(img + 1e-8)
        vmin, vmax = -6.0, 0.0
    else:                        # gamma stretch: vane shadows stay visible
        data = img ** 0.5
        vmin, vmax = 0.0, 1.0
    ax.imshow(data, cmap="inferno", origin="lower", vmin=vmin, vmax=vmax,
              extent=[-extent_um, extent_um, -extent_um, extent_um])
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("focal plane [µm]", fontsize=8)
    ax.tick_params(labelsize=7)


def main():
    mask, rho2 = make_pupil()

    print(f"Aperture D = {D_PRIMARY*1e3:.0f} mm, obstruction = {OBSTRUCTION:.2f}, "
          f"{len(VANE_ANGLES)} vanes x {VANE_WIDTH*1e3:.1f} mm")
    print(f"lambda = {WAVELENGTH*1e9:.0f} nm, f/{F_NUMBER:g}")
    print(f"Image pixel: {PIX_UM:.2f} um = {PIX_ARCSEC:.3f} arcsec")
    print(f"Airy radius: {1.22*WAVELENGTH*F_NUMBER*1e6:.2f} um")
    print(f"Defocus scale: 1 wave (W020) ~ {DZ_PER_WAVE_MM:.2f} mm focus shift")

    # ---------------- figure 1: pupil + defocus gallery ----------------
    ncols = len(DEFOCUS_SET) + 1
    fig, axes = plt.subplots(1, ncols, figsize=(3.1 * ncols, 3.6))

    half_pupil = N_PUPIL // 2 + 40
    R_mm = D_PRIMARY / 2 * 1e3
    axes[0].imshow(crop_center(mask, half_pupil), cmap="gray", origin="lower",
                   extent=[-R_mm * 1.08, R_mm * 1.08, -R_mm * 1.08, R_mm * 1.08])
    axes[0].set_title("Pupil: annulus + 3-vane spider", fontsize=10)
    axes[0].set_xlabel("pupil [mm]", fontsize=8)
    axes[0].tick_params(labelsize=7)

    for ax, w in zip(axes[1:], DEFOCUS_SET):
        psf = compute_psf(mask, rho2, w)
        half = int(max(80, 1.5 * blur_half_px(w)))
        title = ("In focus (log stretch)" if w == 0
                 else f"Defocus {w:g} waves  (Δz≈{w*DZ_PER_WAVE_MM:.1f} mm)")
        show(ax, crop_center(psf, half), half, title, "log" if w == 0 else "gamma")

    fig.suptitle("Cassegrain point source: effect of 3-leg secondary support vs defocus",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "cassegrain_pupil_and_psfs.png", dpi=150)
    print("saved cassegrain_pupil_and_psfs.png")

    # ---------------- figure 2: intra vs extra focal --------------------
    fig2, axes2 = plt.subplots(1, 2, figsize=(8.5, 4.9))
    for ax, w in zip(axes2, [-STAR_TEST_W, STAR_TEST_W]):
        psf = compute_psf(mask, rho2, w)
        half = int(1.5 * blur_half_px(w))
        side = "Intra-focal" if w < 0 else "Extra-focal"
        show(ax, crop_center(psf, half), half,
             f"{side}  (W020 = {w:+g} waves)", "gamma")
    fig2.suptitle("Star test: vane-shadow pattern flips (rotates 60°) through focus",
                  fontsize=11)
    fig2.tight_layout()
    fig2.savefig(FIG_DIR / "cassegrain_star_test.png", dpi=150)
    print("saved cassegrain_star_test.png")

    # ------- figure 3: through-focus series, negative and positive ------
    n_cols = len(GALLERY_SET) + 1                # + geometric-limit column
    fig3, axes3 = plt.subplots(2, n_cols, figsize=(2.5 * n_cols, 6.2))
    ee_data = [(0.0, ee_diameter_um(compute_psf(mask, rho2, 0.0)))]
    for col, aw in enumerate(GALLERY_SET):
        m, r2 = make_pupil(n_pupil_for(aw))
        half = int(max(100, 1.5 * blur_half_px(aw)))
        for row, sign in ((0, -1), (1, +1)):
            w = sign * aw
            psf = compute_psf(m, r2, w)
            if sign > 0:
                ee_data.append((aw, ee_diameter_um(psf)))
            ax = axes3[row, col]
            show(ax, crop_center(psf, half), half,
                 f"W020 = {w:+g} waves", "gamma")
            ax.set_xlabel(f"Δz ≈ {w * DZ_PER_WAVE_MM:+.1f} mm", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
            del psf

    # Geometric limit (|W| -> infinity): scaled pupil shadow; the two focus
    # sides are exact 180-degree rotations of one another.
    half = N_PUPIL // 2 + 30
    shadow = crop_center(mask, half)
    for row, img, tag in ((0, shadow[::-1, ::-1], "−∞"), (1, shadow, "+∞")):
        ax = axes3[row, -1]
        ax.imshow(img, cmap="inferno", origin="lower", vmin=0.0, vmax=1.0)
        ax.set_title(f"W020 → {tag}", fontsize=10)
        ax.set_xlabel("scaled pupil shadow", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])

    axes3[0, 0].set_ylabel("intra-focal  (W020 < 0)", fontsize=10)
    axes3[1, 0].set_ylabel("extra-focal  (W020 > 0)", fontsize=10)
    fig3.suptitle("Through-focus series: intra vs extra focal images are 180° "
                  "rotations of each other — the 3-vane shadow pattern appears "
                  "twisted by 60° across focus", fontsize=12)
    fig3.tight_layout(h_pad=2.0)
    fig3.savefig(FIG_DIR / "cassegrain_defocus_series.png", dpi=150)
    print("saved cassegrain_defocus_series.png")

    # ---------------- figure 4: spot size vs defocus --------------------
    ws = np.array([w for w, _ in ee_data if w > 0])
    ds = np.array([d for w, d in ee_data if w > 0])
    d0 = ee_data[0][1]                        # in-focus EE80 diameter
    w_line = np.logspace(np.log10(0.3), np.log10(300), 200)
    d_geo = 8.0 * w_line * WAVELENGTH * F_NUMBER * 1e6

    fig4, ax4 = plt.subplots(figsize=(6.6, 4.8))
    ax4.loglog(w_line, np.hypot(d_geo, d0), "-", color="tab:blue",
               label="quadrature model")
    ax4.loglog(w_line, d_geo, "--", color="gray",
               label="geometric: ⌀ = 8·W020·λ·F#")
    ax4.axhline(d0, color="gray", ls=":", lw=1,
                label=f"diffraction floor ({d0:.0f} µm)")
    ax4.loglog(ws, ds, "o", color="tab:red", ms=6, label="simulated EE80")
    ax4.set_xlabel("defocus |W020| [waves]  (same for either sign)")
    ax4.set_ylabel("80% encircled-energy diameter [µm]")
    secax = ax4.secondary_xaxis(
        "top", functions=(lambda w: w * DZ_PER_WAVE_MM,
                          lambda z: z / DZ_PER_WAVE_MM))
    secax.set_xlabel("focus shift |Δz| [mm]")
    ax4.grid(True, which="both", alpha=0.3)
    ax4.legend(fontsize=8, loc="upper left")
    fig4.tight_layout()
    fig4.savefig(FIG_DIR / "cassegrain_spot_size.png", dpi=150)
    print("saved cassegrain_spot_size.png")

    # ------- figure 5: star test with spherical aberration --------------
    w40 = SPHERICAL_W040
    fig5, axes5 = plt.subplots(2, len(SA_TEST_SET),
                               figsize=(2.9 * len(SA_TEST_SET), 6.4))
    for col, aw in enumerate(SA_TEST_SET):
        m, r2 = make_pupil(n_pupil_for(aw + 2.0 * w40))
        half = int(max(100, 1.5 * blur_half_px(aw)))
        for row, sign in ((0, -1), (1, +1)):
            # defocus measured from best focus (W020 = -W040 balances rho^4)
            w020 = sign * aw - w40
            psf = compute_psf(m, r2, w020, w40)
            ax = axes5[row, col]
            show(ax, crop_center(psf, half), half,
                 f"{aw * sign:+g} waves from best focus", "gamma")
            ax.set_xticks([]); ax.set_yticks([])
            del psf
    axes5[0, 0].set_ylabel("intra-focal", fontsize=10)
    axes5[1, 0].set_ylabel("extra-focal", fontsize=10)
    fig5.suptitle(f"Star test with W040 = {w40:g} waves spherical aberration:\n"
                  "equal distances either side of best focus no longer match",
                  fontsize=11)
    fig5.tight_layout(h_pad=2.0)
    fig5.savefig(FIG_DIR / "cassegrain_star_test_spherical.png", dpi=150)
    print("saved cassegrain_star_test_spherical.png")

    print("\n  |W020| [waves]   |Δz| [mm]   EE80 diameter [µm]")
    for w, d in ee_data:
        print(f"  {w:14g}   {w * DZ_PER_WAVE_MM:9.1f}   {d:12.0f}")


if __name__ == "__main__":
    main()
