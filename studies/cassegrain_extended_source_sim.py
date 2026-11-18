"""
Defocused image of an EXTENDED (non point-like) source in the Cassegrain
reflector with a 3-vane spider.

Physics
-------
For a spatially incoherent source (star, planet, thermal pinhole) the image
is the convolution of the point-source PSF with the geometric image of the
source:

    I(theta) = PSF(theta) (x) Source(theta)

Here the source is a uniform disk of angular diameter theta_d. The
convolution is done in the frequency domain, where the disk contributes the
classic visibility factor  2*J1(x)/x,  x = pi * theta_d * q  (exact for any
size, no rasterization error).

Question answered
-----------------
Does the source size affect the ring pattern of the defocused donut?  Yes,
but with a twist: every angular frequency q in the image is attenuated by
the uniform-disk visibility 2*J1(pi*theta_d*q)/(pi*theta_d*q), and the
donut carries several distinct spatial scales:
  * the finest fringes hugging the shadow edges reach the cutoff D/lambda
    -> first to vanish, already for theta_d ~ lambda/D;
  * the INTERIOR rings are edge-wave beats between pupil zones separated by
    only a fraction of D, so their angular period is several lambda/D
    -> they survive until theta_d reaches a few lambda/D;
  * vane shadows and the central hole are coarser still and scale UP with
    defocus -> they are the last features to wash out.
Ring periods are fixed in angle (independent of defocus), so the source
size that erases them is roughly defocus-independent. Atmospheric seeing
acts exactly like an extra source diameter.

Outputs
-------
cassegrain_extended_source_gallery.png  : source size x defocus image grid
cassegrain_extended_source_rings.png    : zoom on the ring region of each image
cassegrain_extended_source_profiles.png : radial profiles, ring contrast decay
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1
from scipy.ndimage import gaussian_filter1d, map_coordinates

import cassegrain_defocus_sim as optics   # pupil / PSF machinery + parameters

# ----------------------------------------------------------------------
# User parameters
# ----------------------------------------------------------------------
SOURCE_DIAMS_LOD = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]   # source diameter [lambda/D]
GRID_DEFOCUS     = [4.0, 8.0, 16.0, 32.0]      # W020 [waves] for the gallery
PROFILE_W        = 16.0                         # defocus for the profile figure
METRICS_W        = 32.0                         # defocus for ring metrics (many rings)

LOD_ARCSEC = np.degrees(optics.WAVELENGTH / optics.D_PRIMARY) * 3600.0

_SPEC_CACHE = {}


def disk_spectrum(n, diam_lod):
    """Uniform-disk visibility 2*J1(x)/x on the rfft2 frequency grid."""
    key = (n, diam_lod)
    if key in _SPEC_CACHE:
        return _SPEC_CACHE[key]
    dtheta = optics.WAVELENGTH / (optics.PAD_FACTOR * optics.D_PRIMARY)  # rad/px
    theta_d = diam_lod * optics.WAVELENGTH / optics.D_PRIMARY            # rad
    fy = np.fft.fftfreq(n, d=dtheta)
    fx = np.fft.rfftfreq(n, d=dtheta)
    q = np.hypot(fy[:, None], fx[None, :])          # cycles/radian
    x = np.pi * theta_d * q
    S = np.ones_like(x)
    nz = x > 1e-12
    S[nz] = 2.0 * j1(x[nz]) / x[nz]
    _SPEC_CACHE[key] = S
    return S


def extended_image(psf_spectrum, n, diam_lod):
    """Image of a uniform-disk source from the precomputed PSF spectrum."""
    img = np.fft.irfft2(psf_spectrum * disk_spectrum(n, diam_lod), s=(n, n))
    img = np.fft.fftshift(img)
    np.clip(img, 0.0, None, out=img)
    return img


def fine_radial_profile(img, r_max, dr=0.25, n_theta=720):
    """Azimuthally averaged profile sampled at sub-pixel radial steps.

    The ring period is ~ lambda/D = PAD_FACTOR pixels (2 px), so integer-
    pixel binning would alias it away; cubic-spline sampling keeps it.
    """
    c = img.shape[0] // 2
    rs = np.arange(0.0, r_max, dr)
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    Y = c + rs[:, None] * np.sin(th)[None, :]
    X = c + rs[:, None] * np.cos(th)[None, :]
    vals = map_coordinates(img, [Y, X], order=3)
    return rs, vals.mean(axis=1)


def ring_metrics(img_crop, w, dr=0.25):
    """(contrast, dominant radial frequency [cycles/px]) of the fine rings."""
    r0 = optics.blur_half_px(w)
    rs, prof = fine_radial_profile(img_crop, 0.92 * r0, dr)
    band = (rs > 0.45 * r0) & (rs < 0.85 * r0)
    seg = prof[band]
    # envelope smoothing ~1.5x the interior ring period (~6 lambda/D)
    env = gaussian_filter1d(seg, sigma=8.0 / dr)
    fine = seg - env
    p_lo, p_hi = np.percentile(fine, [2.5, 97.5])
    contrast = 0.5 * (p_hi - p_lo) / env.mean()
    spec = np.abs(np.fft.rfft(fine * np.hanning(fine.size)))
    freqs = np.fft.rfftfreq(fine.size, d=dr)        # cycles/px
    q_peak = freqs[1 + int(np.argmax(spec[1:]))]
    return contrast, q_peak


def vane_contrast(img_crop, w):
    """Azimuthal modulation on a circle crossing the three vane shadows."""
    c = img_crop.shape[0] // 2
    r = 0.62 * optics.blur_half_px(w)
    ang = np.linspace(0.0, 2.0 * np.pi, 721)
    vals = map_coordinates(img_crop, [c + r * np.sin(ang), c + r * np.cos(ang)],
                           order=1)
    vmin, vmax = vals.min(), vals.max()
    return (vmax - vmin) / (vmax + vmin)


def source_label(d):
    return "point source" if d == 0 else f"⌀ {d:g} λ/D ≈ {d * LOD_ARCSEC:.2f}″"


def hole_depth(img_crop, w):
    """Darkness of the obstruction shadow, sampled on a ring at half the
    hole radius to avoid the bright axial (Arago-like) diffraction spot."""
    c = img_crop.shape[0] // 2
    r = 0.5 * optics.OBSTRUCTION * optics.blur_half_px(w)
    ang = np.linspace(0.0, 2.0 * np.pi, 361)
    vals = map_coordinates(img_crop, [c + r * np.sin(ang), c + r * np.cos(ang)],
                           order=1)
    return 1.0 - float(vals.mean())


def main():
    print(f"lambda/D = {LOD_ARCSEC:.3f} arcsec  "
          f"(D = {optics.D_PRIMARY*1e3:.0f} mm, lambda = {optics.WAVELENGTH*1e9:.0f} nm)")
    print(f"source diameters: " +
          ", ".join(f"{d:g} λ/D = {d*LOD_ARCSEC:.2f}\"" for d in SOURCE_DIAMS_LOD))

    # ---------------- gallery: source size x defocus ---------------------
    n_rows, n_cols = len(SOURCE_DIAMS_LOD), len(GRID_DEFOCUS)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(2.8 * n_cols, 2.8 * n_rows + 0.8))
    images = {}
    for col, w in enumerate(GRID_DEFOCUS):
        m, r2 = optics.make_pupil(optics.n_pupil_for(w))
        psf = optics.compute_psf(m, r2, w)
        n = psf.shape[0]
        spectrum = np.fft.rfft2(np.fft.ifftshift(psf))
        half = int(max(100, 1.5 * optics.blur_half_px(w)))
        for row, d in enumerate(SOURCE_DIAMS_LOD):
            img = psf if d == 0 else extended_image(spectrum, n, d)
            crop = optics.crop_center(img, half).astype(np.float64)
            crop /= crop.max()
            images[(d, w)] = crop
            ax = axes[row, col]
            ax.imshow(crop ** 0.5, cmap="inferno", origin="lower",
                      vmin=0.0, vmax=1.0)
            if row == 0:
                ax.set_title(f"W020 = {w:g} waves\nΔz ≈ {w*optics.DZ_PER_WAVE_MM:.1f} mm",
                             fontsize=9)
            if col == 0:
                ax.set_ylabel(source_label(d), fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        del psf, spectrum

    fig.suptitle("Extended incoherent source: interior rings fade for sources "
                 "of a few λ/D; vane shadows and the donut hole survive far "
                 "larger sources", fontsize=12)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_extended_source_gallery.png", dpi=150)
    print("saved cassegrain_extended_source_gallery.png")

    # ---------------- zoom on the ring region ----------------------------
    # (rings are ~2 px wide, invisible in the shrunken overview panels)
    fig3, axes3 = plt.subplots(n_rows, n_cols,
                               figsize=(2.4 * n_cols, 2.4 * n_rows + 0.8))
    zoom_ang = np.deg2rad(150.0)              # mid-annulus, between two vanes
    for col, w in enumerate(GRID_DEFOCUS):
        for row, d in enumerate(SOURCE_DIAMS_LOD):
            img = images[(d, w)]
            c = img.shape[0] // 2
            r_mid = 0.66 * optics.blur_half_px(w)
            cx = int(round(c + r_mid * np.cos(zoom_ang)))
            cy = int(round(c + r_mid * np.sin(zoom_ang)))
            bh = max(10, int(0.30 * optics.blur_half_px(w)))
            box = img[cy - bh:cy + bh, cx - bh:cx + bh]
            ax = axes3[row, col]
            ax.imshow(box ** 0.5, cmap="inferno", origin="lower",
                      vmin=0.0, vmax=1.0, interpolation="bilinear")
            if row == 0:
                ax.set_title(f"W020 = {w:g} waves", fontsize=9)
            if col == 0:
                ax.set_ylabel(source_label(d), fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
    fig3.suptitle("Zoom on the donut ring structure (box in the bright annulus, "
                  "away from the vanes) — ring period is a few λ/D, set by "
                  "pupil-zone spacing, independent of defocus", fontsize=11)
    fig3.tight_layout()
    fig3.savefig(optics.FIG_DIR / "cassegrain_extended_source_rings.png", dpi=150)
    print("saved cassegrain_extended_source_rings.png")

    # ---------------- radial profiles at PROFILE_W -----------------------
    fig2, ax2 = plt.subplots(figsize=(8.2, 5.4))
    r0 = optics.blur_half_px(PROFILE_W)
    for i, d in enumerate(SOURCE_DIAMS_LOD):
        rs, prof = fine_radial_profile(images[(d, PROFILE_W)], 1.25 * r0)
        prof = prof / prof.max() + 0.35 * i          # waterfall offset
        ax2.plot(rs * optics.PIX_ARCSEC, prof, lw=1.0,
                 label=source_label(d))
    outer_as = 4.0 * PROFILE_W * LOD_ARCSEC          # geometric shadow radius
    ax2.axvline(outer_as, color="gray", ls="--", lw=1)
    ax2.axvline(optics.OBSTRUCTION * outer_as, color="gray", ls=":", lw=1)
    top = 1.1 + 0.35 * (len(SOURCE_DIAMS_LOD) - 1)
    ax2.text(outer_as, top, "outer edge", ha="center", fontsize=8, color="gray")
    ax2.text(optics.OBSTRUCTION * outer_as, top, "obstruction edge",
             ha="center", fontsize=8, color="gray")
    ax2.set_xlim(0.0, 1.25 * outer_as)
    ax2.set_ylim(0.0, top + 0.08)
    ax2.set_xlabel("radius [arcsec]")
    ax2.set_ylabel("azimuthal-mean intensity (offset per curve)")
    ax2.set_title(f"Radial profiles at W020 = {PROFILE_W:g} waves — interior "
                  f"ring period ≈ several λ/D  (λ/D = {LOD_ARCSEC:.2f}″)",
                  fontsize=11)
    ax2.legend(fontsize=8, loc="upper left")
    ax2.grid(alpha=0.25)
    fig2.tight_layout()
    fig2.savefig(optics.FIG_DIR / "cassegrain_extended_source_profiles.png", dpi=150)
    print("saved cassegrain_extended_source_profiles.png")

    # ---------------- contrast metrics -----------------------------------
    rc0, q_peak = ring_metrics(images[(0.0, METRICS_W)], METRICS_W)
    period_lod = 1.0 / (q_peak * optics.PAD_FACTOR)
    print(f"\n  ring analysis at W020 = {METRICS_W:g} waves:")
    print(f"  dominant interior ring period = {period_lod:.1f} λ/D "
          f"= {period_lod * LOD_ARCSEC:.2f} arcsec")
    print(f"\n  source ⌀             ring contrast   measured   theory      vane")
    print(f"                                        ratio      2J1(x)/x    contrast")
    for d in SOURCE_DIAMS_LOD:
        rc, _ = ring_metrics(images[(d, METRICS_W)], METRICS_W)
        x = np.pi * d * optics.PAD_FACTOR * q_peak   # at the measured ring freq
        v_th = 1.0 if d == 0 else abs(2.0 * j1(x) / x)
        vc = vane_contrast(images[(d, METRICS_W)], METRICS_W)
        print(f"  {source_label(d):<18s} {rc:10.3f} {rc / rc0:10.2f} "
              f"{v_th:10.2f} {vc:10.2f}")

    # ---------------- vane & hole survival matrices ----------------------
    vane_ratio = optics.VANE_WIDTH / optics.D_PRIMARY
    print("\n  vane-shadow contrast  (feature width = 8·W·w_vane/D)")
    print("  source ⌀           " +
          "".join(f"  W={w:<7g}" for w in GRID_DEFOCUS))
    print("  vane width [λ/D]   " +
          "".join(f"  {8.0 * w * vane_ratio:<9.1f}" for w in GRID_DEFOCUS))
    for d in SOURCE_DIAMS_LOD:
        cells = "".join(f"  {vane_contrast(images[(d, w)], w):<9.2f}"
                        for w in GRID_DEFOCUS)
        print(f"  {source_label(d):<18s}" + cells)

    print("\n  central-hole depth  (hole diameter = obstruction · 8·W)")
    print("  source ⌀           " +
          "".join(f"  W={w:<7g}" for w in GRID_DEFOCUS))
    print("  hole diam [λ/D]    " +
          "".join(f"  {optics.OBSTRUCTION * 8.0 * w:<9.1f}" for w in GRID_DEFOCUS))
    for d in SOURCE_DIAMS_LOD:
        cells = "".join(f"  {hole_depth(images[(d, w)], w):<9.2f}"
                        for w in GRID_DEFOCUS)
        print(f"  {source_label(d):<18s}" + cells)

    # Axial (Arago-like) bright spot at the exact shadow centre.
    print("\n  axial bright-spot intensity I(0) (relative to image peak)")
    print("  source ⌀           " +
          "".join(f"  W={w:<7g}" for w in GRID_DEFOCUS))
    for d in SOURCE_DIAMS_LOD:
        cells = "".join(
            f"  {float(images[(d, w)][images[(d, w)].shape[0] // 2, images[(d, w)].shape[0] // 2]):<9.2f}"
            for w in GRID_DEFOCUS)
        print(f"  {source_label(d):<18s}" + cells)


if __name__ == "__main__":
    main()
