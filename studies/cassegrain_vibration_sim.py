"""
Effect of line-of-sight VIBRATION on the defocused image of the Cassegrain
with a 3-vane spider, for point and extended sources.

Physics
-------
During an exposure much longer than the vibration period, the image is the
time-average of displaced copies of itself: a convolution with the kernel of
the line-of-sight motion. In the frequency domain this multiplies the image
spectrum by a vibration MTF:

  * random jitter (many uncorrelated cycles, RMS sigma):
        isotropic Gaussian kernel  ->  S(q) = exp(-2 pi^2 sigma^2 q^2)
        (monotonic suppression, like seeing)
  * sinusoidal vibration, amplitude a (half peak-to-peak) along one axis:
        arcsine kernel (dwell at the turning points)
        ->  S(q) = J0(2 pi a q_par),  q_par = frequency along the axis.
        DIRECTIONAL, oscillating sign: rings blur only along the vibration
        axis, vanish at the J0 zeros (first: 2 pi a q = 2.405) and REVIVE
        with reversed contrast in the negative lobes.

Both factors multiply the source-size factor 2J1(x)/x — attenuations of
rings by vibration and by source diameter compound multiplicatively.

Outputs
-------
cassegrain_vibration_gallery.png         : sinusoidal vs random jitter vs amplitude
cassegrain_vibration_ring_contrast.png   : ring contrast (par/perp) vs amplitude + theory
cassegrain_vibration_extended_source.png : vibration on top of a finite source
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0, j1
from scipy.ndimage import gaussian_filter1d, map_coordinates

import cassegrain_defocus_sim as optics
from cassegrain_extended_source_sim import disk_spectrum, LOD_ARCSEC

# ----------------------------------------------------------------------
# User parameters
# ----------------------------------------------------------------------
VIB_W          = 16.0                       # defocus [waves] for all figures
SIN_AMPS_LOD   = [0.0, 1.0, 2.0, 4.0, 8.0]  # sinusoidal amplitude a [lambda/D]
GAUSS_SIGS_LOD = [0.0, 0.5, 1.0, 2.0, 4.0]  # random-jitter RMS sigma [lambda/D]
SWEEP_LOD      = np.r_[0.0, np.geomspace(0.25, 10.0, 14)]  # for contrast curves
SRC_DIAM_LOD   = 2.0                        # extended-source case [lambda/D]
VIB_AXIS_DEG   = 0.0                        # vibration along +x

_FREQ_CACHE = {}


def freq_grids(n):
    """(fy, fx) in cycles/radian on the rfft2 grid."""
    if n not in _FREQ_CACHE:
        dtheta = optics.WAVELENGTH / (optics.PAD_FACTOR * optics.D_PRIMARY)
        _FREQ_CACHE[n] = (np.fft.fftfreq(n, dtheta)[:, None],
                          np.fft.rfftfreq(n, dtheta)[None, :])
    return _FREQ_CACHE[n]


def vib_spectrum(n, model, amp_lod, axis_deg=VIB_AXIS_DEG):
    """Vibration MTF on the rfft2 grid; amp in lambda/D units."""
    if amp_lod == 0.0:
        return 1.0
    fy, fx = freq_grids(n)
    a_rad = amp_lod * optics.WAVELENGTH / optics.D_PRIMARY
    if model == "sin":
        ph = np.deg2rad(axis_deg)
        q_par = fx * np.cos(ph) + fy * np.sin(ph)
        return j0(2.0 * np.pi * a_rad * q_par)
    return np.exp(-2.0 * np.pi**2 * a_rad**2 * (fx**2 + fy**2))


def image_from(psf_spectrum, n, factor):
    img = np.fft.fftshift(np.fft.irfft2(psf_spectrum * factor, s=(n, n)))
    np.clip(img, 0.0, None, out=img)
    return img / img.max()


def wedge_profile(img, r_max, phi0_deg, half_deg=12.0, dr=0.25):
    """Radial profile averaged over a narrow azimuthal wedge."""
    c = img.shape[0] // 2
    th = np.deg2rad(phi0_deg + np.linspace(-half_deg, half_deg, 41))
    rs = np.arange(0.0, r_max, dr)
    vals = map_coordinates(img, [c + rs[:, None] * np.sin(th)[None, :],
                                 c + rs[:, None] * np.cos(th)[None, :]],
                           order=3)
    return rs, vals.mean(axis=1)


def ring_contrast_dir(img, w, phi0_deg, dr=0.25):
    """Detrended ring contrast measured along one direction."""
    r0 = optics.blur_half_px(w)
    rs, prof = wedge_profile(img, 0.92 * r0, phi0_deg, dr=dr)
    band = (rs > 0.45 * r0) & (rs < 0.85 * r0)
    seg = prof[band]
    env = gaussian_filter1d(seg, sigma=8.0 / dr)
    fine = seg - env
    p_lo, p_hi = np.percentile(fine, [2.5, 97.5])
    return 0.5 * (p_hi - p_lo) / env.mean()


def ring_freq_px(img, w, phi0_deg, dr=0.25):
    """Dominant interior-ring radial frequency [cycles/px] in one wedge."""
    r0 = optics.blur_half_px(w)
    rs, prof = wedge_profile(img, 0.92 * r0, phi0_deg, dr=dr)
    band = (rs > 0.45 * r0) & (rs < 0.85 * r0)
    seg = prof[band]
    fine = seg - gaussian_filter1d(seg, sigma=8.0 / dr)
    spec = np.abs(np.fft.rfft(fine * np.hanning(fine.size)))
    freqs = np.fft.rfftfreq(fine.size, d=dr)
    return freqs[1 + int(np.argmax(spec[1:]))]


def panel(ax, img, w, amp_lod, title):
    half = int(1.5 * optics.blur_half_px(w) + 3.0 * amp_lod * optics.PAD_FACTOR)
    crop = optics.crop_center(img, half)
    ax.imshow(crop ** 0.5, cmap="inferno", origin="lower", vmin=0.0, vmax=1.0)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])


def main():
    m, r2 = optics.make_pupil(optics.n_pupil_for(VIB_W))
    psf = optics.compute_psf(m, r2, VIB_W)
    n = psf.shape[0]
    spectrum = np.fft.rfft2(np.fft.ifftshift(psf))
    clean = image_from(spectrum, n, 1.0)

    # interior ring frequency (parallel wedge, phi=0 lies between two vanes)
    q_px = ring_freq_px(clean, VIB_W, 0.0)
    period_lod = 1.0 / (q_px * optics.PAD_FACTOR)
    a_null_lod = 2.405 / (2.0 * np.pi * optics.PAD_FACTOR * q_px)
    print(f"W020 = {VIB_W:g} waves; interior ring period = {period_lod:.1f} "
          f"lambda/D = {period_lod * LOD_ARCSEC:.2f} arcsec")
    print(f"predicted first ring null (sinusoidal, along axis): "
          f"a = {a_null_lod:.1f} lambda/D = {a_null_lod * LOD_ARCSEC:.2f} arcsec")

    # ---------------- gallery: sinusoidal vs random jitter ---------------
    n_cols = len(SIN_AMPS_LOD)
    fig, axes = plt.subplots(2, n_cols, figsize=(2.7 * n_cols, 6.0))
    for col, a in enumerate(SIN_AMPS_LOD):
        img = image_from(spectrum, n, vib_spectrum(n, "sin", a))
        t = "no vibration" if a == 0 else f"a = {a:g} λ/D ({a * LOD_ARCSEC:.1f}″)"
        panel(axes[0, col], img, VIB_W, a, t)
    for col, s in enumerate(GAUSS_SIGS_LOD):
        img = image_from(spectrum, n, vib_spectrum(n, "gauss", s))
        t = "no jitter" if s == 0 else f"σ = {s:g} λ/D ({s * LOD_ARCSEC:.1f}″)"
        panel(axes[1, col], img, VIB_W, s, t)
    axes[0, 0].set_ylabel("sinusoidal, along ↔", fontsize=10)
    axes[1, 0].set_ylabel("random jitter (isotropic)", fontsize=10)
    fig.suptitle(f"Vibration blur at W020 = {VIB_W:g} waves (point source): "
                 "directional smear vs isotropic smear", fontsize=12)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_vibration_gallery.png", dpi=150)
    print("saved cassegrain_vibration_gallery.png")

    # ---------------- ring contrast vs amplitude + theory ----------------
    fig2, (axs, axg) = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    c0 = ring_contrast_dir(clean, VIB_W, 0.0)
    results = {}
    for model, ax, sym in (("sin", axs, "a"), ("gauss", axg, "σ")):
        meas_par, meas_perp = [], []
        for amp in SWEEP_LOD:
            img = image_from(spectrum, n, vib_spectrum(n, model, amp))
            meas_par.append(ring_contrast_dir(img, VIB_W, 0.0) / c0)
            meas_perp.append(ring_contrast_dir(img, VIB_W, 270.0) / c0)
        results[model] = (np.array(meas_par), np.array(meas_perp))
        a_line = np.linspace(0.0, SWEEP_LOD.max(), 400)
        x = 2.0 * np.pi * (a_line * optics.PAD_FACTOR) * q_px
        if model == "sin":
            theory = np.abs(j0(x))
            ax.axvline(a_null_lod, color="gray", lw=0.8, ls=":")
            ax.text(a_null_lod, 1.02, "1st null", fontsize=7.5,
                    ha="center", color="gray")
        else:
            theory = np.exp(-x**2 / 4.0)          # exp(-2 pi^2 s^2 q^2)
        ax.plot(a_line, theory, "-", color="tab:blue",
                label="theory at ring frequency")
        ax.plot(SWEEP_LOD, results[model][0], "o", color="tab:red", ms=5,
                label="measured ∥ vibration axis")
        ax.plot(SWEEP_LOD, results[model][1], "s", color="tab:green", ms=5,
                label="measured ⊥ vibration axis")
        ax.set_xlabel(f"{sym} [λ/D]  (top: arcsec)")
        secax = ax.secondary_xaxis(
            "top", functions=(lambda v: v * LOD_ARCSEC,
                              lambda v: v / LOD_ARCSEC))
        secax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
        ax.set_ylim(-0.05, 1.1)
    axs.set_title("sinusoidal: |J0(2π·a·q)| — nulls & revivals", fontsize=10)
    axg.set_title("random jitter: exp(−2π²σ²q²) — monotonic", fontsize=10)
    axs.set_ylabel("interior-ring contrast (normalized)")
    axs.legend(fontsize=8)
    fig2.suptitle("Effect of vibration on the ring pattern "
                  f"(W020 = {VIB_W:g} waves)", fontsize=12)
    fig2.tight_layout()
    fig2.savefig(optics.FIG_DIR / "cassegrain_vibration_ring_contrast.png", dpi=150)
    print("saved cassegrain_vibration_ring_contrast.png")

    print("\n  sinusoidal vibration (along axis):")
    print("  a [λ/D]   meas ∥    |J0| theory   meas ⊥")
    for amp, mp, mq in zip(SWEEP_LOD, *results["sin"]):
        th = abs(j0(2.0 * np.pi * amp * optics.PAD_FACTOR * q_px))
        print(f"  {amp:7.2f}   {mp:6.2f}   {th:10.2f}   {mq:6.2f}")

    # ---------------- vibration on an extended source --------------------
    src = disk_spectrum(n, SRC_DIAM_LOD)
    cases = [("no vibration", 1.0),
             (f"random σ = 1 λ/D", vib_spectrum(n, "gauss", 1.0)),
             (f"sinusoidal a = 2 λ/D", vib_spectrum(n, "sin", 2.0))]
    fig3, axes3 = plt.subplots(2, 3, figsize=(8.4, 6.0))
    for col, (label, fac) in enumerate(cases):
        panel(axes3[0, col], image_from(spectrum, n, fac), VIB_W, 2.0, label)
        panel(axes3[1, col], image_from(spectrum, n, fac * src), VIB_W, 2.0, "")
    axes3[0, 0].set_ylabel("point source", fontsize=10)
    axes3[1, 0].set_ylabel(f"source ⌀ {SRC_DIAM_LOD:g} λ/D "
                           f"({SRC_DIAM_LOD * LOD_ARCSEC:.1f}″)", fontsize=10)
    fig3.suptitle("Vibration and source size compound multiplicatively "
                  f"(W020 = {VIB_W:g} waves)", fontsize=11)
    fig3.tight_layout()
    fig3.savefig(optics.FIG_DIR / "cassegrain_vibration_extended_source.png", dpi=150)
    print("saved cassegrain_vibration_extended_source.png")

    # multiplicativity check at the ring frequency
    img_s = image_from(spectrum, n, src)
    img_sv = image_from(spectrum, n, src * vib_spectrum(n, "sin", 2.0))
    c_s = ring_contrast_dir(img_s, VIB_W, 0.0) / c0
    c_sv = ring_contrast_dir(img_sv, VIB_W, 0.0) / c0
    c_v = np.interp(2.0, SWEEP_LOD, results["sin"][0])
    print(f"\n  multiplicativity: source(2λ/D) {c_s:.2f} x vib(a=2) {c_v:.2f} "
          f"= {c_s * c_v:.2f}  vs measured combined {c_sv:.2f}")


if __name__ == "__main__":
    main()
