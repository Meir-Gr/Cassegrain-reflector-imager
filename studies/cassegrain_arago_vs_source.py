"""
Can the Arago (axial) spot explain a LARGE circular spot inside the donut
hole, once the source is extended?

For an incoherent source the central axial feature is the point-source
Arago spot convolved with the source disk: a central IMAGE of the source.
Prediction: spot diameter = source angular diameter (independent of
defocus), with surface brightness that dilutes as the source grows.
This script measures both, at two defocus values.

Outputs
-------
cassegrain_arago_vs_source_gallery.png : donut + hole zooms vs source size
cassegrain_arago_vs_source_metrics.png : spot diameter & brightness vs source
"""

import numpy as np
import matplotlib.pyplot as plt

import cassegrain_defocus_sim as optics
from cassegrain_extended_source_sim import disk_spectrum, LOD_ARCSEC
from cassegrain_parameter_study import azprofile

W_SET   = [8.0, 16.0]                       # defocus [waves]
SRC_SET = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]   # source diameter [lambda/D]


def images_for(w):
    m, r2 = optics.make_pupil(optics.n_pupil_for(w))
    psf = optics.compute_psf(m, r2, w)
    n = psf.shape[0]
    spec = np.fft.rfft2(np.fft.ifftshift(psf))
    out = {}
    for d in SRC_SET:
        if d == 0.0:
            img = psf
        else:
            img = np.fft.fftshift(np.fft.irfft2(spec * disk_spectrum(n, d),
                                                s=(n, n)))
            np.clip(img, 0.0, None, out=img)
        out[d] = img / img.max()
    return out


def rim_level(img, w):
    r_edge = optics.blur_half_px(w)
    rs, prof = azprofile(img, 1.02 * r_edge, dr=0.5)
    return float(prof[(rs > 0.5 * r_edge) & (rs < 0.85 * r_edge)].mean())


def spot_metrics(img, w):
    """(centre / rim brightness, spot diameter [lambda/D]).

    Diameter = half-height crossing measured against the local hole floor
    (annulus 0.80-0.95 of the hole radius), immune to the axial rings.
    """
    hole_r = optics.OBSTRUCTION * optics.blur_half_px(w)
    rs, prof = azprofile(img, 0.98 * hole_r, dr=0.25)
    rim = rim_level(img, w)
    floor = float(np.percentile(
        prof[(rs > 0.80 * hole_r) & (rs < 0.95 * hole_r)], 20))
    center = float(prof[0])
    peak_rel = center / rim
    if center - floor < 0.02 * rim:
        return peak_rel, np.nan          # no distinct central spot
    half = floor + 0.5 * (center - floor)
    below = np.nonzero(prof < half)[0]
    diam = 2.0 * rs[below[0]] / optics.PAD_FACTOR if below.size else np.nan
    return peak_rel, diam


def main():
    galleries = {w: images_for(w) for w in W_SET}

    # ---------------- gallery ------------------------------------------
    n_cols = len(SRC_SET)
    fig, axes = plt.subplots(3, n_cols, figsize=(2.2 * n_cols, 7.0))
    w0 = W_SET[0]
    for col, d in enumerate(SRC_SET):
        img = galleries[w0][d]
        half = int(1.3 * optics.blur_half_px(w0))
        crop = optics.crop_center(img, half)
        axes[0, col].imshow((crop / crop.max()) ** 0.5, cmap="inferno",
                            origin="lower", vmin=0, vmax=1)
        axes[0, col].set_title("point" if d == 0 else
                               f"⌀ {d:g} λ/D ({d * LOD_ARCSEC:.1f}″)",
                               fontsize=8.5)
        for row, w in enumerate(W_SET, start=1):
            img = galleries[w][d]
            hole_r = optics.OBSTRUCTION * optics.blur_half_px(w)
            zoom = optics.crop_center(img, int(hole_r + 6)) / rim_level(img, w)
            axes[row, col].imshow(np.clip(zoom, 0, 1.4) ** 0.5, cmap="inferno",
                                  origin="lower", vmin=0, vmax=1.2,
                                  interpolation="bilinear")
        for row in range(3):
            axes[row, col].set_xticks([]); axes[row, col].set_yticks([])
    axes[0, 0].set_ylabel(f"donut, W020 = {w0:g} λ", fontsize=9)
    axes[1, 0].set_ylabel(f"hole zoom, W020 = {W_SET[0]:g} λ\n"
                          "(brightness relative to rim)", fontsize=8)
    axes[2, 0].set_ylabel(f"hole zoom, W020 = {W_SET[1]:g} λ\n"
                          "(brightness relative to rim)", fontsize=8)
    fig.suptitle("Arago spot vs source size: the central feature is an image "
                 "of the source — same size at both defocuses", fontsize=11)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_arago_vs_source_gallery.png", dpi=150)
    print("saved cassegrain_arago_vs_source_gallery.png")

    # ---------------- metrics ------------------------------------------
    fig2, (axd, axb) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    marks = {W_SET[0]: "o", W_SET[1]: "s"}
    print("\n  W020   source ⌀ [λ/D]   spot ⌀ [λ/D]   spot/hole   centre/rim")
    for w in W_SET:
        hole_d = 2.0 * optics.OBSTRUCTION * 4.0 * w      # in lambda/D
        ds, peaks = [], []
        for d in SRC_SET:
            peak, diam = spot_metrics(galleries[w][d], w)
            ds.append(diam); peaks.append(peak)
            frac = "—" if np.isnan(diam) else f"{diam / hole_d:8.2f}"
            dd = "—" if np.isnan(diam) else f"{diam:8.1f}"
            print(f"  {w:4g}   {d:12g}   {dd:>12s}   {frac:>8s}   {peak:8.2f}")
        axd.plot(SRC_SET, ds, marks[w] + "-", ms=5, label=f"W020 = {w:g} λ")
        axb.plot(SRC_SET, peaks, marks[w] + "-", ms=5, label=f"W020 = {w:g} λ")
    lim = max(SRC_SET)
    axd.plot([0, lim], [0, lim], "k--", lw=1, label="spot ⌀ = source ⌀")
    axd.set_xlabel("source diameter [λ/D]")
    axd.set_ylabel("central-spot diameter [λ/D]")
    axd.set_title("spot size = source size,\nindependent of defocus", fontsize=10)
    axd.legend(fontsize=8); axd.grid(alpha=0.3)
    axb.set_xlabel("source diameter [λ/D]")
    axb.set_ylabel("spot peak / rim brightness")
    axb.set_title("spot brightness relative to the rim", fontsize=10)
    axb.legend(fontsize=8); axb.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(optics.FIG_DIR / "cassegrain_arago_vs_source_metrics.png", dpi=150)
    print("saved cassegrain_arago_vs_source_metrics.png")


if __name__ == "__main__":
    main()
