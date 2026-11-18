"""
Parameter study: what sets the SIZE of the defocused donut, and what do the
central obstruction and the spider vane width actually change?

Analytic expectations (geometric optics, verified here by simulation)
----------------------------------------------------------------------
Donut OUTER diameter:   d_out = 8 * W020 * lambda * F#  =  Delta_z / F#
    -> grows linearly with defocus; scales with aperture/focal length via
       F# = f/D;  for a FIXED physical defocus Delta_z it does not depend
       on wavelength at all (lambda only converts Delta_z into waves and
       sets the fringe fine-structure).
Hole diameter:          d_hole = OBSTRUCTION * d_out
Vane shadow width:      geometric (VANE_WIDTH / D) * d_out, BUT the lane
    cannot be narrower than the Fresnel diffraction scale sqrt(lambda*dz);
    for practical defocus the geometric width is below that floor, so the
    lane width saturates near sqrt(lambda*dz) and the vane width instead
    controls the shadow DEPTH (deficit ~ vane area spread over the floor).
    Obstruction and vane width NEVER change the donut size, they only
    re-shape the intensity inside it.
In focus, the obstruction pushes energy from the Airy core into the rings,
and wider vanes make brighter (~ width^2) but angularly shorter spikes.

Outputs
-------
cassegrain_donut_size_study.png   : measured vs predicted donut/hole size
cassegrain_obstruction_study.png  : obstruction ratio gallery (focus + donut)
cassegrain_vane_width_study.png   : vane width gallery (focus + donut)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d, map_coordinates

import cassegrain_defocus_sim as optics

# ----------------------------------------------------------------------
# User parameters
# ----------------------------------------------------------------------
OBSTRUCTION_SET = [0.15, 0.25, 0.35, 0.50]     # secondary/primary ratio
VANE_WIDTH_SET  = [0.001, 0.002, 0.004, 0.008]  # [m]
SIZE_W_SET      = [4.0, 8.0, 16.0, 32.0]        # defocus for size scaling
GALLERY_W       = 16.0                          # defocus for the galleries


def defocused_psf(w, obstruction=None, vane_width=None):
    m, r2 = optics.make_pupil(optics.n_pupil_for(max(w, 1.0)),
                              obstruction, vane_width)
    return optics.compute_psf(m, r2, w)


def azprofile(img, r_max, dr=0.5):
    """Azimuthally averaged profile at sub-pixel radial steps."""
    c = img.shape[0] // 2
    rs = np.arange(0.0, r_max, dr)
    th = np.linspace(0.0, 2.0 * np.pi, 720, endpoint=False)
    vals = map_coordinates(img, [c + rs[:, None] * np.sin(th)[None, :],
                                 c + rs[:, None] * np.cos(th)[None, :]],
                           order=1)
    return rs, vals.mean(axis=1)


def donut_radii_px(psf, w):
    """(hole_radius, outer_radius) in px from 50%-of-plateau crossings."""
    r0 = optics.blur_half_px(w)
    rs, prof = azprofile(psf, 1.35 * r0)
    prof = gaussian_filter1d(prof, sigma=4)          # ~2 px smoothing
    band = (rs > 0.55 * r0) & (rs < 0.85 * r0)
    plateau = np.median(prof[band])
    above = prof > 0.5 * plateau
    outer = rs[np.nonzero(above)[0][-1]]
    interior = np.nonzero(~above & (rs < 0.55 * r0))[0]
    hole = rs[interior[-1]] if interior.size else 0.0
    return hole, outer


def vane_shadow_fwhm_px(psf, w):
    """(FWHM in px, relative depth) of the 210-deg vane shadow lane."""
    c = psf.shape[0] // 2
    r = 0.62 * optics.blur_half_px(w)
    phi = np.deg2rad(np.linspace(150.0, 270.0, 1921))
    vals = map_coordinates(psf, [c + r * np.sin(phi), c + r * np.cos(phi)],
                           order=1)
    vals = gaussian_filter1d(vals, 3)
    base = np.percentile(vals, 80)
    win = np.abs(phi - np.deg2rad(210.0)) < np.deg2rad(15.0)
    i0 = np.nonzero(win)[0][0] + int(np.argmin(vals[win]))
    half = 0.5 * (base + vals[i0])
    lo = i0
    while lo > 0 and vals[lo - 1] < half:
        lo -= 1
    hi = i0
    while hi < vals.size - 1 and vals[hi + 1] < half:
        hi += 1
    return (phi[hi] - phi[lo]) * r, 1.0 - vals[i0] / base


def spike_excess(psf_infocus):
    """In-focus spike intensity 6-12 lambda/D out along phi=0, minus the
    inter-spike background at the same radii (Airy rings cancel)."""
    c = psf_infocus.shape[0] // 2
    rs = np.arange(6.0, 12.0, 0.25) * optics.PAD_FACTOR
    on = map_coordinates(psf_infocus, [np.full_like(rs, float(c)), c + rs],
                         order=1)
    a = np.deg2rad(30.0)
    off = map_coordinates(psf_infocus, [c + rs * np.sin(a), c + rs * np.cos(a)],
                          order=1)
    return float((on - off).mean())


def panel(ax, img, w, title, mode):
    half = int(max(100, 1.5 * optics.blur_half_px(w)))
    optics.show(ax, optics.crop_center(img, half), half, title, mode)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("")


def main():
    lam_fn_um = optics.WAVELENGTH * optics.F_NUMBER * 1e6
    print("donut outer diameter: 8*W020*lambda*F# = Delta_z/F# "
          "(independent of obstruction & vanes)")
    print("hole = OBSTRUCTION x outer;  vane shadow width: max(geometric, "
          "Fresnel floor sqrt(lambda*dz))\n")

    # ---------------- donut size: measured vs predicted -----------------
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(OBSTRUCTION_SET)))
    print("  eps    W020   outer meas/pred [um]     hole meas/pred [um]")
    for eps, col in zip(OBSTRUCTION_SET, colors):
        outs, holes = [], []
        for w in SIZE_W_SET:
            psf = defocused_psf(w, obstruction=eps)
            hole, outer = donut_radii_px(psf, w)
            d_out = 2.0 * outer * optics.PIX_UM
            d_hole = 2.0 * hole * optics.PIX_UM
            pred = 8.0 * w * lam_fn_um
            outs.append(d_out); holes.append(d_hole)
            print(f"  {eps:.2f} {w:6g}   {d_out:8.0f} /{pred:6.0f}      "
                  f"{d_hole:8.0f} /{eps * pred:6.0f}")
        ax.loglog(SIZE_W_SET, outs, "o-", color=col, ms=5)
        ax.loglog(SIZE_W_SET, holes, "s--", color=col, ms=5)
    w_line = np.array([3.0, 40.0])
    ax.loglog(w_line, 8.0 * w_line * lam_fn_um, "k:", lw=1.5)
    handles = ([Line2D([], [], color=c, label=f"obstruction = {e:.2f}")
                for e, c in zip(OBSTRUCTION_SET, colors)] +
               [Line2D([], [], color="k", marker="o", ls="-", label="outer ⌀"),
                Line2D([], [], color="k", marker="s", ls="--", label="hole ⌀"),
                Line2D([], [], color="k", ls=":", label="8·W020·λ·F#")])
    ax.legend(handles=handles, fontsize=8)
    ax.set_xlabel("defocus W020 [waves]")
    ax.set_ylabel("diameter [µm]")
    ax.set_title("Donut size: outer ⌀ set only by defocus; hole ⌀ = ε × outer",
                 fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_donut_size_study.png", dpi=150)
    print("saved cassegrain_donut_size_study.png")

    # ---------------- obstruction gallery --------------------------------
    fig2, axes2 = plt.subplots(2, len(OBSTRUCTION_SET),
                               figsize=(2.9 * len(OBSTRUCTION_SET), 6.4))
    for col_i, eps in enumerate(OBSTRUCTION_SET):
        panel(axes2[0, col_i], defocused_psf(0.0, obstruction=eps), 0.0,
              f"ε = {eps:.2f}", "log")
        panel(axes2[1, col_i], defocused_psf(GALLERY_W, obstruction=eps),
              GALLERY_W, "", "gamma")
    axes2[0, 0].set_ylabel("in focus (log)", fontsize=10)
    axes2[1, 0].set_ylabel(f"W020 = {GALLERY_W:g} waves", fontsize=10)
    fig2.suptitle("Obstruction ratio: hole ⌀ = ε × donut ⌀ (outer size fixed);\n"
                  "in focus: energy moves from Airy core into the rings",
                  fontsize=11)
    fig2.tight_layout()
    fig2.savefig(optics.FIG_DIR / "cassegrain_obstruction_study.png", dpi=150)
    print("saved cassegrain_obstruction_study.png")

    # ---------------- vane width gallery ---------------------------------
    fig3, axes3 = plt.subplots(2, len(VANE_WIDTH_SET),
                               figsize=(2.9 * len(VANE_WIDTH_SET), 6.4))
    dz_m = GALLERY_W * optics.DZ_PER_WAVE_MM * 1e-3
    fresnel_um = np.sqrt(optics.WAVELENGTH * dz_m) * 1e6
    print(f"\n  Fresnel width floor at W020={GALLERY_W:g}: "
          f"sqrt(lambda*dz) = {fresnel_um:.0f} um")
    print("  vane width | shadow FWHM meas [µm] (geom pred) | depth | spike excess")
    for col_i, wv in enumerate(VANE_WIDTH_SET):
        psf0 = defocused_psf(0.0, vane_width=wv)
        psfw = defocused_psf(GALLERY_W, vane_width=wv)
        panel(axes3[0, col_i], psf0, 0.0, f"vanes {wv * 1e3:g} mm", "log")
        panel(axes3[1, col_i], psfw, GALLERY_W, "", "gamma")
        fwhm, depth = vane_shadow_fwhm_px(psfw, GALLERY_W)
        pred_px = 8.0 * GALLERY_W * optics.PAD_FACTOR * wv / optics.D_PRIMARY
        print(f"  {wv * 1e3:6.0f} mm  |  {fwhm * optics.PIX_UM:8.0f} /"
              f"{pred_px * optics.PIX_UM:5.0f}        |  {depth:.2f} | "
              f"{spike_excess(psf0):.2e}")
    axes3[0, 0].set_ylabel("in focus (log)", fontsize=10)
    axes3[1, 0].set_ylabel(f"W020 = {GALLERY_W:g} waves", fontsize=10)
    fig3.suptitle("Vane width: donut size unchanged — shadow lane width is "
                  "diffraction-limited (≈√(λ·Δz)),\nso wider vanes deepen the "
                  "shadow; in focus: spike brightness grows ≈ w²", fontsize=11)
    fig3.tight_layout()
    fig3.savefig(optics.FIG_DIR / "cassegrain_vane_width_study.png", dpi=150)
    print("saved cassegrain_vane_width_study.png")


if __name__ == "__main__":
    main()
