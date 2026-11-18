"""
Demonstration: strong spherical aberration puts a LARGE bright spot inside
the hole of the defocused donut.

Geometric picture: with undercorrection W040 > 0, pupil zone rho focuses at
defocus W020 = -2*W040*rho^2. A detector plane inside this caustic range
( -2*W040 <= W020 <= 0, relative to paraxial focus ) receives the zone
rho0 = sqrt(-W020/(2*W040)) ON AXIS -> a bright central spot coexisting
with the donut made by the other zones. Outside the caustic (and on the
whole opposite side of focus) there is no spot: the hole is dark.

Output: cassegrain_sa_central_spot.png + spot metrics on stdout.
"""

import numpy as np
import matplotlib.pyplot as plt

import cassegrain_defocus_sim as optics
from cassegrain_parameter_study import azprofile
from cassegrain_extended_source_sim import LOD_ARCSEC

W040_SCAN = 12.0                                  # strong undercorrection
SCAN_W020 = [-36.0, -24.0, -12.0, -6.0, 0.0, 12.0]  # planes vs paraxial focus
SCAN_TAGS = ["beyond marginal focus", "marginal focus", "mid-caustic",
             "inside caustic", "paraxial focus", "other side of focus"]
MAG_SET   = [1.0, 2.0, 4.0, 8.0, 16.0]            # W040 series at mid-caustic


def sa_psf(w020, w040):
    n_pup = optics.n_pupil_for(abs(w020) + 2.0 * abs(w040))
    m, r2 = optics.make_pupil(n_pup)
    return optics.compute_psf(m, r2, w020, w040)


def edge_radius_px(w020, w040):
    """Max geometric image radius over the annulus zones [px]."""
    rho = np.linspace(optics.OBSTRUCTION, 1.0, 200)
    s = 2.0 * w020 * rho + 4.0 * w040 * rho**3
    return 2.0 * optics.PAD_FACTOR * np.abs(s).max()


def spot_metrics(psf):
    """(I(0)/Imax, FWHM, diameter at 10% of centre) — sizes in lambda/D."""
    rs, prof = azprofile(psf, 60.0 * optics.PAD_FACTOR, dr=0.25)
    i0 = prof[0] / psf.max()
    if i0 < 0.02:
        return i0, np.nan, np.nan
    half_idx = np.nonzero(prof < 0.5 * prof[0])[0]
    ten_idx = np.nonzero(prof < 0.1 * prof[0])[0]
    fwhm = 2.0 * rs[half_idx[0]] / optics.PAD_FACTOR if half_idx.size else np.nan
    d10 = 2.0 * rs[ten_idx[0]] / optics.PAD_FACTOR if ten_idx.size else np.nan
    return i0, fwhm, d10


def central_fill(psf, w020, w040):
    """Mean intensity of the inner hole region (axial spike excluded)
    relative to the outer-rim zone — 'how filled is the hole'."""
    r_edge = edge_radius_px(w020, w040)
    rs, prof = azprofile(psf, 1.05 * r_edge, dr=0.5)
    rim_band = (rs > 0.82 * r_edge) & (rs < r_edge)
    inner_band = (rs > 3.0 * optics.PAD_FACTOR) & (rs < 0.45 * r_edge)
    if not (rim_band.any() and inner_band.any()):
        return np.nan
    return float(prof[inner_band].mean() / prof[rim_band].mean())


def show(ax, psf, half, title, tag):
    crop = optics.crop_center(psf, max(half, 80))
    ax.imshow((crop / crop.max()) ** 0.5, cmap="inferno", origin="lower",
              vmin=0.0, vmax=1.0)
    ax.set_title(title, fontsize=8.5)
    ax.set_xlabel(tag, fontsize=7.5)
    ax.set_xticks([]); ax.set_yticks([])


def main():
    fig, axes = plt.subplots(2, len(SCAN_W020), figsize=(2.45 * len(SCAN_W020), 6.1))

    # Row A: fixed strong SA, detector plane scanned through the caustic.
    print(f"W040 = {W040_SCAN:g} waves; caustic spans W020 in "
          f"[{-2 * W040_SCAN:g}, 0] (paraxial reference)")
    ref = central_fill(sa_psf(16.0, 0.0), 16.0, 0.0)
    print(f"reference hole fill, pure defocus W020=16: {ref:.2f}")
    print("  W020    on-axis zone rho0   I(0)/Imax   ⌀10% [λ/D]   hole fill")
    for ax, w020, tag in zip(axes[0], SCAN_W020, SCAN_TAGS):
        psf = sa_psf(w020, W040_SCAN)
        half = int(1.35 * edge_radius_px(w020, W040_SCAN))
        show(ax, psf, half, f"W020 = {w020:+g} λ", tag)
        i0, fwhm, d10 = spot_metrics(psf)
        fill = central_fill(psf, w020, W040_SCAN)
        r0sq = -w020 / (2.0 * W040_SCAN)
        zone = (f"{np.sqrt(r0sq):.2f}"
                if optics.OBSTRUCTION**2 <= r0sq <= 1.0 else "  — ")
        d_s = "—" if np.isnan(d10) else f"{d10:6.1f}"
        print(f"  {w020:+5g}        {zone:>5s}          {i0:6.2f}   "
              f"{d_s:>9s}   {fill:8.2f}")
    axes[0, 0].set_ylabel(f"W040 = {W040_SCAN:g} waves,\nscan through focus",
                          fontsize=9)

    # Row B: SA magnitude series at mid-caustic (balanced focus) W020 = -W040.
    for ax, w040 in zip(axes[1], MAG_SET):
        psf = sa_psf(-w040, w040)
        half = int(1.35 * edge_radius_px(-w040, w040))
        fill = central_fill(psf, -w040, w040)
        tag = "hole fill —" if np.isnan(fill) else f"hole fill {fill:.2f}"
        show(ax, psf, half, f"W040 = {w040:g} λ", tag)
    axes[1, 0].set_ylabel("mid-caustic,\nincreasing W040", fontsize=9)

    # last cell: zoom on the paraxial-plane central blob (the fattest one)
    psf = sa_psf(0.0, W040_SCAN)
    show(axes[1, -1], psf, int(0.45 * edge_radius_px(0.0, W040_SCAN)),
         f"zoom: W040 = {W040_SCAN:g} λ, paraxial plane",
         "broad blob from the inner pupil zones")

    fig.suptitle("Spherical aberration filling the donut hole: the on-axis "
                 "caustic exists only inside [marginal, paraxial] focus",
                 fontsize=11.5)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_sa_central_spot.png", dpi=150)
    print("saved cassegrain_sa_central_spot.png")


if __name__ == "__main__":
    main()
