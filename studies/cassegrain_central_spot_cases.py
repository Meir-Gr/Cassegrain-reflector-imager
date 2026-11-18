"""
Candidate explanations for a LARGE circular spot inside the hole of a
defocused Cassegrain donut, simulated side by side at W020 = 16 waves.

  (a) ghost reflection      : incoherent second pupil image with a different
                              defocus (sensor window / filter / corrector);
                              concentric for an on-axis star.
  (b) spherical aberration  : caustic concentrates light on axis on ONE side
                              of focus; blob size scales with the donut.
  (c) baffle leak           : unfocused direct light through the primary
                              hole; uniform disk, size independent of focus.
  (d) resolved source       : Arago spot convolved with the source = central
                              image of the source; size = source diameter,
                              independent of defocus.

Output: cassegrain_central_spot_cases.png + scaling-law table on stdout.
"""

import numpy as np
import matplotlib.pyplot as plt

import cassegrain_defocus_sim as optics
from cassegrain_extended_source_sim import disk_spectrum, LOD_ARCSEC

W_MAIN   = 16.0
W_GHOST  = 5.0      # ghost's own defocus (different optical path)
GHOST_F  = 0.12     # ghost energy fraction (exaggerated ~3x for visibility)
W040_SA  = 1.5      # strong undercorrection [waves]
LEAK_R_FRAC = 0.65  # leak disk radius as fraction of the hole radius
SRC_LOD  = 8.0      # resolved-source diameter [lambda/D]


def norm_energy(psf):
    return psf / psf.sum()


def main():
    m, r2 = optics.make_pupil(optics.n_pupil_for(W_MAIN))
    psf_main = optics.compute_psf(m, r2, W_MAIN)
    n = psf_main.shape[0]
    half = int(1.4 * optics.blur_half_px(W_MAIN))
    hole_px = optics.OBSTRUCTION * optics.blur_half_px(W_MAIN)

    cases = []

    # (a) ghost: incoherent sum of the main donut and a second defocused image
    ghost = optics.compute_psf(m, r2, W_GHOST)
    img_a = norm_energy(psf_main) + GHOST_F * norm_energy(ghost)
    cases.append((img_a, f"(a) ghost reflection\n(second image at "
                  f"W020={W_GHOST:g}, {GHOST_F*100:.0f}% energy)",
                  "size tracks focus differently than the donut"))

    # (b) strong spherical aberration, extra-focal side (best-focus referenced)
    img_b = optics.compute_psf(m, r2, W_MAIN - W040_SA, W040_SA)
    cases.append((img_b, f"(b) spherical aberration\n(W040 = {W040_SA:g} waves, "
                  "outside focus)", "appears on ONE side of focus; size ∝ donut"))

    # (c) baffle leak: additive uniform disk, fixed angular size
    c = n // 2
    ax0 = np.arange(n, dtype=np.float32) - c
    rr = np.hypot(ax0[None, :], ax0[:, None])
    leak = np.clip(LEAK_R_FRAC * hole_px - rr + 0.5, 0.0, 1.0)
    img_c = psf_main / psf_main.max() + 0.30 * leak
    cases.append((img_c, "(c) baffle / stray-light leak",
                  "size FIXED as focus changes; stays at field centre"))

    # (d) resolved source: donut + central Arago image of the source
    spec = np.fft.rfft2(np.fft.ifftshift(psf_main)) * disk_spectrum(n, SRC_LOD)
    img_d = np.fft.fftshift(np.fft.irfft2(spec, s=(n, n)))
    np.clip(img_d, 0.0, None, out=img_d)
    cases.append((img_d, f"(d) resolved source, ⌀ {SRC_LOD:g} λ/D "
                  f"({SRC_LOD*LOD_ARCSEC:.1f}″)\n(central Arago image)",
                  "size = source diameter, independent of defocus"))

    fig, axes = plt.subplots(1, 4, figsize=(12.8, 4.0))
    for ax, (img, title, sig) in zip(axes, cases):
        crop = optics.crop_center(img, half)
        ax.imshow((crop / crop.max()) ** 0.5, cmap="inferno", origin="lower",
                  vmin=0.0, vmax=1.0)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel(sig, fontsize=7.5)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"A large circular spot inside the donut hole: four causes "
                 f"(W020 = {W_MAIN:g} waves)", fontsize=11)
    fig.tight_layout()
    fig.savefig(optics.FIG_DIR / "cassegrain_central_spot_cases.png", dpi=150)
    print("saved cassegrain_central_spot_cases.png")

    print(f"\nhole radius here: {hole_px * optics.PIX_UM:.0f} um "
          f"({optics.OBSTRUCTION * 8 * W_MAIN:.0f} lambda/D)")
    print("\nhow the central spot scales when you rack the focuser:")
    print("  ghost           : both disks change size, at DIFFERENT rates")
    print("  spherical aberr.: spot only on one side of focus, size ~ donut")
    print("  baffle leak     : spot size constant, indifferent to focus")
    print("  resolved source : spot size constant = source angular diameter")
    print("and when you move the star across the field:")
    print("  ghost moves differently (mirrored/2x); leak stays at field centre;")
    print("  SA blob and source image stay concentric with the donut.")


if __name__ == "__main__":
    main()
