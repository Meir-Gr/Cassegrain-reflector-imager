"""
Interactive Cassegrain defocus simulator.

Left panel : system parameters (applied with the button) + live sliders.
Right panel: pupil drawing and the resulting image of the source.

Image model (same physics as the study scripts):
    image = SUM_lambda |FFT(pupil_lambda * defocus phase)|^2   (flat spectrum,
            pupil rescaled per sub-wavelength -> common angular grid)
            (x)  source disk  (x)  jitter  (x)  linear motion blur
applied in the frequency domain (uniform-disk visibility 2J1(x)/x,
Gaussian jitter MTF exp(-2 pi^2 sigma^2 q^2), drift MTF sinc(L q_par)).
Mirrors are achromatic: one physical defocus for all wavelengths;
W020(lambda) scales as 1/lambda, diffraction scales as lambda.

Run:  python cassegrain_gui.py            (opens the window)
      python cassegrain_gui.py --selftest (headless compute check)
"""

import math
import sys
import time
import tkinter as tk
from tkinter import ttk

import numpy as np
from scipy import fft as sfft
from scipy.special import j1
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

if __package__:
    from . import cassegrain_spectrum as spectrum
else:
    import cassegrain_spectrum as spectrum

PAD = 2                # samples per lambda_min F#
MAX_GRID = 3072
# Pupil-sampling tiers (defocus phase steps <= pi/2 per sample: n >= 16|W|).
TIERS = ((32.0, 512), (64.0, 1024))
BASE_N = TIERS[0][1]
W_SLIDER_MAX = 64.0

DEFAULTS = dict(d_mm=200.0, obstruction=0.35, vane_mm=4.0,
                angles_deg=(90.0, 210.0, 330.0), lambda_nm=550.0, fnum=10.0)


class OpticalModel:
    """Pupil + image computation for one set of system parameters."""

    def __init__(self, **kw):
        self.set_params(**kw)

    def set_params(self, d_mm, obstruction, vane_mm, angles_deg,
                   lambda_nm, fnum):
        if d_mm <= 0 or lambda_nm <= 0 or fnum <= 0:
            raise ValueError("D, wavelength and F# must be positive")
        if not 0.0 <= obstruction <= 0.95:
            raise ValueError("obstruction must be in 0 .. 0.95")
        if vane_mm < 0:
            raise ValueError("vane width cannot be negative")
        self.D = d_mm * 1e-3
        self.eps = obstruction
        self.wv = vane_mm * 1e-3
        self.angles = tuple(float(a) for a in angles_deg)
        self.lam = lambda_nm * 1e-9
        self.fnum = fnum
        self._trig = [(np.float32(np.cos(a)), np.float32(np.sin(a)))
                      for a in np.deg2rad(self.angles)]
        self._tiers = {}
        self._spec_key = self._spec = None
        self._preview_key = self._preview = None

    @staticmethod
    def n_pupil_for(w):
        for lim, n in TIERS:
            if abs(w) <= lim:
                return n
        return TIERS[-1][1]

    def _tier(self, n_pupil, pad=PAD):
        """Lazily built coordinate / frequency grids for one tier."""
        t = self._tiers.get(n_pupil)
        if t is not None and t["pad"] == pad:
            return t
        n = 64 * math.ceil(n_pupil * pad / 64)
        if n > MAX_GRID:
            raise ValueError(f"Band and defocus need a {n}² grid (limit {MAX_GRID}); "
                             "narrow the band or reduce defocus")
        dx = self.D * pad / n     # rounding the FFT size must not change image pitch
        c = ((np.arange(n) - n // 2) * dx).astype(np.float32)
        fy = np.fft.fftfreq(n)[:, None].astype(np.float32)
        fx = np.fft.rfftfreq(n)[None, :].astype(np.float32)
        q = np.hypot(fy, fx)                     # cycles/px
        t = dict(n=n, pad=pad, coords=c, q=q, q2=q * q, fy=fy, fx=fx)
        self._tiers[n_pupil] = t
        return t

    def _mask_rho2(self, coords, scale):
        """Aperture mask and rho^2 for a pupil scaled by lam0/lam_i.

        Enlarging the pupil in grid samples by `scale` keeps the angular
        image pixel identical for every sub-wavelength (no interpolation).
        """
        x, y = coords[None, :], coords[:, None]
        rr2 = x * x + y * y
        R = np.float32(self.D / 2.0 * scale)
        mask = (rr2 <= R * R) & (rr2 >= (self.eps * R) ** 2)
        if self.wv > 0:
            hw = np.float32(self.wv / 2.0 * scale)
            for ca, sa in self._trig:
                along = x * ca + y * sa
                across = y * ca - x * sa
                mask &= ~((along >= 0) & (np.abs(across) <= hw))
        return mask.astype(np.float32), rr2 * np.float32(1.0 / (R * R))

    n_lambda = staticmethod(spectrum.n_lambda)

    # ---- derived scales -------------------------------------------------
    @property
    def lod_arcsec(self):
        return np.degrees(self.lam / self.D) * 3600.0

    @property
    def pix_um(self):
        """Monochromatic image pitch at the reference wavelength."""
        return self.lam * self.fnum / PAD * 1e6

    def image_pixel_um(self, bandwidth):
        return self.pix_um * spectrum.short_rel(bandwidth)

    @property
    def dz_per_wave_mm(self):
        return 8.0 * self.lam * self.fnum ** 2 * 1e3

    # ---- computation ----------------------------------------------------
    def _sampling(self, w020, bandwidth):
        nl = self.n_lambda(w020, bandwidth)
        return self._tier(self.n_pupil_for(w020), PAD / spectrum.short_rel(bandwidth)), nl

    def _pupils(self, t, w020, bandwidth, nl):
        """Yield each transmitting pupil in its bounding box, and its power."""
        n = t["n"]
        for k in range(nl):
            rel = 1.0 - bandwidth / 2.0 + bandwidth * (k + 0.5) / nl
            radius = math.ceil(n / (2 * t["pad"] * rel)) + 1  # include edge roundoff
            coords = t["coords"][n//2 - radius:n//2 + radius + 1]
            mask, rho2 = self._mask_rho2(coords, 1.0 / rel)
            power = mask.sum()
            if power == 0:
                raise ValueError("No transmitting pupil samples; check the obstruction and vanes")
            if w020 == 0:
                yield mask, power
            else:
                phase = np.float32(2.0 * np.pi * w020 / rel) * rho2
                yield mask * np.exp(1j * phase), power

    def _psf_spectrum(self, w020, bandwidth):
        t, nl = self._sampling(w020, bandwidth)
        key = (w020, bandwidth)
        if self._spec_key != key:
            acc = np.zeros((t["n"], t["n"]), np.float32)
            for pupil, power in self._pupils(t, w020, bandwidth, nl):
                # Zero-pad to the original grid. Translating a pupil changes
                # Fourier phase, not intensity; keep the PSF in native FFT order.
                field = sfft.fft2(pupil, s=acc.shape, workers=-1)
                acc += np.abs(field)**2 / power
            self._spec = sfft.rfft2(acc, workers=-1)
            self._spec_key = key
        return self._spec, t

    def image(self, w020, src_lod, jitter_lod, motion_lod=0.0, motion_deg=0.0,
              bandwidth=0.0):
        F, t = self._psf_spectrum(w020, bandwidth)
        if src_lod > 0:
            x = np.float32(np.pi * src_lod * t["pad"]) * t["q"]
            S = np.ones_like(x)
            nz = x > 1e-12
            S[nz] = 2.0 * j1(x[nz]) / x[nz]
            F = F * S
        if jitter_lod > 0:
            F = F * np.exp(np.float32(-2.0 * np.pi**2 * (jitter_lod * t["pad"])**2)
                           * t["q2"])
        if motion_lod > 0:
            ph = np.deg2rad(motion_deg)
            q_par = (t["fx"] * np.float32(np.cos(ph))
                     + t["fy"] * np.float32(np.sin(ph)))
            F = F * np.sinc(np.float32(motion_lod * t["pad"]) * q_par)
        img = np.fft.fftshift(sfft.irfft2(F, s=(t["n"], t["n"]), workers=-1))
        np.clip(img, 0.0, None, out=img)
        m = img.max()
        return img / m if m > 0 else img

    def crop_half_px(self, w020, src_lod, jitter_lod, motion_lod, n_grid, bandwidth=0.0):
        pad = PAD / spectrum.short_rel(bandwidth)
        geom = 4.0 * abs(w020) * pad
        extra = pad * (src_lod / 2.0 + 3.0 * jitter_lod + motion_lod / 2.0)
        return int(min(n_grid // 2 - 2, max(80, 1.5 * geom + extra + 20)))

    def _point_crop(self, t, w020, bandwidth, nl, half):
        """Compute just the central DFT bins, or request the full field if its peak is uncertain."""
        n = t["n"]
        bins = np.arange(-half, half)  # fftshift order without shifting a full array
        outside = slice(half, n - half)
        acc = np.zeros((2 * half, 2 * half), np.float32)
        outside_bound = 0.0
        for pupil, power in self._pupils(t, w020, bandwidth, nl):
            rows = sfft.fft(pupil, n=n, axis=1, workers=-1)
            cols = sfft.fft(rows[:, bins], n=n, axis=0, workers=-1)
            acc += np.abs(cols[bins, :])**2 / power
            if w020 != 0:  # in focus, the peak of a nonnegative pupil is at DC
                # In omitted columns, every y-FFT amplitude is <= sum_y |row FFT|.
                # In retained columns, the omitted y bins are already known.
                bx = np.abs(rows[:, outside]).sum(axis=0, dtype=np.float64).max()
                by = np.abs(cols[outside, :]).max()
                outside_bound += float(max(bx, by))**2 / float(power)
        peak = acc.max()
        if outside_bound >= peak:
            return None  # do not normalize to a local maximum
        return acc / peak

    def preview_image(self, w020, src_lod, jitter_lod, motion_lod=0.0, motion_deg=0.0,
                      bandwidth=0.0):
        """Display crop and full grid size, with the same scale and normalization as image()."""
        t, nl = self._sampling(w020, bandwidth)
        n = t["n"]
        half = self.crop_half_px(w020, src_lod, jitter_lod, motion_lod, n, bandwidth)
        point = src_lod == jitter_lod == motion_lod == 0
        key = (w020, bandwidth, half)
        if point and self._preview_key == key:
            return self._preview, n

        crop = None
        if point and self._spec_key != (w020, bandwidth):
            crop = self._point_crop(t, w020, bandwidth, nl, half)
        if crop is None:
            # Blur on the whole field: flux outside the viewport can enter it.
            img = self.image(w020, src_lod, jitter_lod, motion_lod, motion_deg, bandwidth)
            view = slice(n//2 - half, n//2 + half)
            crop = img[view, view]
            if point:
                crop = crop.copy()  # do not pin the full image in the preview cache
        if point:
            crop.setflags(write=False)
            self._preview_key, self._preview = key, crop
        return crop, n

    def pupil_view(self):
        """(mask crop, half-extent in mm) for display, at the base tier."""
        t = self._tier(BASE_N)
        half = BASE_N // 2 + 24
        c = t["n"] // 2
        crop, _ = self._mask_rho2(t["coords"][c - half:c + half], 1.0)
        return crop, half * (self.D / BASE_N) * 1e3


class App:
    def __init__(self, root):
        self.root = root
        root.title("Cassegrain defocus simulator — 3-vane spider")
        self.model = OpticalModel(**DEFAULTS)
        self.bandwidth = 0.0
        self._job = None

        panel = ttk.Frame(root, padding=8)
        panel.grid(row=0, column=0, sticky="ns")
        self._build_entries(panel)
        self._build_sliders(panel)

        self.fig = Figure(figsize=(9.2, 4.7), dpi=100)
        self.ax_pupil = self.fig.add_subplot(1, 2, 1)
        self.ax_img = self.fig.add_subplot(1, 2, 2)
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=0, column=1, sticky="nsew")
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)

        self.status = ttk.Label(root, anchor="w")
        self.status.grid(row=1, column=0, columnspan=2, sticky="ew", padx=8)

        self._draw_pupil()
        self._update_image()

    # ---- widgets ---------------------------------------------------------
    def _build_entries(self, parent):
        box = ttk.LabelFrame(parent, text="System parameters", padding=6)
        box.pack(fill="x")
        self.entries = {}
        rows = [("Primary diameter D [mm]", "d_mm", DEFAULTS["d_mm"]),
                ("Obstruction ratio d/D", "obstruction", DEFAULTS["obstruction"]),
                ("Vane width [mm]", "vane_mm", DEFAULTS["vane_mm"]),
                ("Vane angle 1 [deg]", "a1", DEFAULTS["angles_deg"][0]),
                ("Vane angle 2 [deg]", "a2", DEFAULTS["angles_deg"][1]),
                ("Vane angle 3 [deg]", "a3", DEFAULTS["angles_deg"][2]),
                ("Lower wavelength [µm]", "lower_um", DEFAULTS["lambda_nm"] / 1000),
                ("Upper wavelength [µm]", "upper_um", DEFAULTS["lambda_nm"] / 1000),
                ("Focal ratio F#", "fnum", DEFAULTS["fnum"])]
        for i, (label, key, val) in enumerate(rows):
            ttk.Label(box, text=label).grid(row=i, column=0, sticky="w")
            e = ttk.Entry(box, width=9, justify="right")
            e.insert(0, f"{val:g}")
            e.grid(row=i, column=1, padx=4, pady=1)
            self.entries[key] = e
        ttk.Button(box, text="Apply system parameters",
                   command=self._apply_params).grid(
            row=len(rows), column=0, columnspan=2, pady=(6, 2), sticky="ew")
        self.err = ttk.Label(box, foreground="red", wraplength=210)
        self.err.grid(row=len(rows) + 1, column=0, columnspan=2, sticky="w")

    def _build_sliders(self, parent):
        box = ttk.LabelFrame(parent, text="Live controls", padding=6)
        box.pack(fill="x", pady=(8, 0))
        self.sliders, self.slabels = {}, {}
        specs = [("defocus", "Defocus W020 [waves]",
                  -W_SLIDER_MAX, W_SLIDER_MAX, 0.5, 8.0),
                 ("source", "Source diameter [λ/D]", 0.0, 16.0, 0.25, 0.0),
                 ("jitter", "Jitter RMS σ [λ/D]", 0.0, 4.0, 0.1, 0.0),
                 ("motion", "Motion blur per integration [µm on sensor]",
                  0.0, 300.0, 1.0, 0.0),
                 ("motion_dir", "Motion direction [deg]", 0.0, 180.0, 5.0, 0.0)]
        for key, label, lo, hi, res, init in specs:
            ttk.Label(box, text=label).pack(anchor="w", pady=(4, 0))
            s = tk.Scale(box, from_=lo, to=hi, resolution=res,
                         orient="horizontal", length=230, showvalue=False,
                         command=lambda _v, k=key: self._on_slider(k))
            s.set(init)
            s.pack(fill="x")
            self.sliders[key] = s
            lab = ttk.Label(box, text="")
            lab.pack(anchor="w")
            self.slabels[key] = lab
        self.log_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(box, text="Log stretch (in-focus detail)",
                        variable=self.log_var,
                        command=lambda: self._on_slider("defocus")).pack(
            anchor="w", pady=(6, 0))

    # ---- callbacks --------------------------------------------------------
    def _apply_params(self):
        try:
            g = lambda k: float(self.entries[k].get())
            lam, band = spectrum.from_um(self.entries["lower_um"].get(),
                                         self.entries["upper_um"].get())
            self.model.set_params(
                d_mm=g("d_mm"), obstruction=g("obstruction"),
                vane_mm=g("vane_mm"),
                angles_deg=(g("a1"), g("a2"), g("a3")),
                lambda_nm=lam, fnum=g("fnum"))
            self.bandwidth = band
            self.err.config(text="")
        except ValueError as exc:
            self.err.config(text=str(exc))
            return
        self._draw_pupil()
        self._update_image()

    def _on_slider(self, _key):
        if self._job is None:
            self._job = self.root.after(25, self._update_image)

    # ---- drawing ----------------------------------------------------------
    def _draw_pupil(self):
        m = self.model
        ax = self.ax_pupil
        ax.clear()
        crop, r_mm = m.pupil_view()
        ax.imshow(crop, cmap="gray", origin="lower",
                  extent=[-r_mm, r_mm, -r_mm, r_mm])
        ax.set_title(f"Pupil   (λ/D = {m.lod_arcsec:.2f}″)", fontsize=10)
        ax.set_xlabel("pupil [mm]", fontsize=8)
        ax.tick_params(labelsize=7)

    def _update_image(self):
        self._job = None
        m = self.model
        w = float(self.sliders["defocus"].get())
        b = self.bandwidth
        src = float(self.sliders["source"].get())
        jit = float(self.sliders["jitter"].get())
        mot_um = float(self.sliders["motion"].get())
        mot = mot_um / (m.lam * m.fnum * 1e6)    # focal-plane um -> lambda/D
        mdir = float(self.sliders["motion_dir"].get())

        t0 = time.perf_counter()
        try:
            crop, n_grid = m.preview_image(w, src, jit, mot, mdir, b)
        except ValueError as exc:
            self.ax_img.clear()
            self.ax_img.set_title("No image for these settings")
            self.canvas.draw_idle()
            self.status.config(text=str(exc))
            return
        dt = time.perf_counter() - t0

        ax = self.ax_img
        ax.clear()
        pitch = m.image_pixel_um(b)
        ext = (crop.shape[0] // 2) * pitch
        if self.log_var.get():
            ax.imshow(np.log10(crop + 1e-8), cmap="inferno", origin="lower",
                      vmin=-6, vmax=0, extent=[-ext, ext, -ext, ext])
        else:
            ax.imshow(crop ** 0.5, cmap="inferno", origin="lower",
                      vmin=0, vmax=1, extent=[-ext, ext, -ext, ext])
        donut_um = 8.0 * abs(w) * m.lam * m.fnum * 1e6
        ax.set_title(f"Image   (donut ⌀ ≈ {donut_um:.0f} µm)", fontsize=10)
        ax.set_xlabel("focal plane [µm]", fontsize=8)
        ax.tick_params(labelsize=7)
        self.canvas.draw_idle()

        self.slabels["defocus"].config(
            text=f"{w:+.1f} waves   (Δz ≈ {w * m.dz_per_wave_mm:+.1f} mm)")
        nl = m.n_lambda(w, b)
        self.slabels["source"].config(
            text=f"{src:.2f} λ/D   (= {src * m.lod_arcsec:.2f}″)")
        self.slabels["jitter"].config(
            text=f"{jit:.2f} λ/D   (= {jit * m.lod_arcsec:.2f}″ RMS)")
        self.slabels["motion"].config(
            text=f"{mot_um:.0f} µm smear   (= {mot:.2f} λ/D "
                 f"= {mot * m.lod_arcsec:.2f}″ on sky)")
        self.slabels["motion_dir"].config(
            text=f"{mdir:.0f}° from +x" if mot > 0 else "(inactive)")
        self.status.config(
            text=f"λ₀ {m.lam * 1e6:g} µm · Δλ/λ₀ {b:.2%} · Nλ {nl} · "
                 f"grid {n_grid}² · image pixel {pitch:.2f} µm · render {dt * 1e3:.0f} ms")


def selftest():
    from unittest.mock import patch

    model = OpticalModel(**DEFAULTS)
    for w, s, j, mo, md, b in [(0, 0, 0, 0, 0, 0), (8, 0, 0, 0, 0, 0),
                               (32, 4, 1, 0, 0, 0), (-16, 2, 0.5, 0, 0, 0),
                               (16, 0, 0, 8, 45, 0), (16, 0, 0, 0, 0, 0.2),
                               (64, 4, 1, 4, 90, 0.4)]:
        t0 = time.perf_counter()
        img = model.image(w, s, j, mo, md, b)
        dt = (time.perf_counter() - t0) * 1e3
        print(f"W020={w:+6.1f}  src={s:4.1f}  jit={j:3.1f}  "
              f"mot={mo:4.1f}@{md:3.0f}°  b={b:4.2f} "
              f"Nλ={model.n_lambda(w, b):2d}  grid={img.shape[0]:4d}²  "
              f"max={img.max():.2f}  sum={img.sum():.3e}  {dt:5.0f} ms")
    model.set_params(d_mm=300, obstruction=0.2, vane_mm=1.0,
                     angles_deg=(0, 120, 240), lambda_nm=650, fnum=8)
    img = model.image(10, 1, 0.2, 2, 30, 0.1)
    print(f"param change OK, image max={img.max():.2f}")

    # Independent in-focus check: integrate the analytic circular-aperture
    # Airy intensity on a physical radial grid. Unit-flux peaks scale as 1/lambda^2.
    lam, b = spectrum.from_um(3.5, 5.5)
    model.set_params(**{**DEFAULTS, "lambda_nm": lam, "obstruction": 0, "vane_mm": 0})
    img = model.image(0, 0, 0, bandwidth=b)
    pitch = model.image_pixel_um(b)
    assert math.isclose(pitch, 17.5)
    wavelengths = 3.5 + (np.arange(2048) + 0.5) * 2 / 2048
    u = np.pi * np.arange(1, 9)[:, None] * pitch / (10 * wavelengths)
    airy = np.average((2 * j1(u) / u)**2, axis=1, weights=wavelengths**-2)
    c = img.shape[0] // 2
    error = np.max(np.abs(img[c, c+1:c+9] - airy))
    assert error < 0.001, f"MWIR Airy profile error {error:.3g}"
    print(f"3.5–5.5 µm: pitch {pitch:g} µm, Airy profile max error {error:.2e} of peak")

    # Compare both fast paths against a dense, uncropped diffraction calculation.
    w = -8.0
    model.set_params(**{**DEFAULTS, "lambda_nm": lam, "angles_deg": (12, 131, 233)})
    crop, n = model.preview_image(w, 0, 0, bandwidth=b)
    t, nl = model._sampling(w, b)
    dense = np.zeros((n, n), np.float32)
    for k in range(nl):
        rel = 1 - b/2 + b*(k + 0.5)/nl
        mask, rho2 = model._mask_rho2(t["coords"], 1/rel)
        pupil = mask * np.exp(1j * np.float32(2*np.pi*w/rel) * rho2)
        dense += np.abs(sfft.fft2(pupil, workers=-1))**2 / mask.sum()
    dense = np.fft.fftshift(dense / dense.max())
    half = model.crop_half_px(w, 0, 0, 0, n, b)
    view = slice(n//2 - half, n//2 + half)
    np.testing.assert_allclose(crop, dense[view, view], atol=3e-6, rtol=0)
    np.testing.assert_allclose(model.image(w, 0, 0, bandwidth=b), dense, atol=3e-6, rtol=0)
    for src, jit, mot in [(2, 0, 0), (0, 0.5, 0), (0, 0, 4)]:
        crop, n = model.preview_image(w, src, jit, mot, 35, b)
        half = model.crop_half_px(w, src, jit, mot, n, b)
        view = slice(n//2 - half, n//2 + half)
        full = model.image(w, src, jit, mot, 35, b)
        np.testing.assert_allclose(crop, full[view, view], atol=3e-6, rtol=0)
    with patch.object(model, "crop_half_px", return_value=2):
        crop, n = model.preview_image(-w, 0, 0, bandwidth=b)
        full = model.image(-w, 0, 0, bandwidth=b)
        view = slice(n//2 - 2, n//2 + 2)
        assert full[view, view].max() < 1  # the global peak is outside this window
        np.testing.assert_allclose(crop, full[view, view], atol=3e-6, rtol=0)
    print("MWIR dense FFT, point/blur previews and off-screen peak checks passed")
    for lower, upper in [(5.5, 3.5), (0, 5.5), (3.5, np.inf), (np.nan, 5.5)]:
        np.testing.assert_raises(ValueError, spectrum.from_um, lower, upper)
    np.testing.assert_raises_regex(ValueError, "grid", model.image, 0, 0, 0, bandwidth=1.5)
    np.testing.assert_raises_regex(ValueError, "wavelength samples", model.image, 64, 0, 0, bandwidth=1)
    print("selftest passed")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        try:  # crisp text on Windows high-DPI displays
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
        root = tk.Tk()
        App(root)
        root.mainloop()
