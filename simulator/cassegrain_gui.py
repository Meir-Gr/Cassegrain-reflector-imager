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

import sys
import time
import tkinter as tk
from tkinter import ttk

import numpy as np
from scipy import fft as sfft
from scipy.special import j1
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

PAD = 2                # image pixel = lambda*F#/2 (Nyquist)
# pupil-sampling tiers (defocus phase steps <= pi/2 per sample: n >= 16|W|);
# capped at 2048^2 grids to keep slider response fluid.
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
        self._spec_key = None

    @staticmethod
    def n_pupil_for(w):
        for lim, n in TIERS:
            if abs(w) <= lim:
                return n
        return TIERS[-1][1]

    def _tier(self, n_pupil):
        """Lazily built coordinate / frequency grids for one tier."""
        t = self._tiers.get(n_pupil)
        if t is not None:
            return t
        n = n_pupil * PAD
        dx = self.D / n_pupil
        c = ((np.arange(n) - n // 2) * dx).astype(np.float32)
        X, Y = np.meshgrid(c, c)
        fy = np.fft.fftfreq(n)[:, None].astype(np.float32)
        fx = np.fft.rfftfreq(n)[None, :].astype(np.float32)
        q = np.hypot(fy, fx)                     # cycles/px
        t = dict(n=n, X=X, Y=Y, rr2=X * X + Y * Y,
                 q=q, q2=q * q, fy=fy, fx=fx)
        self._tiers[n_pupil] = t
        return t

    def _mask_rho2(self, t, scale):
        """Aperture mask and rho^2 for a pupil scaled by lam0/lam_i.

        Enlarging the pupil in grid samples by `scale` keeps the angular
        image pixel identical for every sub-wavelength (no interpolation).
        """
        R = np.float32(self.D / 2.0 * scale)
        mask = (t["rr2"] <= R * R) & (t["rr2"] >= (self.eps * R) ** 2)
        if self.wv > 0:
            hw = np.float32(self.wv / 2.0 * scale)
            for ca, sa in self._trig:
                along = t["X"] * ca + t["Y"] * sa
                across = t["Y"] * ca - t["X"] * sa
                mask &= ~((along >= 0) & (np.abs(across) <= hw))
        return mask.astype(np.float32), t["rr2"] * np.float32(1.0 / (R * R))

    @staticmethod
    def n_lambda(w020, bandwidth):
        """Spectral samples so adjacent sub-images shift < ~1 ring width."""
        if bandwidth <= 0:
            return 1
        return int(np.clip(3 + 20.0 * bandwidth
                           + 0.35 * bandwidth * abs(w020) * PAD, 3, 15))

    # ---- derived scales -------------------------------------------------
    @property
    def lod_arcsec(self):
        return np.degrees(self.lam / self.D) * 3600.0

    @property
    def pix_um(self):
        return self.lam * self.fnum / PAD * 1e6

    @property
    def dz_per_wave_mm(self):
        return 8.0 * self.lam * self.fnum ** 2 * 1e3

    # ---- computation ----------------------------------------------------
    def _psf_spectrum(self, w020, bandwidth):
        n_pupil = self.n_pupil_for(w020)     # criterion is lambda-invariant
        key = (w020, round(bandwidth, 4), n_pupil)
        if self._spec_key != key:
            t = self._tier(n_pupil)
            nl = self.n_lambda(w020, bandwidth)
            acc = None
            for k in range(nl):
                rel = (1.0 if nl == 1 else
                       1.0 - bandwidth / 2.0 + bandwidth * (k + 0.5) / nl)
                mask, rho2 = self._mask_rho2(t, 1.0 / rel)
                phase = np.float32(2.0 * np.pi * w020 / rel) * rho2
                pupil = mask * np.exp(1j * phase)
                field = np.fft.fftshift(sfft.fft2(np.fft.ifftshift(pupil),
                                                  workers=-1))
                psf = np.abs(field).astype(np.float32) ** 2
                psf /= mask.sum()                # equal energy per wavelength
                acc = psf if acc is None else acc + psf
            self._spec = sfft.rfft2(np.fft.ifftshift(acc), workers=-1)
            self._spec_key = key
        return self._spec, self._tier(self.n_pupil_for(w020))

    def image(self, w020, src_lod, jitter_lod, motion_lod=0.0, motion_deg=0.0,
              bandwidth=0.0):
        F, t = self._psf_spectrum(w020, bandwidth)
        if src_lod > 0:
            x = np.float32(np.pi * src_lod * PAD) * t["q"]
            S = np.ones_like(x)
            nz = x > 1e-12
            S[nz] = 2.0 * j1(x[nz]) / x[nz]
            F = F * S
        if jitter_lod > 0:
            F = F * np.exp(np.float32(-2.0 * np.pi**2 * (jitter_lod * PAD)**2)
                           * t["q2"])
        if motion_lod > 0:
            ph = np.deg2rad(motion_deg)
            q_par = (t["fx"] * np.float32(np.cos(ph))
                     + t["fy"] * np.float32(np.sin(ph)))
            F = F * np.sinc(np.float32(motion_lod * PAD) * q_par)
        img = np.fft.fftshift(sfft.irfft2(F, s=(t["n"], t["n"]), workers=-1))
        np.clip(img, 0.0, None, out=img)
        m = img.max()
        return img / m if m > 0 else img

    def crop_half_px(self, w020, src_lod, jitter_lod, motion_lod, n_grid):
        geom = 4.0 * abs(w020) * PAD
        extra = PAD * (src_lod / 2.0 + 3.0 * jitter_lod + motion_lod / 2.0)
        return int(min(n_grid // 2 - 2, max(80, 1.5 * geom + extra + 20)))

    def pupil_view(self):
        """(mask crop, half-extent in mm) for display, at the base tier."""
        t = self._tier(BASE_N)
        mask, _ = self._mask_rho2(t, 1.0)
        half = BASE_N // 2 + 24
        c = mask.shape[0] // 2
        crop = mask[c - half:c + half, c - half:c + half]
        return crop, half * (self.D / BASE_N) * 1e3


class App:
    def __init__(self, root):
        self.root = root
        root.title("Cassegrain defocus simulator — 3-vane spider")
        self.model = OpticalModel(**DEFAULTS)
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
                ("Wavelength [nm]", "lambda_nm", DEFAULTS["lambda_nm"]),
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
                 ("bandwidth", "Bandwidth Δλ/λ (flat spectrum)",
                  0.0, 0.4, 0.01, 0.0),
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
            self.model.set_params(
                d_mm=g("d_mm"), obstruction=g("obstruction"),
                vane_mm=g("vane_mm"),
                angles_deg=(g("a1"), g("a2"), g("a3")),
                lambda_nm=g("lambda_nm"), fnum=g("fnum"))
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
        b = float(self.sliders["bandwidth"].get())
        src = float(self.sliders["source"].get())
        jit = float(self.sliders["jitter"].get())
        mot_um = float(self.sliders["motion"].get())
        mot = mot_um / (m.lam * m.fnum * 1e6)    # focal-plane um -> lambda/D
        mdir = float(self.sliders["motion_dir"].get())

        t0 = time.perf_counter()
        img = m.image(w, src, jit, mot, mdir, b)
        n_grid = img.shape[0]
        half = m.crop_half_px(w, src, jit, mot, n_grid)
        c = n_grid // 2
        crop = img[c - half:c + half, c - half:c + half]
        dt = time.perf_counter() - t0

        ax = self.ax_img
        ax.clear()
        ext = half * m.pix_um
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
        lam0 = m.lam * 1e9
        nl = m.n_lambda(w, b)
        self.slabels["bandwidth"].config(
            text=(f"Δλ ≈ {b * lam0:.0f} nm   "
                  f"({lam0 * (1 - b / 2):.0f}–{lam0 * (1 + b / 2):.0f} nm, "
                  f"Nλ = {nl})") if b > 0 else "(monochromatic)")
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
            text=f"grid {n_grid}² · Nλ {nl} · render {dt * 1e3:.0f} ms · "
                 f"image pixel {m.pix_um:.2f} µm")


def selftest():
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
