"""
Interactive front end of the Cassegrain video simulator (cassegrain_video.py).

Left : parameters in tabs (Optics, Camera, Noise, Motion, Output) with live
       readouts (scales, speed in all units, blur, photometry); buttons to
       run the whole sequence in the background (saved to output/<name>/ as
       .npy / .tif / .mat / raw binary, chosen in the Output tab), cancel it,
       and save / load the JSON configuration that the command line also uses.
Right: the full sensor frame with the true track and the current exposure
       streak, a zoom on the spot, an exact pixel view (every camera pixel
       drawn as a uniform k x k block of screen pixels, 1:1 up to 32x, with
       pixel grid, scale bar and the value of the pixel under the mouse),
       and playback (frame slider, play, display stretch).

Every edit re-renders the frame under the slider. Noise is seeded per frame,
so a preview is identical to the same frame of a run.

Run:  python cassegrain_video_gui.py [config.json]
"""

import queue
import sys
import threading
import time
import tkinter as tk
from dataclasses import fields
from pathlib import Path
from tkinter import filedialog, ttk

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk

import cassegrain_video as cv

# choice fields: (stored values, labels shown)
CHOICES = {
    "speed_unit": (cv.SPEED_UNITS, ("px/s", "µm/s", "arcsec/s", "deg/s")),
    "intensity_mode": (cv.INTENSITY_MODES,
                       ("Static peak [% of full scale]", "Photoelectrons [e⁻/s]")),
    "preview": (cv.PREVIEW_KINDS, ("MP4 video", "Animated GIF", "none")),
    "raw": (cv.RAW_KINDS, ("none", "one file per frame (raw/)", "one file, all frames")),
}

# (field, label, kind): f float, i int, b bool, t three floats, c choice, s text
TABS = (
    ("Optics", (
        ("d_mm", "Primary diameter D [mm]", "f"),
        ("obstruction", "Obstruction ratio d/D", "f"),
        ("vane_mm", "Vane width [mm]", "f"),
        ("vane_angles_deg", "Vane angles [deg]", "t"),
        ("lower_um", "Lower wavelength [µm]", "f"),
        ("upper_um", "Upper wavelength [µm]", "f"),
        ("fnum", "Focal ratio F#", "f"),
        ("defocus_waves", "Defocus W020 [waves]", "f"),
        ("source_lod", "Source diameter [λ/D]", "f"),
        ("jitter_lod", "In-frame jitter RMS [λ/D]", "f"))),
    ("Camera", (
        ("width", "Image width W [px]", "i"),
        ("height", "Image height H [px]", "i"),
        ("pixel_um", "Pixel size [µm]", "f"),
        ("frames", "Sequence length [frames]", "i"),
        ("fps", "Frame rate [fps]", "f"),
        ("exposure_ms", "Exposure [ms]", "f"),
        ("rolling_shutter", "Rolling shutter", "b"),
        ("line_time_us", "Line time [µs] (rolling)", "f"))),
    ("Noise", (
        ("bit_depth", "Bit depth", "i"),
        ("gain_e_per_dn", "Gain [e⁻/DN]", "f"),
        ("offset_dn", "Black level [DN]", "f"),
        ("full_well_e", "Full well [e⁻]", "f"),
        ("read_noise_e", "Read noise [e⁻ RMS]", "f"),
        ("dark_e_per_s", "Dark current [e⁻/s/px]", "f"),
        ("background_e_per_s", "Background [e⁻/s/px]", "f"),
        ("shot_noise", "Photon shot noise", "b"),
        ("seed", "Noise seed", "i"))),
    ("Motion", (
        ("start_x", "Start x (column) [px]", "f"),
        ("start_y", "Start y (row) [px]", "f"),
        ("speed", "Speed", "f"),
        ("speed_unit", "Speed unit", "c"),
        ("direction_deg", "Direction [deg] (0 right, 90 up)", "f"),
        ("frame_jitter_px", "Frame-to-frame jitter [px RMS]", "f"),
        ("intensity_mode", "Intensity given as", "c"),
        ("intensity", "Intensity value", "f"))),
    ("Output", (
        ("name", "Run name (blank: time stamp)", "s"),
        ("folder", "Parent folder (blank: output/)", "s"),
        ("npy", "frames.npy (lossless DN)", "b"),
        ("tiff", "frames.tif (multi-page)", "b"),
        ("mat", "frames.mat (MATLAB, H×W×N)", "b"),
        ("raw", "Raw binary frames (.raw)", "c"),
        ("preview", "Preview video", "c"),
        ("truth", "truth.csv (true positions)", "b"))),
)
STRETCH = ("sqrt", "linear", "log")
SIM_KEYS = {f.name for f in fields(cv.VideoParams)}
ZOOM_STEPS = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32)       # pixel-view magnifications
PIXEL_ZOOMS = ("fit spot",) + tuple("1:1" if z == 1 else f"{z}×" for z in ZOOM_STEPS)
GRID_MIN_ZOOM = 6                                   # pixel grid drawn from 6× up
GRID_RGB = np.array((70, 120, 190), np.uint16)      # mixed 50 % into grid lines
OUTSIDE_RGB = (38, 44, 66)                          # pixel view beyond the sensor
PV_HINT = "wheel: zoom · drag: pan · double-click: re-centre"


def spot_extent_px(p):
    """Side [px] of a box holding the spot: donut + motion blur + source,
    plus a margin of 3 lambda/D + 8 px on each side."""
    lod = p.lf_um / p.pixel_um
    return (8 * abs(p.defocus_waves) * lod + p.blur_px + p.source_lod * lod
            + 6 * lod + 16)


class App:
    def __init__(self, root, config=None):
        self.root = root
        root.title("Cassegrain video simulator — moving spot on a real camera")
        self.vars = {}
        self.q = queue.Queue()
        self.cancel = threading.Event()
        self.worker = None            # the single background thread
        self.pending = None           # frame index waiting to be rendered
        self.run_request = None       # (params, output) waiting to run
        self.params = self.output = None
        self.cur = None               # (k, dn, sim) on screen
        self.last_sim = None          # simulator of the latest preview
        self.run = None               # last run: frames memmap, params, sim, folder
        self.playing = False
        self._edit_job = None
        self._loading = False
        self._track = (None, None)
        self._gen = 0                 # bumped on every parameter change
        self._want = None             # (frame, generation) last requested
        self._view = None             # (stretched frame, dn, sim, truth) on screen
        self._pv_map = None           # pixel view: (x0, y0, zoom, dn, params)
        self._pv_offset = (0, 0)      # pixel view panning [px] from the spot
        self._pv_drag = None          # (x, y, offset) at the start of a drag
        self._pv_mouse = None         # mouse position over the pixel view
        self._pv_photo = None         # Tk image (must stay referenced)
        self._pv_rgb = None           # the RGB array shown in the pixel view

        left = ttk.Frame(root, padding=6)
        left.grid(row=0, column=0, sticky="ns")
        nb = ttk.Notebook(left)
        nb.pack(fill="x")
        for title, specs in TABS:
            self._build_tab(nb, title, specs)
        self.readout = ttk.Label(left, justify="left", font="TkFixedFont",
                                 padding=(2, 6))
        self.readout.pack(fill="x")
        row = ttk.Frame(left)
        row.pack(fill="x", pady=(2, 0))
        self.b_run = ttk.Button(row, text="Run sequence", command=self.start_run)
        self.b_run.pack(side="left", expand=True, fill="x")
        self.b_cancel = ttk.Button(row, text="Cancel", state="disabled",
                                   command=self.cancel.set)
        self.b_cancel.pack(side="left", padx=(4, 0))
        self.progress = ttk.Progressbar(left, maximum=1.0)
        self.progress.pack(fill="x", pady=(4, 0))
        row = ttk.Frame(left)
        row.pack(fill="x", pady=(6, 0))
        ttk.Button(row, text="Load config…", command=self.load_cfg).pack(
            side="left", expand=True, fill="x")
        ttk.Button(row, text="Save config…", command=self.save_cfg).pack(
            side="left", expand=True, fill="x", padx=(4, 0))
        self.err = ttk.Label(left, foreground="red", wraplength=330, justify="left")
        self.err.pack(fill="x", pady=(4, 0))

        right = ttk.Frame(root)
        right.grid(row=0, column=1, sticky="nsew")
        self.fig = Figure(figsize=(7.8, 4.8), dpi=100)
        gs = self.fig.add_gridspec(1, 2, width_ratios=(1.75, 1))
        self.ax = self.fig.add_subplot(gs[0])
        self.axz = self.fig.add_subplot(gs[1])
        self.fig.subplots_adjust(left=0.06, right=0.98, bottom=0.07, top=0.93, wspace=0.2)
        self._init_artists()
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self._build_pixel_view(right)
        right.columnconfigure(0, weight=2)
        right.columnconfigure(1, weight=1)
        right.rowconfigure(0, weight=1)
        ctl = ttk.Frame(right, padding=(6, 0))
        ctl.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.b_play = ttk.Button(ctl, text="▶ Play", width=9, command=self.toggle_play)
        self.b_play.pack(side="left")
        self.slider = tk.Scale(ctl, from_=0, to=1, orient="horizontal", resolution=1,
                               label="frame", command=self._on_slider)
        self.slider.pack(side="left", fill="x", expand=True, padx=6)
        ttk.Label(ctl, text="stretch").pack(side="left")
        self.stretch = tk.StringVar(value="sqrt")
        cb = ttk.Combobox(ctl, textvariable=self.stretch, values=STRETCH,
                          state="readonly", width=7)
        cb.pack(side="left", padx=(4, 8))
        cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())
        self.show_truth = tk.BooleanVar(value=True)
        ttk.Checkbutton(ctl, text="true track", variable=self.show_truth,
                        command=self._redraw).pack(side="left")
        self.status = ttk.Label(root, anchor="w", padding=(8, 2))
        self.status.grid(row=1, column=0, columnspan=2, sticky="ew")
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        root.protocol("WM_DELETE_WINDOW", self.close)

        params, out = (cv.load_config(config) if config
                       else (cv.VideoParams(), cv.OutputSpec()))
        self._set_values(params, out)
        self.root.after(40, self._poll)
        self._edited()

    # ---- parameter widgets ----------------------------------------------------
    def _build_tab(self, nb, title, specs):
        frm = ttk.Frame(nb, padding=6)
        nb.add(frm, text=title)
        for i, (key, label, kind) in enumerate(specs):
            ttk.Label(frm, text=label).grid(row=i, column=0, sticky="w", pady=1)
            if kind == "b":
                var = tk.BooleanVar()
                w = ttk.Checkbutton(frm, variable=var)
            elif kind == "c":
                var = tk.StringVar()
                w = ttk.Combobox(frm, textvariable=var, values=CHOICES[key][1],
                                 state="readonly", width=26)
            else:
                var = tk.StringVar()
                w = ttk.Entry(frm, textvariable=var, justify="right",
                              width=26 if kind == "s" else 14)
            w.grid(row=i, column=1, sticky="ew" if kind != "b" else "w",
                   padx=(6, 0), pady=1)
            if key == "folder":
                ttk.Button(frm, text="…", width=3, command=self._browse).grid(
                    row=i, column=2, padx=(2, 0))
            var.trace_add("write", self._on_change)
            self.vars[key] = (var, kind)
        frm.columnconfigure(1, weight=1)

    def _set_values(self, params, output):
        self._loading = True
        try:
            lower, upper = cv.spectrum.edges_um(params.lambda_nm, params.bandwidth)
            edges = {"lower_um": lower, "upper_um": upper}
            for key, (var, kind) in self.vars.items():
                if key in edges:
                    var.set(str(edges[key]))
                    continue
                v = getattr(params if key in SIM_KEYS else output, key)
                if kind == "b":
                    var.set(bool(v))
                elif kind == "c":
                    values, labels = CHOICES[key]
                    var.set(labels[values.index(v)])
                elif kind == "t":
                    var.set(", ".join(f"{a:g}" for a in v))
                elif kind == "f":
                    var.set(f"{v:g}")
                else:
                    var.set(str(v))
        finally:
            self._loading = False

    def _read(self):
        sim, out = {}, {}
        sim["lambda_nm"], sim["bandwidth"] = cv.spectrum.from_um(
            self.vars["lower_um"][0].get(), self.vars["upper_um"][0].get())
        for key, (var, kind) in self.vars.items():
            if key in ("lower_um", "upper_um"):
                continue
            v = var.get()
            if kind == "c":
                values, labels = CHOICES[key]
                v = values[labels.index(v)]
            (sim if key in SIM_KEYS else out)[key] = v
        params = cv.VideoParams.from_dict(sim)
        params.validate()
        return params, cv.OutputSpec.from_dict(out)

    def _browse(self):
        path = filedialog.askdirectory(initialdir=cv.DEFAULT_OUTPUT.parent)
        if path:
            self.vars["folder"][0].set(path)

    def _on_change(self, *_):
        if self._loading:
            return
        if self._edit_job:
            self.root.after_cancel(self._edit_job)
        self._edit_job = self.root.after(350, self._edited)

    def _edited(self):
        self._edit_job = None
        try:
            params, out = self._read()
        except ValueError as exc:
            self.params = None
            self.err.config(text=str(exc))
            self._readouts()
            return
        self.err.config(text="")
        same_sim = self.params is not None and params.to_dict() == self.params.to_dict()
        self.params, self.output = params, out
        if same_sim:                      # only output settings changed
            return
        self._gen += 1
        self.slider.config(to=max(params.frames - 1, 0))
        self._readouts()
        self._goto(int(self.slider.get()))

    # ---- readouts ---------------------------------------------------------------
    def _readouts(self):
        p = self.params
        if p is None:
            self.readout.config(text="(fix the parameter above)")
            return
        lod = p.lf_um / p.pixel_um
        s = p.speed_px_s
        nyq = 0.5 * p.lf_min_um / p.pixel_um
        bytes_px = 1 if p.bit_depth <= 8 else 2
        lines = [
            f"λ/D      {lod:.2f} px = {lod * p.pixel_arcsec:.3f}″",
            f"pixel    {p.pixel_arcsec:.3f}″  ({nyq:.2f} × Nyquist at λmin)",
            f"donut ⌀  {8 * abs(p.defocus_waves) * lod:.0f} px",
            f"speed    {s:.1f} px/s = {s * p.pixel_um:.0f} µm/s",
            f"         = {s * p.pixel_arcsec:.2f}″/s = {s * p.pixel_arcsec / 3600:.4f}°/s",
            f"step     {s / p.fps:.2f} px/frame",
            f"blur     {p.blur_px:.2f} px per exposure (duty {p.exposure_s * p.fps:.0%})",
            f"length   {p.frames / p.fps:.2f} s, travel {s * p.frames / p.fps:.0f} px",
            f"raw data {p.frames * p.width * p.height * bytes_px / 1e6:.0f} MB per file"]
        sim = self.last_sim
        if sim is not None and sim.p.to_dict() == p.to_dict():
            d = sim.derived()
            lines += [
                f"signal   {d['signal_e_per_frame']:.3g} e⁻/frame",
                f"         = {d['signal_e_per_s']:.3g} e⁻/s",
                f"peak     {d['static_peak_percent_fs']:.0f}% FS static → "
                f"{d['moving_peak_percent_fs']:.0f}% moving",
                f"backgr.  {d['background_e_per_frame']:.1f} e⁻/px/frame",
                f"grid     {d['fine_grid']}² (pixel/{d['oversampling']})"]
        else:
            lines += ["signal   …", "", "peak     …", "backgr.  …", "grid     …"]
        self.readout.config(text="\n".join(lines))

    # ---- rendering / background work ----------------------------------------
    def _goto(self, k):
        """Show frame k unless it is already requested for these parameters.
        (Called directly, not only from the slider callback: Tk runs a
        scale's command only when the widget redraws.)"""
        if self._want == (k, self._gen):
            return
        self._want = (k, self._gen)
        self.request_frame(k)

    def request_frame(self, k):
        """Show frame k: read it from the last run when that run used the
        current parameters, otherwise render it in the background."""
        run = self.run
        if (run is not None and run["frames"] is not None and self.params is not None
                and run["params"] == self.params.to_dict() and k < len(run["frames"])):
            self._display(k, np.asarray(run["frames"][k]), run["sim"])
            return
        self.pending = k
        self._kick()

    def _kick(self):
        if self.worker is not None:
            return
        if self.run_request is not None:
            params, out = self.run_request
            self.run_request = None
            self.cancel.clear()
            self._busy(True)
            self.worker = threading.Thread(target=self._run_worker,
                                           args=(params, out), daemon=True)
            self.worker.start()
        elif self.pending is not None and self.params is not None:
            k, self.pending = self.pending, None
            self.worker = threading.Thread(target=self._frame_worker,
                                           args=(self.params, k), daemon=True)
            self.worker.start()

    def _frame_worker(self, params, k):
        try:
            t0 = time.perf_counter()
            sim = cv.VideoSimulator(params)
            k = min(k, params.frames - 1)
            dn, _ = sim.render_frame(k)
            self.q.put(("frame", k, dn, sim, time.perf_counter() - t0))
        except Exception as exc:          # report any failure in the GUI
            self.q.put(("error", f"{type(exc).__name__}: {exc}"))

    def _run_worker(self, params, out):
        try:
            sim = cv.VideoSimulator(params)
            res = sim.run(out, progress=lambda d, t: self.q.put(("progress", d, t)),
                          cancel=self.cancel, log=lambda m: self.q.put(("log", m)))
            self.q.put(("run_done", res, sim))
        except Exception as exc:
            self.q.put(("run_error", f"{type(exc).__name__}: {exc}"))

    def _poll(self):
        try:
            while True:
                msg = self.q.get_nowait()
                kind = msg[0]
                if kind == "frame":
                    _, k, dn, sim, dt = msg
                    self.worker = None
                    self.last_sim = sim
                    self._display(k, dn, sim)
                    self._readouts()
                    self.status.config(text=f"frame {k} rendered in {dt * 1e3:.0f} ms "
                                            f"(fine grid {sim.model.N}²)")
                elif kind == "progress":
                    self.progress["value"] = msg[1] / msg[2]
                    self.status.config(text=f"rendering frame {msg[1]} / {msg[2]} …")
                elif kind == "log":
                    self.status.config(text=msg[1])
                elif kind == "run_done":
                    self._run_finished(*msg[1:])
                elif kind in ("error", "run_error"):
                    self.worker = None
                    if kind == "run_error":
                        self._busy(False)
                    self.err.config(text=msg[1])
        except queue.Empty:
            pass
        self._kick()
        self.root.after(40, self._poll)

    # ---- run / playback ---------------------------------------------------------
    def start_run(self):
        if self.params is None:
            return
        self.playing = False
        self.b_play.config(text="▶ Play")
        self._release_run()
        self.run_request = (self.params, self.output)
        self.status.config(text="starting the run …")
        self._kick()

    def _busy(self, busy):
        self.b_run.config(state="disabled" if busy else "normal")
        self.b_cancel.config(state="normal" if busy else "disabled")
        if busy:
            self.progress["value"] = 0.0

    def _run_finished(self, res, sim):
        self.worker = None
        self._busy(False)
        folder = Path(res["folder"])
        npy = folder / "frames.npy"
        frames = np.load(npy, mmap_mode="r") if npy.exists() else None
        self.run = dict(frames=frames, params=sim.p.to_dict(), sim=sim, folder=folder)
        self.last_sim = sim
        self._readouts()
        what = "cancelled after" if res["cancelled"] else "saved"
        self.status.config(text=f"{what} {res['frames_written']} frames in "
                                f"{res['render_seconds']:.1f} s → {folder}")
        k = int(self.slider.get())
        self._want = (k, self._gen)
        self.request_frame(k)            # now served from the saved frames

    def _release_run(self):
        if self.run is not None:
            self.run["frames"] = None     # close the memmap (file handle)
        self.run = None

    def toggle_play(self):
        self.playing = not self.playing
        self.b_play.config(text="❚❚ Pause" if self.playing else "▶ Play")
        if self.playing:
            self._play_tick()

    def _play_tick(self):
        """Advance one frame once the current one is on screen, then leave
        the event loop some slack (a timer that is always due would starve
        Tk's idle drawing)."""
        if not self.playing or self.params is None:
            return
        t0 = time.perf_counter()
        k = int(self.slider.get())
        if self.cur is not None and self.cur[0] == k:
            k = (k + 1) % self.params.frames
            self.slider.set(k)
            self._goto(k)
        spent = (time.perf_counter() - t0) * 1e3
        interval = 1e3 / min(self.params.fps, 25.0)
        self.root.after(max(10, int(interval - spent)), self._play_tick)

    def _on_slider(self, value):
        self._goto(int(float(value)))

    # ---- drawing ----------------------------------------------------------------
    def _init_artists(self):
        """Persistent artists: frames update with set_data (fast playback)."""
        blank = np.zeros((2, 2))
        self.im = self.ax.imshow(blank, cmap="gray", vmin=0.0, vmax=1.0,
                                 interpolation="nearest")
        (self.track_line,) = self.ax.plot([], [], "-", color="tab:red", lw=0.7, alpha=0.8)
        (self.streak,) = self.ax.plot([], [], "-", color="cyan", lw=1.5)
        self.imz = self.axz.imshow(blank, cmap="inferno", vmin=0.0, vmax=1.0,
                                   interpolation="nearest")
        (self.streak_z,) = self.axz.plot([], [], "-", color="cyan", lw=1.2)
        (self.mark_z,) = self.axz.plot([], [], "+", color="cyan", ms=10, mew=1.2)
        self.msg_z = self.axz.text(0.5, 0.5, "", ha="center", va="center",
                                   transform=self.axz.transAxes)
        for a in (self.ax, self.axz):
            a.tick_params(labelsize=7)

    def _display(self, k, dn, sim):
        self.cur = (k, dn, sim)
        self._redraw()

    def _stretch(self, dn, lo, hi):
        v = np.clip((dn.astype(np.float32) - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
        mode = self.stretch.get()
        if mode == "sqrt":
            return np.sqrt(v)
        if mode == "log":
            return np.clip(1.0 + np.log10(np.maximum(v, 1e-4)) / 4.0, 0.0, 1.0)
        return v

    def _track_of(self, sim):
        if self._track[0] is not sim:
            pts = [sim.truth_geometry(i) for i in range(sim.p.frames)]
            self._track = (sim, np.array([[g["x_apparent_px"], g["y_apparent_px"]]
                                          for g in pts]))
        return self._track[1]

    def _redraw(self):
        if self.cur is None:
            return
        k, dn, sim = self.cur
        p = sim.p
        img = self._stretch(dn, *sim.display_range())
        geo = sim.truth_geometry(k)
        xa, ya = geo["x_apparent_px"], geo["y_apparent_px"]
        truth = self.show_truth.get()
        streak = ([geo["x_start_px"], geo["x_end_px"]], [geo["y_start_px"], geo["y_end_px"]])
        self.im.set_data(img)
        self.im.set_extent((-0.5, p.width - 0.5, p.height - 0.5, -0.5))
        tr = self._track_of(sim) if truth else np.empty((0, 2))
        self.track_line.set_data(tr[:, 0], tr[:, 1])
        self.streak.set_data(*(streak if truth else ([], [])))
        self.ax.set_xlim(-0.5, p.width - 0.5)
        self.ax.set_ylim(p.height - 0.5, -0.5)
        self.ax.set_title(f"frame {k} / {p.frames - 1}    t = {geo['t_mid_s']:.4f} s    "
                          f"{p.width}×{p.height} px", fontsize=9)

        half = int(0.5 * spot_extent_px(p))
        x0, x1 = max(int(xa) - half, 0), min(int(xa) + half + 1, p.width)
        y0, y1 = max(int(ya) - half, 0), min(int(ya) + half + 1, p.height)
        inside = x0 < x1 and y0 < y1
        self.imz.set_visible(inside)
        self.msg_z.set_text("" if inside else "spot outside the frame")
        self.streak_z.set_data(*(streak if truth and inside else ([], [])))
        self.mark_z.set_data(*(([xa], [ya]) if truth and inside else ([], [])))
        if inside:
            self.imz.set_data(img[y0:y1, x0:x1])
            self.imz.set_extent((x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5))
            self.axz.set_xlim(x0 - 0.5, x1 - 0.5)
            self.axz.set_ylim(y1 - 0.5, y0 - 0.5)
            self.axz.set_title(f"spot zoom    peak {int(dn[y0:y1, x0:x1].max())} DN "
                               f"of {p.dn_max}", fontsize=9)
        else:
            self.axz.set_title("spot zoom", fontsize=9)
        self._view = (img, dn, sim, geo)
        self._draw_pixels()
        if self.playing:
            self.canvas.draw()            # synchronous: one frame per tick
        else:
            self.canvas.draw_idle()

    # ---- exact pixel view -------------------------------------------------------
    def _build_pixel_view(self, parent):
        """Tk canvas right of the figure: images on a Tk canvas are shown
        pixel for pixel, so magnified camera pixels stay exact squares."""
        box = ttk.LabelFrame(parent, text="Pixel view", padding=(4, 2))
        box.grid(row=0, column=1, sticky="nsew", padx=(0, 6), pady=(4, 2))
        self.pv_box = box
        side = int(round(320 * self.root.winfo_fpixels("1i") / 96.0))  # DPI-scaled
        self.pv = c = tk.Canvas(box, width=side, height=side, background="#1e1e1e",
                                highlightthickness=0, cursor="crosshair")
        c.grid(row=0, column=0, columnspan=4, sticky="nsew")
        ttk.Label(box, text="zoom").grid(row=1, column=0, sticky="w", pady=(4, 0))
        self.pv_zoom = tk.StringVar(value=PIXEL_ZOOMS[0])
        cb = ttk.Combobox(box, textvariable=self.pv_zoom, values=PIXEL_ZOOMS,
                          state="readonly", width=8)
        cb.grid(row=1, column=1, sticky="w", padx=4, pady=(4, 0))
        cb.bind("<<ComboboxSelected>>", lambda _e: self._draw_pixels())
        self.pv_grid = tk.BooleanVar(value=True)
        ttk.Checkbutton(box, text=f"pixel grid (from {GRID_MIN_ZOOM}×)",
                        variable=self.pv_grid, command=self._draw_pixels).grid(
            row=1, column=2, sticky="w", pady=(4, 0))
        # width=1 + sticky "ew": long texts are clipped instead of widening the panel
        self.pv_info = ttk.Label(box, width=1, anchor="w")
        self.pv_info.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(2, 0))
        ttk.Label(box, text=PV_HINT, width=1, anchor="w", foreground="gray40").grid(
            row=3, column=0, columnspan=4, sticky="ew")
        box.columnconfigure(3, weight=1)
        box.rowconfigure(0, weight=1)
        c.bind("<Configure>", lambda _e: self._draw_pixels())
        c.bind("<Motion>", self._pv_hover)
        c.bind("<Leave>", self._pv_leave)
        c.bind("<MouseWheel>", self._pv_wheel)
        c.bind("<Button-4>", lambda e: self._pv_wheel(e, +1))     # X11 wheel
        c.bind("<Button-5>", lambda e: self._pv_wheel(e, -1))
        c.bind("<ButtonPress-1>", self._pv_press)
        c.bind("<B1-Motion>", self._pv_drag_to)
        c.bind("<ButtonRelease-1>", self._pv_release)
        c.bind("<Double-Button-1>", self._pv_recentre)

    def _pv_zoom_factor(self, p, cw, ch):
        """Screen pixels per camera pixel (an integer, so blocks are exact)."""
        label = self.pv_zoom.get()
        if label in PIXEL_ZOOMS[1:]:
            return ZOOM_STEPS[PIXEL_ZOOMS.index(label) - 1]
        fit = min(cw, ch) // max(spot_extent_px(p), 1.0)   # largest zoom that fits
        return int(min(max(fit, 1), ZOOM_STEPS[-1]))

    def _draw_pixels(self):
        """Every camera pixel becomes a uniform k x k block of screen pixels
        (integer magnification, never resampled), centred on the spot plus
        the panning offset; the same display stretch as the other views."""
        c = self.pv
        c.delete("all")
        self._pv_map = None
        if self._view is None:
            return
        img, dn, sim, geo = self._view
        p = sim.p
        cw = c.winfo_width() if c.winfo_width() > 1 else int(c["width"])
        ch = c.winfo_height() if c.winfo_height() > 1 else int(c["height"])
        k = self._pv_zoom_factor(p, cw, ch)
        nx, ny = -(-cw // k), -(-ch // k)            # blocks covering the canvas
        x0 = int(round(geo["x_apparent_px"] + self._pv_offset[0])) - nx // 2
        y0 = int(round(geo["y_apparent_px"] + self._pv_offset[1])) - ny // 2
        rgb = np.empty((ny, nx, 3), np.uint8)
        rgb[:] = OUTSIDE_RGB
        xs0, xs1 = max(x0, 0), min(x0 + nx, p.width)
        ys0, ys1 = max(y0, 0), min(y0 + ny, p.height)
        inside = xs0 < xs1 and ys0 < ys1
        if inside:
            rgb[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0] = (
                255.0 * img[ys0:ys1, xs0:xs1] + 0.5).astype(np.uint8)[..., None]
        big = np.ascontiguousarray(rgb.repeat(k, axis=0).repeat(k, axis=1)[:ch, :cw])
        if self.pv_grid.get() and k >= GRID_MIN_ZOOM:
            for lines in (big[::k], big[:, ::k]):  # first row / column of each block
                lines[...] = ((lines + GRID_RGB) // 2).astype(np.uint8)
        self._pv_rgb = big
        self._pv_photo = ImageTk.PhotoImage(Image.fromarray(big), master=c)
        c.create_image(0, 0, anchor="nw", image=self._pv_photo)
        title = ("Pixel view  ·  1:1 (actual pixels)" if k == 1
                 else f"Pixel view  ·  1 px = {k}×{k} screen px")
        if str(self.pv_box["text"]) != title:
            self.pv_box.config(text=title)

        def at(x, y):                                # sensor px -> canvas
            return (x - x0 + 0.5) * k, (y - y0 + 0.5) * k

        if not inside:
            c.create_text(cw / 2, ch / 2, fill="white", justify="center",
                          text="spot outside the sensor" if self._pv_offset == (0, 0)
                          else "outside the sensor\n(double-click to re-centre)")
        elif self.show_truth.get():
            c.create_line(*at(geo["x_start_px"], geo["y_start_px"]),
                          *at(geo["x_end_px"], geo["y_end_px"]), fill="cyan", width=2)
            xm, ym = at(geo["x_apparent_px"], geo["y_apparent_px"])
            c.create_line(xm - 7, ym, xm + 8, ym, fill="cyan")
            c.create_line(xm, ym - 7, xm, ym + 8, fill="cyan")
        # scale bar: n camera pixels with their size on the sensor and the sky
        n = max((s for s in (1, 2, 5, 10, 20, 50, 100, 200, 500) if s * k <= 0.4 * cw),
                default=1)
        y_bar = ch - 10
        bar = c.create_line(10, y_bar, 10 + n * k, y_bar, fill="white", width=3)
        label = c.create_text(10, y_bar - 5, anchor="sw", fill="white",
                              text=f"{n} px = {n * p.pixel_um:.4g} µm = "
                                   f"{n * p.pixel_arcsec:.3g}″")
        bx0, by0, bx1, _ = c.bbox(label)
        back = c.create_rectangle(bx0 - 4, by0 - 2, max(bx1, 10 + n * k) + 4, ch - 3,
                                  fill="black", outline="")
        c.tag_lower(back, bar)
        self._pv_map = (x0, y0, k, dn, p)
        self._pv_readout()

    def _pv_readout(self):
        """Value of the camera pixel under the mouse."""
        m, pos = self._pv_map, self._pv_mouse
        if m is None or pos is None:
            self.pv_info.config(text="hover over a pixel to read its value")
            return
        x0, y0, k, dn, p = m
        ix, iy = x0 + int(pos[0]) // k, y0 + int(pos[1]) // k
        if 0 <= ix < p.width and 0 <= iy < p.height:
            v = int(dn[iy, ix])
            text = (f"pixel ({ix}, {iy}):  {v} DN  ≈ "
                    f"{(v - p.offset_dn) * p.gain_e_per_dn:.0f} e⁻ above black")
        else:
            text = f"({ix}, {iy}): outside the sensor"
        self.pv_info.config(text=text)

    def _pv_hover(self, event):
        self._pv_mouse = (event.x, event.y)
        self._pv_readout()

    def _pv_leave(self, _event):
        self._pv_mouse = None
        self._pv_readout()

    def _pv_wheel(self, event, step=0):
        if self._pv_map is None:
            return
        up = step > 0 if step else event.delta > 0
        k = self._pv_map[2]
        steps = [z for z in ZOOM_STEPS if (z > k if up else z < k)]
        if steps:
            k = steps[0] if up else steps[-1]
            self.pv_zoom.set(PIXEL_ZOOMS[ZOOM_STEPS.index(k) + 1])
            self._draw_pixels()

    def _pv_press(self, event):
        self._pv_drag = (event.x, event.y, self._pv_offset)

    def _pv_drag_to(self, event):
        self._pv_mouse = (event.x, event.y)
        if self._pv_drag is None or self._pv_map is None:
            return
        x, y, (ox, oy) = self._pv_drag
        k = self._pv_map[2]
        offset = (ox - round((event.x - x) / k), oy - round((event.y - y) / k))
        if offset != self._pv_offset:             # whole pixels: blocks stay aligned
            self._pv_offset = offset
            self._draw_pixels()
        else:
            self._pv_readout()

    def _pv_release(self, _event):
        self._pv_drag = None

    def _pv_recentre(self, _event=None):
        self._pv_offset = (0, 0)
        self._draw_pixels()

    # ---- configuration files --------------------------------------------------
    def save_cfg(self):
        try:
            params, out = self._read()
        except ValueError as exc:
            self.err.config(text=str(exc))
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".json", filetypes=[("JSON", "*.json")],
            initialdir=cv.ROOT, title="Save configuration")
        if path:
            cv.save_config(path, params, out)
            self.status.config(text=f"saved {path}")

    def load_cfg(self):
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json")],
                                          initialdir=cv.ROOT, title="Load configuration")
        if not path:
            return
        try:
            params, out = cv.load_config(path)
            self._set_values(params, out)
        except (ValueError, OSError) as exc:
            self.err.config(text=f"{Path(path).name}: {exc}")
            return
        self._edited()
        self.status.config(text=f"loaded {path}")

    def close(self):
        self.cancel.set()
        self.playing = False
        self._release_run()
        self.root.destroy()


def main():
    try:  # crisp text on Windows high-DPI displays
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass
    root = tk.Tk()
    App(root, sys.argv[1] if len(sys.argv) > 1 else None)
    root.update_idletasks()
    if (root.winfo_reqwidth() > 0.95 * root.winfo_screenwidth()
            or root.winfo_reqheight() > 0.9 * root.winfo_screenheight()):
        try:
            root.state("zoomed")          # larger than the screen: start maximized
        except tk.TclError:
            pass
    root.mainloop()


if __name__ == "__main__":
    main()
