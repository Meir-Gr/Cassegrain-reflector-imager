"""
Cassegrain video simulator: a camera recording of a spot in linear motion.

Every frame is the spot image of cassegrain_gui.py (same scalar-diffraction
model: annular pupil with a 3-vane spider, defocus W020, flat-spectrum
bandwidth, uniform-disk source, Gaussian in-frame jitter) as a real sensor
records it:

  * constant-velocity linear motion: during one exposure the spot travels
    L = speed x exposure (an exact box blur, MTF sinc(L q_par) along the
    motion); between frames it advances speed / fps;
  * pixel integration (box MTF of the pixel pitch) and sampling on the
    sensor grid, with the aliasing of undersampled sensors included exactly
    (spectrum folding); sub-pixel spot positions by Fourier shift;
  * photon shot noise, dark current, uniform background, full-well
    saturation, read noise, gain / offset / bit-depth quantization;
  * optional rolling shutter (row r exposed r x line time later) and random
    frame-to-frame pointing jitter.

Conventions: x = column, y = row in pixels (origin at the centre of the
top-left pixel, y down). (start_x, start_y) is the spot centre at the middle
of the first frame's exposure (rolling shutter: of row 0). The motion
direction is counter-clockwise from +x as seen on screen (0 = right,
90 = up). Frame k is exposed from k/fps to k/fps + exposure. The PSF is the
same everywhere on the sensor (on-axis model).

Intensity: either the photoelectron rate of the whole spot [e-/s] (signal
per frame = rate x exposure, so motion blur lowers the peak), or the peak a
static spot would reach [% of full scale], converted to the same rate.

Outputs (output/<name>/, chosen in OutputSpec): frames.npy, frames.tif,
frames.mat (MATLAB: frames H x W x N plus truth / params), raw frames
(headerless little-endian binary, one file per frame in raw/ or all frames
in frames.raw, layout in raw_format.txt), preview.mp4 / .gif, truth.csv
(true position per frame) and params.json (all settings).

Run:  python cassegrain_video.py [run.json] [--set key=value ...]
      python cassegrain_video.py --write-config run.json   (defaults, to edit)
      python cassegrain_video.py --selftest
      python cassegrain_video_gui.py                       (interactive)
"""

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import asdict, dataclass, fields
from pathlib import Path

import numpy as np
from scipy import fft as sfft
from scipy.special import j1

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = ROOT / "output"
RAD_TO_ARCSEC = 180.0 / math.pi * 3600.0
SPEED_UNITS = ("px/s", "um/s", "arcsec/s", "deg/s")
INTENSITY_MODES = ("static_peak", "e_per_s")
PREVIEW_KINDS = ("mp4", "gif", "none")
RAW_KINDS = ("none", "per_frame", "single")
MAT_PART_BYTES = 2 ** 30     # frames per .mat file stay below this (MAT v5: < 2 GB)
MAX_GRID = 8192              # largest fine-grid FFT (complex64 N^2 = 512 MB)
MIN_PUPIL_RADIUS = 96        # pupil samples: keeps grid-induced centroid bias < 0.001 px
GIF_MAX_FRAMES = 300
PREVIEW_MAX_FPS = 60.0


# ---------------------------------------------------------------- parameters
@dataclass
class VideoParams:
    """All simulation inputs (units in the comments)."""
    # optics (same meaning as in cassegrain_gui.py)
    d_mm: float = 200.0                  # primary diameter
    obstruction: float = 0.35            # secondary / primary diameter
    vane_mm: float = 4.0                 # spider vane width
    vane_angles_deg: tuple = (90.0, 210.0, 330.0)
    lambda_nm: float = 550.0             # central wavelength
    fnum: float = 10.0                   # focal ratio
    defocus_waves: float = 8.0           # W020 at the central wavelength
    bandwidth: float = 0.0               # d(lambda)/lambda, flat spectrum
    source_lod: float = 0.0              # source diameter [lambda/D]
    jitter_lod: float = 0.0              # in-frame jitter RMS [lambda/D]
    # sensor
    width: int = 640                     # W [px]
    height: int = 512                    # H [px]
    pixel_um: float = 3.45               # pixel pitch
    bit_depth: int = 12
    gain_e_per_dn: float = 2.5
    offset_dn: float = 64.0              # black level
    full_well_e: float = 10000.0
    read_noise_e: float = 2.5            # RMS
    dark_e_per_s: float = 5.0            # per pixel
    background_e_per_s: float = 0.0      # per pixel (sky / scene)
    shot_noise: bool = True
    rolling_shutter: bool = False
    line_time_us: float = 10.0           # row-to-row delay (rolling shutter)
    # sequence and motion
    frames: int = 60
    fps: float = 30.0
    exposure_ms: float = 20.0
    start_x: float = 120.0               # spot centre in frame 0 [px]
    start_y: float = 300.0
    speed: float = 200.0                 # in speed_unit
    speed_unit: str = "px/s"             # px/s, um/s, arcsec/s, deg/s
    direction_deg: float = 20.0          # counter-clockwise from +x on screen
    frame_jitter_px: float = 0.0         # frame-to-frame pointing RMS per axis
    # intensity
    intensity_mode: str = "static_peak"  # "static_peak" [% full scale] or "e_per_s"
    intensity: float = 50.0
    seed: int = 1

    def __post_init__(self):
        self.vane_angles_deg = tuple(float(a) for a in self.vane_angles_deg)

    # ---- derived scales ----------------------------------------------------
    @property
    def lf_um(self):
        """lambda F# [um]: the focal-plane size of lambda/D."""
        return self.lambda_nm * 1e-3 * self.fnum

    @property
    def focal_mm(self):
        return self.d_mm * self.fnum

    @property
    def pixel_arcsec(self):
        return self.pixel_um * 1e-3 / self.focal_mm * RAD_TO_ARCSEC

    @property
    def speed_px_s(self):
        to_px = {"px/s": 1.0, "um/s": 1.0 / self.pixel_um,
                 "arcsec/s": 1.0 / self.pixel_arcsec,
                 "deg/s": 3600.0 / self.pixel_arcsec}
        return self.speed * to_px[self.speed_unit]

    @property
    def velocity_px_s(self):
        """(vx, vy) [px/s] with y pointing down the rows."""
        a = math.radians(self.direction_deg)
        return self.speed_px_s * math.cos(a), -self.speed_px_s * math.sin(a)

    @property
    def exposure_s(self):
        return self.exposure_ms * 1e-3

    @property
    def blur_px(self):
        return self.speed_px_s * self.exposure_s

    @property
    def dn_max(self):
        return 2 ** int(self.bit_depth) - 1

    @property
    def full_scale_e(self):
        """Signal that saturates a pixel: the full well or the ADC range."""
        return min(self.full_well_e,
                   (self.dn_max - self.offset_dn) * self.gain_e_per_dn)

    # ---- checks and (de)serialization --------------------------------------
    def validate(self):
        def need(ok, msg):
            if not ok:
                raise ValueError(msg)
        need(self.d_mm > 0 and self.lambda_nm > 0 and self.fnum > 0,
             "D, wavelength and F# must be positive")
        need(0.0 <= self.obstruction <= 0.95, "obstruction must be in 0 .. 0.95")
        need(self.vane_mm >= 0, "vane width cannot be negative")
        need(len(self.vane_angles_deg) == 3, "three vane angles are required")
        need(0.0 <= self.bandwidth <= 0.5, "bandwidth must be in 0 .. 0.5")
        need(self.source_lod >= 0 and self.jitter_lod >= 0,
             "source size and jitter cannot be negative")
        need(8 <= self.width <= 16384 and 8 <= self.height <= 16384,
             "image width and height must be 8 .. 16384 px")
        need(self.pixel_um > 0, "pixel size must be positive")
        need(1 <= self.bit_depth <= 16, "bit depth must be 1 .. 16")
        need(self.gain_e_per_dn > 0 and self.full_well_e > 0,
             "gain and full well must be positive")
        need(0 <= self.offset_dn < self.dn_max, "offset must lie inside the ADC range")
        need(min(self.read_noise_e, self.dark_e_per_s, self.background_e_per_s) >= 0,
             "noise and background levels cannot be negative")
        need(self.frames >= 1, "the sequence needs at least one frame")
        need(self.fps > 0 and self.exposure_ms > 0,
             "frame rate and exposure must be positive")
        need(self.exposure_s <= (1.0 + 1e-9) / self.fps,
             f"exposure {self.exposure_ms:g} ms is longer than the frame period "
             f"{1e3 / self.fps:.4g} ms")
        need(self.speed >= 0, "speed cannot be negative (set the direction instead)")
        need(self.speed_unit in SPEED_UNITS, f"speed unit must be one of {SPEED_UNITS}")
        need(self.frame_jitter_px >= 0, "frame jitter cannot be negative")
        need(self.intensity_mode in INTENSITY_MODES,
             f"intensity mode must be one of {INTENSITY_MODES}")
        need(self.intensity >= 0, "intensity cannot be negative")
        if self.rolling_shutter:
            need(self.line_time_us > 0, "line time must be positive")
            readout = self.height * self.line_time_us * 1e-6
            need(readout <= (1.0 + 1e-9) / self.fps,
                 f"rolling-shutter readout ({readout * 1e3:.4g} ms for {self.height} "
                 f"rows) is longer than the frame period")
            need(abs(self.velocity_px_s[1] * self.line_time_us * 1e-6) < 0.5,
                 "the spot moves vertically faster than half a row per line time")

    def to_dict(self):
        d = asdict(self)
        d["vane_angles_deg"] = list(self.vane_angles_deg)
        return d

    @classmethod
    def from_dict(cls, data):
        return cls(**_coerce_all(cls, data))


@dataclass
class OutputSpec:
    """What a run writes, and where: <folder>/<name>/."""
    name: str = ""          # "" -> video_YYYYmmdd_HHMMSS
    folder: str = ""        # "" -> <project>/output
    npy: bool = True        # frames.npy: lossless DN, shape (frames, H, W)
    tiff: bool = True       # frames.tif: multi-page 8/16-bit stack
    mat: bool = False       # frames.mat: MATLAB frames (H x W x N) + truth + params
    raw: str = "none"       # raw binary: "per_frame" (raw/*.raw), "single" (frames.raw)
    preview: str = "mp4"    # preview.mp4 / preview.gif / "none"
    truth: bool = True      # truth.csv: true spot position per frame

    def run_dir(self):
        base = Path(self.folder) if self.folder else DEFAULT_OUTPUT
        return base / (self.name or time.strftime("video_%Y%m%d_%H%M%S"))

    @classmethod
    def from_dict(cls, data):
        out = cls(**_coerce_all(cls, data))
        if out.preview not in PREVIEW_KINDS:
            raise ValueError(f"preview must be one of {PREVIEW_KINDS}")
        if out.raw not in RAW_KINDS:
            raise ValueError(f"raw must be one of {RAW_KINDS}")
        return out


def _coerce(f, value):
    """Convert a JSON / command-line value to the type of field f."""
    default = f.default
    if isinstance(default, bool):
        if isinstance(value, str):
            s = value.strip().lower()
            if s in ("1", "true", "yes", "on"):
                return True
            if s in ("0", "false", "no", "off"):
                return False
            raise ValueError(f"{f.name}: expected true/false, got {value!r}")
        return bool(value)
    if isinstance(default, int):
        x = float(value)
        if not x.is_integer():
            raise ValueError(f"{f.name}: expected an integer, got {value!r}")
        return int(x)
    if isinstance(default, float):
        return float(value)
    if isinstance(default, tuple):
        if isinstance(value, str):
            value = value.split(",")
        return tuple(float(v) for v in value)
    return str(value)


def _coerce_all(cls, data):
    known = {f.name: f for f in fields(cls)}
    unknown = sorted(set(data) - set(known))
    if unknown:
        raise ValueError(f"unknown {cls.__name__} key(s) {unknown}; "
                         f"valid keys: {sorted(known)}")
    return {k: _coerce(known[k], v) for k, v in data.items()}


def save_config(path, params, output=None):
    data = {"simulation": params.to_dict(), "output": asdict(output or OutputSpec())}
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")


def load_config(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    extra = sorted(set(data) - {"simulation", "output"})
    if extra:
        raise ValueError(f"{path}: unknown section(s) {extra}; "
                         "expected 'simulation' and 'output'")
    return (VideoParams.from_dict(data.get("simulation", {})),
            OutputSpec.from_dict(data.get("output", {})))


def apply_overrides(params, output, pairs):
    """Apply 'key=value' strings to a VideoParams / OutputSpec pair."""
    sim_f = {f.name: f for f in fields(VideoParams)}
    out_f = {f.name: f for f in fields(OutputSpec)}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        key = key.strip()
        if not sep:
            raise ValueError(f"--set expects key=value, got {pair!r}")
        if key in sim_f:
            setattr(params, key, _coerce(sim_f[key], value))
        elif key in out_f:
            setattr(output, key, _coerce(out_f[key], value))
        else:
            raise ValueError(f"unknown key {key!r}")


# ----------------------------------------------------------------- spot model
def n_lambda(w020, bandwidth):
    """Spectral samples for a flat band (same rule as cassegrain_gui.py)."""
    if bandwidth <= 0:
        return 1
    return int(np.clip(3 + 20.0 * bandwidth + 0.7 * bandwidth * abs(w020), 3, 15))


def _even_fast_len(n):
    k = sfft.next_fast_len(max(2, int(math.ceil(n))))
    while k % 2:
        k = sfft.next_fast_len(k + 1)
    return k


_PSF_CACHE = {}
_MODEL_CACHE = {}


def _cover(d, a, b):
    """Area of a unit pixel lying on the side u <= d of a straight edge
    (u = distance from the pixel centre along the edge normal, whose
    components have magnitudes a, b): the CDF of U(-a/2, a/2) + U(-b/2, b/2),
    exact at any edge angle."""
    a, b = np.maximum(a, b), np.minimum(a, b)
    h1, h2 = 0.5 * (a + b), 0.5 * (a - b)
    ab2 = np.maximum(2.0 * a * b, 1e-9)
    return np.where(d <= -h1, 0.0,
           np.where(d < -h2, (d + h1) ** 2 / ab2,
           np.where(d <= h2, (d + 0.5 * a) / a,
           np.where(d < h1, 1.0 - (h1 - d) ** 2 / ab2, 1.0))))


def _psf_fine(p, pad0, n):
    """Polychromatic PSF (sum 1, centre n//2) sampled at pad0 samples per
    lambda0 F#. The pupil Y axis is flipped so that row 0 is the top of the
    image: displayed with row 0 at the top, it matches cassegrain_gui.py.
    Pupil samples transmit power equal to their exact open area (anti-
    aliased edges), so thin vanes keep their true area and the spider
    geometry sets the centroid, not the sampling grid."""
    key = (p.obstruction, p.vane_mm / p.d_mm, tuple(p.vane_angles_deg),
           p.defocus_waves, p.bandwidth, round(pad0, 12), n)
    if key in _PSF_CACHE:
        return _PSF_CACHE[key]
    c = (np.arange(n) - n // 2).astype(np.float32)
    X, Y = np.meshgrid(c, -c)
    rr2 = X * X + Y * Y
    r = np.sqrt(rr2)
    rs = np.where(r > 0, r, 1.0)
    na, nb = np.where(r > 0, np.abs(X) / rs, 1.0), np.abs(Y) / rs  # radial normal
    trig = [(np.float32(math.cos(a)), np.float32(math.sin(a)))
            for a in np.deg2rad(p.vane_angles_deg)]
    nl = n_lambda(p.defocus_waves, p.bandwidth)
    acc = np.zeros((n, n), np.float32)
    for k in range(nl):
        rel = (1.0 if nl == 1 else
               1.0 - p.bandwidth / 2.0 + p.bandwidth * (k + 0.5) / nl)
        R = n / (2.0 * pad0 * rel)          # pupil radius [samples], ~ 1/lambda
        T = _cover(R - r, na, nb)
        if p.obstruction > 0:
            T *= 1.0 - _cover(p.obstruction * R - r, na, nb)
        hw = R * p.vane_mm / p.d_mm
        if hw > 0:
            for ca, sa in trig:
                across = Y * ca - X * sa
                strip = (_cover(hw - across, abs(sa), abs(ca))
                         - _cover(-hw - across, abs(sa), abs(ca)))
                T *= 1.0 - strip * (X * ca + Y * sa >= 0)
        T = T.astype(np.float32)
        # amplitude = sqrt(open area): each sample transmits power equal to
        # its open area, so energy and centroid follow the true pupil shape
        T = np.sqrt(T)
        # achromatic mirrors: same physical defocus, W020 ~ 1/lambda
        phase = np.float32(2.0 * np.pi * p.defocus_waves / rel / (R * R)) * rr2
        field = sfft.fft2(np.fft.ifftshift(T * np.exp(1j * phase)), workers=-1)
        psf = np.abs(np.fft.fftshift(field)) ** 2
        acc += psf / np.float32((T * T).sum())   # equal energy per wavelength
    # fade the outermost halo smoothly to zero (circular cosine taper over the
    # outer quarter of the field): otherwise the field edge shows as a faint
    # square step, and the taper also keeps the PSF point-symmetric on the
    # even grid (its unpaired edge row / column falls outside the taper)
    r_edge = n / 2.0 - 1.0
    fade = np.clip((r_edge - r) / (0.25 * r_edge), 0.0, 1.0)
    acc *= 0.5 - 0.5 * np.cos(np.pi * fade)
    acc /= acc.sum()
    _PSF_CACHE.clear()
    _PSF_CACHE[key] = acc
    return acc


class SpotModel:
    """Spot image on the sensor: PSF x source x jitter x motion x pixel.

    The intensity is sampled on a fine grid of pitch delta = pixel/m, with
    m chosen so that delta <= lambda_min F#/2 (alias-free). The whole PSF
    field (>= 192 lambda/D wide, so the faint diffraction halo is kept the
    same for any motion) sits on an N x N canvas (N = m n_s) centred at
    N//2, enlarged by the blur, source and jitter extents so nothing wraps
    around. A crop is n_s x n_s sensor pixels holding the fraction of the
    spot energy collected by each pixel.
    """

    def __init__(self, p):
        lam_min_f = p.lf_um * (1.0 - p.bandwidth / 2.0)
        self.m = m = max(1, math.ceil(p.pixel_um / (0.5 * lam_min_f) - 1e-9))
        self.delta = p.pixel_um / m
        self.pad0 = p.lf_um / self.delta            # fine samples per lambda0 F#
        # defocus phase step <= pi/2 per pupil sample: n_pupil >= 16 |W020|
        n_psf = _even_fast_len(max(16.0 * abs(p.defocus_waves),
                                   2.0 * MIN_PUPIL_RADIUS) * self.pad0)
        vx, vy = p.velocity_px_s
        speed = math.hypot(vx, vy)
        kernel_px = (p.blur_px                       # extents convolved onto the PSF
                     + (p.source_lod + 8.0 * p.jitter_lod) * p.lf_um / p.pixel_um)
        if p.rolling_shutter:
            kernel_px += speed * p.line_time_us * 1e-6 * n_psf / m
        n_s = _even_fast_len(n_psf / m + 1.15 * kernel_px + 8)
        N = m * n_s
        if max(N, n_psf) > MAX_GRID:
            raise ValueError(
                f"the spot plus its motion blur needs a {max(N, n_psf)}^2 grid "
                f"(limit {MAX_GRID}); reduce defocus, speed x exposure or "
                f"bandwidth, or use larger pixels")
        self.n_psf, self.n_s, self.N = n_psf, n_s, N

        psf = _psf_fine(p, self.pad0, n_psf)
        canvas = np.zeros((N, N), np.float32)
        o = N // 2 - n_psf // 2
        canvas[o:o + n_psf, o:o + n_psf] = psf
        F = sfft.fft2(canvas, workers=-1)
        del canvas
        f = sfft.fftfreq(N).astype(np.float32)          # cycles / fine sample
        pix = np.sinc(np.float32(m) * f)                 # pixel = m x m samples
        F *= pix[:, None]
        F *= pix[None, :]
        src = np.float32(p.source_lod * self.pad0)      # diameter [samples]
        sig = np.float32(p.jitter_lod * self.pad0)      # RMS [samples]
        for i in range(0, N, 512):                       # bounded temporaries
            blk = F[i:i + 512]
            if src > 0 or sig > 0:
                q = np.hypot(f[i:i + 512, None], f[None, :])
                if src > 0:                              # uniform disk 2J1(x)/x
                    x = np.float32(np.pi) * src * q
                    S = np.ones_like(x)
                    nz = x > 1e-12
                    S[nz] = 2.0 * j1(x[nz]) / x[nz]
                    blk *= S
                if sig > 0:                              # Gaussian jitter
                    blk *= np.exp(np.float32(-2.0 * np.pi ** 2) * sig * sig * q * q)
        self.kc = sfft.fftfreq(N) * N                    # centred frequencies
        self.F = F
        self.static_peak = float(self.crop_global(0.0, 0.0)[0].max())
        L = np.float32(p.blur_px * m)                    # blur length [samples]
        if L > 0:                                        # box blur along v
            ux, uy = np.float32(vx / speed), np.float32(vy / speed)
            for i in range(0, N, 512):
                F[i:i + 512] *= np.sinc(L * (f[None, :] * ux + f[i:i + 512, None] * uy))
        self.moving_peak = float(self.crop_global(0.0, 0.0)[0].max())

    def crop_global(self, xc, yc):
        """(crop, r0, c0): energy fraction per pixel of a spot centred at
        (xc, yc) [px]; crop[0, 0] is frame pixel (r0, c0)."""
        m, n_s, N = self.m, self.n_s, self.N
        c0 = int(round(xc)) - n_s // 2
        r0 = int(round(yc)) - n_s // 2
        # shift the canvas so that every m-th fine sample is a pixel centre
        ry = np.exp(-2j * np.pi * self.kc * (m * (yc - r0 - n_s / 2) / N))
        rx = np.exp(-2j * np.pi * self.kc * (m * (xc - c0 - n_s / 2) / N))
        G = self.F * ry.astype(np.complex64)[:, None]
        G *= rx.astype(np.complex64)[None, :]
        G = G.reshape(m, n_s, m, n_s).sum(axis=(0, 2))   # decimate = fold aliases
        return sfft.ifft2(G, workers=-1).real, r0, c0

    def crop_rolling(self, xc, yc, vx, vy, line_s):
        """Rolling shutter: row r is exposed r*line_s later than row 0 and
        sees the spot centred at (xc + vx r line_s, yc + vy r line_s).
        Returns (crop, r0, c0)."""
        m, n_s, N = self.m, self.n_s, self.N
        ax_row, ay_row = vx * line_s, vy * line_s       # px of motion per row
        x_app, y_app = rolling_apparent(xc, yc, vx, vy, line_s)
        c0 = int(round(x_app)) - n_s // 2
        r0 = int(round(y_app)) - n_s // 2
        rows = r0 + np.arange(n_s)
        ny = N / 2 + m * (rows - (yc + ay_row * rows))  # canvas row of each crop row
        nx0 = N / 2 + m * (c0 - (xc + ax_row * rows))   # canvas column of crop column 0
        Ey = np.exp(2j * np.pi * np.outer(ny, self.kc) / N).astype(np.complex64)
        Hm = Ey @ self.F                                 # evaluate rows: (n_s, N)
        Hm *= np.exp(2j * np.pi * np.outer(nx0, self.kc) / N).astype(np.complex64)
        Hm = Hm.reshape(n_s, m, n_s).sum(axis=1)         # fold column aliases
        crop = sfft.ifft(Hm, axis=1, workers=-1).real / n_s
        return crop, r0, c0


def rolling_apparent(xc, yc, vx, vy, line_s):
    """Centre of the spot recorded with a rolling shutter: the row r* that
    sees the spot centre at its own exposure time, r* = yc + vy line_s r*."""
    y_app = yc / (1.0 - vy * line_s)
    return xc + vx * line_s * y_app, y_app


_MODEL_FIELDS = ("d_mm", "obstruction", "vane_mm", "vane_angles_deg", "lambda_nm",
                 "fnum", "defocus_waves", "bandwidth", "source_lod", "jitter_lod",
                 "pixel_um", "exposure_ms", "rolling_shutter", "line_time_us")


def get_model(p):
    """SpotModel for p, reusing the previous one when only non-optical
    settings (noise, intensity, position, frame count ...) changed."""
    key = (tuple(getattr(p, k) for k in _MODEL_FIELDS)
           + (p.speed_px_s, p.direction_deg))
    model = _MODEL_CACHE.get(key)
    if model is None:
        _MODEL_CACHE.clear()
        model = _MODEL_CACHE[key] = SpotModel(p)
    return model


# ------------------------------------------------------------------ simulator
def _paste(frame, crop, r0, c0):
    """Add crop into frame at (r0, c0), clipped to the frame."""
    H, W = frame.shape
    n_r, n_c = crop.shape
    ra, rb = max(r0, 0), min(r0 + n_r, H)
    ca, cb = max(c0, 0), min(c0 + n_c, W)
    if ra < rb and ca < cb:
        frame[ra:rb, ca:cb] += crop[ra - r0:rb - r0, ca - c0:cb - c0]


def to_display(frame, lo, hi, gamma=0.5):
    """DN frame -> uint8 display image (gamma stretch between lo and hi)."""
    v = np.clip((frame.astype(np.float32) - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
    return (255.0 * v ** gamma + 0.5).astype(np.uint8)


class VideoSimulator:
    """Renders frames (DN) and their ground truth for one VideoParams."""

    def __init__(self, params):
        params.validate()
        self.p = p = params
        self.model = get_model(p)
        self.vx, self.vy = p.velocity_px_s
        if p.intensity_mode == "e_per_s":
            self.rate = p.intensity
        else:
            self.rate = (p.intensity / 100.0 * p.full_scale_e
                         / (self.model.static_peak * p.exposure_s))
        self.dtype = np.uint8 if p.bit_depth <= 8 else np.uint16

    def position(self, k):
        """(x, y, jx, jy): true spot centre at the middle of frame k's
        exposure (row 0), including the frame's pointing jitter (jx, jy)."""
        p = self.p
        jx = jy = 0.0
        if p.frame_jitter_px > 0:
            jx, jy = np.random.default_rng([p.seed, k, 1]).normal(
                0.0, p.frame_jitter_px, 2)
        t = k / p.fps
        return p.start_x + self.vx * t + jx, p.start_y + self.vy * t + jy, jx, jy

    def truth_geometry(self, k):
        """Ground-truth geometry of frame k (no rendering): spot centre at
        exposure start / middle / end, where the recorded spot appears
        (differs only with a rolling shutter), jitter offset and blur."""
        p = self.p
        xc, yc, jx, jy = self.position(k)
        xa, ya = xc, yc
        if p.rolling_shutter:
            xa, ya = rolling_apparent(xc, yc, self.vx, self.vy, p.line_time_us * 1e-6)
        half = 0.5 * p.exposure_s
        return dict(frame=k, t_start_s=k / p.fps, t_mid_s=k / p.fps + half,
                    x_mid_px=xc, y_mid_px=yc,
                    x_start_px=xc - self.vx * half, y_start_px=yc - self.vy * half,
                    x_end_px=xc + self.vx * half, y_end_px=yc + self.vy * half,
                    x_apparent_px=xa, y_apparent_px=ya,
                    jitter_dx_px=jx, jitter_dy_px=jy, blur_px=p.blur_px)

    def signal_frame(self, k):
        """(expected spot e- per pixel, crop, truth geometry) for frame k."""
        p, mdl = self.p, self.model
        geo = self.truth_geometry(k)
        xc, yc = geo["x_mid_px"], geo["y_mid_px"]
        if p.rolling_shutter:
            crop, r0, c0 = mdl.crop_rolling(xc, yc, self.vx, self.vy,
                                            p.line_time_us * 1e-6)
        else:
            crop, r0, c0 = mdl.crop_global(xc, yc)
        sig = np.zeros((p.height, p.width))
        _paste(sig, crop, r0, c0)
        return sig * (self.rate * p.exposure_s), crop, geo

    def render_frame(self, k):
        """(frame in DN, ground-truth dict) for frame k (0-based). Noise is
        seeded per frame, so a frame is identical whether rendered alone or
        within a run."""
        p = self.p
        spot_e, crop, truth = self.signal_frame(k)
        e_spot = self.rate * p.exposure_s
        mean_e = spot_e + (p.dark_e_per_s + p.background_e_per_s) * p.exposure_s
        rng = np.random.default_rng([p.seed, k])
        e = rng.poisson(mean_e).astype(np.float64) if p.shot_noise else mean_e
        full = e >= p.full_well_e
        e = np.minimum(e, p.full_well_e)
        if p.read_noise_e > 0:
            e = e + rng.normal(0.0, p.read_noise_e, e.shape)
        dn = np.clip(np.rint(e / p.gain_e_per_dn + p.offset_dn), 0, p.dn_max)
        dn = dn.astype(self.dtype)
        total = float(crop.sum()) * e_spot
        truth.update(signal_e=float(spot_e.sum()),
                     in_frame_fraction=float(spot_e.sum() / total) if total > 0 else 0.0,
                     peak_e=float(spot_e.max()),
                     saturated_px=int(np.count_nonzero(full | (dn >= p.dn_max))))
        return dn, truth

    def display_range(self):
        """(lo, hi) DN for previews: black level .. expected moving peak."""
        p = self.p
        bg = (p.dark_e_per_s + p.background_e_per_s) * p.exposure_s
        peak = self.model.moving_peak * self.rate * p.exposure_s
        top = min(peak + bg, p.full_well_e) + 3.0 * math.sqrt(peak + bg + p.read_noise_e ** 2)
        lo = float(p.offset_dn)
        return lo, max(lo + 1.0, min(float(p.dn_max), lo + top / p.gain_e_per_dn))

    def derived(self):
        """Scales and photometry of this configuration (for readouts)."""
        p, mdl = self.p, self.model
        s, tau = p.speed_px_s, p.exposure_s
        vx, vy = p.velocity_px_s
        e_frame = self.rate * tau
        fs = p.full_scale_e
        d = dict(
            focal_length_mm=p.focal_mm,
            pixel_arcsec=p.pixel_arcsec,
            lambda_over_d_px=p.lf_um / p.pixel_um,
            lambda_over_d_arcsec=p.lf_um / p.pixel_um * p.pixel_arcsec,
            nyquist_ratio=0.5 * p.lf_um / p.pixel_um,
            donut_diameter_px=8.0 * abs(p.defocus_waves) * p.lf_um / p.pixel_um,
            speed_px_s=s, speed_um_s=s * p.pixel_um,
            speed_arcsec_s=s * p.pixel_arcsec, speed_deg_s=s * p.pixel_arcsec / 3600.0,
            vx_px_s=vx, vy_px_s=vy,
            step_px_per_frame=s / p.fps, blur_px=p.blur_px,
            duty_cycle=tau * p.fps,
            signal_e_per_s=self.rate, signal_e_per_frame=e_frame,
            static_peak_e=mdl.static_peak * e_frame,
            static_peak_percent_fs=100.0 * mdl.static_peak * e_frame / fs,
            moving_peak_e=mdl.moving_peak * e_frame,
            moving_peak_percent_fs=100.0 * mdl.moving_peak * e_frame / fs,
            full_scale_e=fs, dn_max=p.dn_max,
            background_e_per_frame=(p.dark_e_per_s + p.background_e_per_s) * tau,
            oversampling=mdl.m, fine_grid=mdl.N, psf_grid=mdl.n_psf, crop_px=mdl.n_s,
            spectral_samples=n_lambda(p.defocus_waves, p.bandwidth))
        if p.rolling_shutter:
            d["rolling_readout_ms"] = p.height * p.line_time_us * 1e-3
        return d

    # ---- whole sequence ------------------------------------------------------
    def frames(self):
        """Yield (k, frame DN, truth) for the whole sequence."""
        for k in range(self.p.frames):
            dn, truth = self.render_frame(k)
            yield k, dn, truth

    def run(self, out=None, progress=None, cancel=None, log=print):
        """Render the sequence into <out folder>/<name>/ and return a summary
        dict. progress(done, total) is called per frame; setting the
        threading.Event `cancel` stops after the current frame."""
        out = out or OutputSpec()
        p = self.p
        folder = _unique_dir(out.run_dir())
        folder.mkdir(parents=True)
        lo, hi = self.display_range()
        writers = []
        if out.npy:
            writers.append(_NpyWriter(folder / "frames.npy", p, self.dtype))
        if out.tiff:
            writers.append(_TiffWriter(folder / "frames.tif"))
        preview = out.preview
        if preview == "mp4":
            try:
                writers.append(_Mp4Writer(folder / "preview.mp4", p, lo, hi))
            except ImportError:
                log("imageio-ffmpeg is not installed: writing preview.gif instead")
                preview = "gif"
        if preview == "gif":
            writers.append(_GifWriter(folder / "preview.gif", p, lo, hi))
        if out.raw != "none":
            writers.append(_RawWriter(folder, out.raw, self.dtype, p.frames))
        mat = _MatWriter(folder, p, self.dtype) if out.mat else None
        if mat is not None:
            writers.append(mat)
        truth, t0, cancelled = [], time.perf_counter(), False
        try:
            for k in range(p.frames):
                if cancel is not None and cancel.is_set():
                    cancelled = True
                    break
                dn, info = self.render_frame(k)
                for w in writers:
                    w.add(dn)
                truth.append(info)
                if progress is not None:
                    progress(k + 1, p.frames)
        finally:
            for w in writers:
                w.close()
        seconds = time.perf_counter() - t0
        meta = dict(simulation=p.to_dict(), output=asdict(out), derived=self.derived(),
                    created=time.strftime("%Y-%m-%d %H:%M:%S"))
        if out.truth:
            _write_truth(folder / "truth.csv", truth)
        if out.raw != "none":
            _write_raw_info(folder, p, out.raw, len(truth), self.dtype)
        if mat is not None:
            mat.finish(truth, meta)
        files = [f"{f.name}/ ({sum(1 for _ in f.iterdir())} files)" if f.is_dir() else f.name
                 for f in sorted(folder.iterdir())] + ["params.json"]
        summary = dict(folder=str(folder), files=files, frames_written=len(truth),
                       cancelled=cancelled, render_seconds=round(seconds, 3),
                       preview_display_range_dn=[lo, hi])
        meta.update(summary)
        (folder / "params.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return summary


def _unique_dir(path):
    """path, or path_2, path_3 ... if it already exists (never overwrite)."""
    if not path.exists():
        return path
    i = 2
    while path.with_name(f"{path.name}_{i}").exists():
        i += 1
    return path.with_name(f"{path.name}_{i}")


def _write_truth(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        if not rows:
            return
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v)
                        for k, v in r.items()})


class _NpyWriter:
    def __init__(self, path, p, dtype):
        self.path, self.k = path, 0
        self.mm = np.lib.format.open_memmap(path, mode="w+", dtype=dtype,
                                            shape=(p.frames, p.height, p.width))

    def add(self, frame):
        self.mm[self.k] = frame
        self.k += 1

    def close(self):
        self.mm.flush()
        if self.k < self.mm.shape[0]:          # cancelled: keep written frames
            part = np.array(self.mm[:self.k])
            del self.mm
            np.save(self.path, part)
        else:
            del self.mm


class _TiffWriter:
    def __init__(self, path):
        from PIL import TiffImagePlugin
        self.tf = TiffImagePlugin.AppendingTiffWriter(str(path), True)

    def add(self, frame):
        from PIL import Image
        Image.fromarray(frame).save(self.tf, format="TIFF")
        self.tf.newFrame()

    def close(self):
        self.tf.close()


class _Mp4Writer:
    def __init__(self, path, p, lo, hi):
        import imageio_ffmpeg
        self.lo, self.hi = lo, hi
        self.shape = (p.height + p.height % 2, p.width + p.width % 2)  # yuv420p
        self.gen = imageio_ffmpeg.write_frames(
            str(path), (self.shape[1], self.shape[0]), pix_fmt_in="gray",
            fps=min(p.fps, PREVIEW_MAX_FPS), quality=8, macro_block_size=1)
        self.gen.send(None)

    def add(self, frame):
        img = to_display(frame, self.lo, self.hi)
        if img.shape != self.shape:
            img = np.pad(img, ((0, self.shape[0] - img.shape[0]),
                               (0, self.shape[1] - img.shape[1])))
        self.gen.send(np.ascontiguousarray(img).tobytes())

    def close(self):
        self.gen.close()


class _GifWriter:
    def __init__(self, path, p, lo, hi):
        self.path, self.lo, self.hi = path, lo, hi
        self.step = max(1, math.ceil(p.frames / GIF_MAX_FRAMES))
        self.ms = max(20, int(round(1000.0 * self.step / min(p.fps, 50.0))))
        self.images, self.k = [], 0

    def add(self, frame):
        from PIL import Image
        if self.k % self.step == 0:
            self.images.append(Image.fromarray(to_display(frame, self.lo, self.hi)))
        self.k += 1

    def close(self):
        if self.images:
            self.images[0].save(self.path, save_all=True, append_images=self.images[1:],
                                duration=self.ms, loop=0)


class _RawWriter:
    """Headerless little-endian binary frames: one file per frame in raw/,
    or all frames back to back in frames.raw."""

    def __init__(self, folder, kind, dtype, n_frames):
        self.dt = np.dtype(dtype).newbyteorder("<")
        self.kind, self.k = kind, 0
        self.digits = max(5, len(str(max(n_frames - 1, 0))))
        if kind == "per_frame":
            self.dir = folder / "raw"
            self.dir.mkdir()
        else:
            self.f = open(folder / "frames.raw", "wb")

    def add(self, frame):
        data = np.ascontiguousarray(frame, dtype=self.dt)
        if self.kind == "per_frame":
            data.tofile(self.dir / f"frame_{self.k:0{self.digits}d}.raw")
        else:
            self.f.write(data.tobytes())
        self.k += 1

    def close(self):
        if self.kind != "per_frame":
            self.f.close()


def _write_raw_info(folder, p, kind, n, dtype):
    """raw_format.txt: layout of the raw files and how to read them."""
    size = np.dtype(dtype).itemsize
    bits, W, H = 8 * size, p.width, p.height
    mtype, ptype = f"uint{bits}", f"<u{size}"
    if kind == "per_frame":
        d = max(5, len(str(max(p.frames - 1, 0))))
        first = f"raw/frame_{0:0{d}d}.raw"
        where = (f"files  : {first} ... raw/frame_{max(n - 1, 0):0{d}d}.raw "
                 f"({n} files, one frame of {W * H * size} bytes each)")
        matlab = (f"fid = fopen('{first}'); img = fread(fid, [{W} {H}], "
                  f"'{mtype}=>{mtype}')'; fclose(fid);")
        python = f"img = np.fromfile('{first}', dtype='{ptype}').reshape({H}, {W})"
        imagej = (f"File > Import > Raw...: {bits}-bit Unsigned, width {W}, "
                  f"height {H}, little-endian byte order")
    else:
        where = (f"file   : frames.raw ({n} frames back to back, "
                 f"{W * H * size} bytes each, no header)")
        matlab = (f"fid = fopen('frames.raw'); v = fread(fid, Inf, '{mtype}=>{mtype}'); "
                  f"fclose(fid); frames = permute(reshape(v, {W}, {H}, []), [2 1 3]);")
        python = f"frames = np.fromfile('frames.raw', dtype='{ptype}').reshape(-1, {H}, {W})"
        imagej = (f"File > Import > Raw...: {bits}-bit Unsigned, width {W}, height {H}, "
                  f"number of images {n}, little-endian byte order")
    lines = [
        "Raw frames written by cassegrain_video.py",
        f"format : headerless binary, unsigned {bits}-bit integers, little-endian,",
        "         row-major (each row left to right, row 0 = top of the image)",
        f"size   : width {W} x height {H} pixels",
        where,
        f"values : camera counts (DN), {p.bit_depth}-bit right-aligned "
        f"(0 .. {p.dn_max}), black level {p.offset_dn:g} DN",
        "",
        "Read with",
        f"  MATLAB : {matlab}",
        f"  Python : {python}",
        f"  ImageJ : {imagej}",
        "",
        "Ground truth per frame: truth.csv; all settings: params.json."]
    (folder / "raw_format.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


MAT_INFO = ("frames(:,:,k) is frame first_frame+k-1 in camera counts (DN). truth has "
            "one value per frame (truth.frame is the 0-based frame index). Pixel "
            "coordinates are 0-based with (0,0) at the centre of the top-left pixel: "
            "MATLAB row = y + 1, column = x + 1. params and derived hold all settings.")


class _MatWriter:
    """Streams frames to a temporary file; finish() writes frames.mat in
    MATLAB layout (H x W x N) with the truth table and the settings. Above
    MAT_PART_BYTES the frames are split over frames_part001.mat, ... (the
    MAT v5 format that scipy writes holds < 2 GB per variable)."""

    def __init__(self, folder, p, dtype):
        self.folder, self.shape, self.dtype = folder, (p.height, p.width), np.dtype(dtype)
        self.tmp = folder / ".frames_mat.tmp"
        self.f = open(self.tmp, "wb")
        self.k = 0

    def add(self, frame):
        self.f.write(np.ascontiguousarray(frame, dtype=self.dtype).tobytes())
        self.k += 1

    def close(self):
        self.f.close()

    def finish(self, truth, meta):
        from scipy.io import savemat
        H, W = self.shape
        per_part = max(1, MAT_PART_BYTES // (H * W * self.dtype.itemsize))
        common = dict(info=MAT_INFO, params=meta["simulation"], derived=meta["derived"])
        if truth:
            common["truth"] = {key: np.array([row[key] for row in truth]) for key in truth[0]}
        data = None
        try:
            if self.k:
                data = np.memmap(self.tmp, dtype=self.dtype, mode="r", shape=(self.k, H, W))
                starts = range(0, self.k, per_part)
                for i, a in enumerate(starts):
                    name = "frames.mat" if len(starts) == 1 else f"frames_part{i + 1:03d}.mat"
                    savemat(self.folder / name,
                            dict(common, frames=np.moveaxis(data[a:a + per_part], 0, -1),
                                 first_frame=float(a)),
                            do_compression=False, oned_as="column")
        finally:
            del data                      # release the mapping before deleting
            try:
                self.tmp.unlink()
            except OSError:
                pass


# ------------------------------------------------------------------- selftest
def selftest():
    """Headless physics and I/O checks; returns 0 if all pass."""
    import tempfile
    t_start = time.perf_counter()
    fails = []

    def check(name, ok, detail):
        print(f"  {'OK ' if ok else 'FAIL'} {name:34s} {detail}")
        if not ok:
            fails.append(name)

    def moments(img):
        y, x = np.indices(img.shape, dtype=np.float64)
        s = img.sum()
        cx, cy = (img * x).sum() / s, (img * y).sum() / s
        cov = np.array([[(img * (x - cx) ** 2).sum(), (img * (x - cx) * (y - cy)).sum()],
                        [(img * (x - cx) * (y - cy)).sum(), (img * (y - cy) ** 2).sum()]]) / s
        return cx, cy, cov

    quiet = dict(shot_noise=False, read_noise_e=0.0, dark_e_per_s=0.0,
                 background_e_per_s=0.0, width=256, height=256)

    # energy and sub-pixel centroid of a static spot
    sim = VideoSimulator(VideoParams(speed=0.0, defocus_waves=6.0, **quiet))
    rng = np.random.default_rng(7)
    errs, sums = [], []
    for _ in range(6):
        xc, yc = rng.uniform(100, 156, 2)
        crop, r0, c0 = sim.model.crop_global(xc, yc)
        cx, cy, _ = moments(crop)
        errs.append(math.hypot(cx + c0 - xc, cy + r0 - yc))
        sums.append(crop.sum())
    check("energy conservation", max(abs(s - 1) for s in sums) < 1e-4,
          f"sum of pixel fractions = {min(sums):.6f} .. {max(sums):.6f}")
    check("sub-pixel centroid", max(errs) < 0.01, f"max error {max(errs):.4f} px")

    # motion blur adds L^2/12 along the motion direction to the covariance
    base = VideoParams(speed=0.0, defocus_waves=6.0, **quiet)
    _, _, cov0 = moments(VideoSimulator(base).model.crop_global(0.3, -0.2)[0])
    for direction in (0.0, 90.0, 30.0):
        pm = VideoParams(**{**base.to_dict(), "speed": 2000.0, "direction_deg": direction})
        L, (vx, vy) = pm.blur_px, pm.velocity_px_s
        u = np.array([vx, vy]) / math.hypot(vx, vy)
        _, _, cov = moments(VideoSimulator(pm).model.crop_global(0.3, -0.2)[0])
        expect = L * L / 12.0 * np.outer(u, u)
        err = np.abs(cov - cov0 - expect).max() / (L * L / 12.0)
        check(f"blur covariance, dir {direction:g} deg", err < 0.01,
              f"L = {L:.0f} px, relative error {err:.4f}")

    # rolling shutter: zero line time equals the global shutter; otherwise
    # the recorded spot is centred where the analytic model puts it
    pr = VideoParams(**{**base.to_dict(), "speed": 3000.0, "direction_deg": 60.0,
                        "rolling_shutter": True, "line_time_us": 30.0})
    mdl = SpotModel(pr)
    vx, vy = pr.velocity_px_s
    g, r0, c0 = mdl.crop_global(128.4, 120.7)
    r, r1, c1 = mdl.crop_rolling(128.4, 120.7, vx, vy, 0.0)
    check("rolling (0 line time) == global", (r0, c0) == (r1, c1) and
          np.abs(g - r).max() < 1e-6 * g.max(), f"max diff {np.abs(g - r).max():.2e}")
    crop, r0, c0 = mdl.crop_rolling(128.4, 120.7, vx, vy, 30e-6)
    xa, ya = rolling_apparent(128.4, 120.7, vx, vy, 30e-6)
    cx, cy, _ = moments(crop)
    err = math.hypot(cx + c0 - xa, cy + r0 - ya)
    check("rolling-shutter apparent centre", err < 0.05,
          f"shift ({xa - 128.4:+.2f}, {ya - 120.7:+.2f}) px, error {err:.4f} px")

    # static-peak intensity calibration and ADC
    ps = VideoParams(speed=0.0, intensity_mode="static_peak", intensity=40.0,
                     start_x=128.0, start_y=128.0, **quiet)
    s = VideoSimulator(ps)
    dn, info = s.render_frame(0)
    want = round(0.40 * ps.full_scale_e / ps.gain_e_per_dn + ps.offset_dn)
    check("static-peak calibration", abs(int(dn.max()) - want) <= 1,
          f"peak {int(dn.max())} DN, expected {want} DN")

    # noise statistics on a flat field (spot outside the frame)
    pn = VideoParams(start_x=-5000.0, speed=0.0, width=200, height=200,
                     background_e_per_s=20000.0, dark_e_per_s=1000.0,
                     read_noise_e=5.0, gain_e_per_dn=2.0, exposure_ms=20.0)
    dn, _ = VideoSimulator(pn).render_frame(3)
    mean_e = (20000.0 + 1000.0) * pn.exposure_s
    m_th = mean_e / pn.gain_e_per_dn + pn.offset_dn
    v_th = (mean_e + pn.read_noise_e ** 2) / pn.gain_e_per_dn ** 2 + 1.0 / 12.0
    m, v = dn.mean(), dn.astype(float).var()
    check("noise mean / variance", abs(m - m_th) < 0.2 and abs(v / v_th - 1) < 0.03,
          f"mean {m:.2f} (th {m_th:.2f}), var {v:.2f} (th {v_th:.2f}) DN")

    # writers and per-frame reproducibility
    pw = VideoParams(width=96, height=80, frames=6, speed=400.0, start_x=30.0,
                     start_y=40.0, defocus_waves=3.0, frame_jitter_px=0.5)
    sw = VideoSimulator(pw)
    with tempfile.TemporaryDirectory() as tmp:
        quiet_log = lambda *a: None
        res = sw.run(OutputSpec(name="t", folder=tmp, preview="gif", mat=True,
                                raw="per_frame"), log=quiet_log)
        run_dir = Path(res["folder"])
        arr = np.load(run_dir / "frames.npy")
        from PIL import Image
        im = Image.open(run_dir / "frames.tif")   # explicit close: a seeked
        try:                                       # TIFF stays open after 'with'
            n_tif = im.n_frames
            im.seek(4)
            tif4 = np.array(im)
        finally:
            im.close()
        with open(run_dir / "truth.csv", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        meta = json.loads((run_dir / "params.json").read_text(encoding="utf-8"))
        again, _ = sw.render_frame(4)
        ok = (arr.shape == (6, 80, 96) and arr.dtype == np.uint16 and n_tif == 6
              and np.array_equal(tif4, arr[4]) and len(rows) == 6
              and (run_dir / "preview.gif").exists()
              and meta["frames_written"] == 6 and np.array_equal(again, arr[4]))
        check("writers + reproducibility", ok,
              f"npy {arr.shape}, tif {n_tif} pages, csv {len(rows)} rows, gif")

        from scipy.io import loadmat
        mat = loadmat(run_dir / "frames.mat")
        mx = mat["truth"]["x_mid_px"][0, 0][:, 0]
        ok = (mat["frames"].shape == (80, 96, 6) and mat["frames"].dtype == np.uint16
              and np.array_equal(np.moveaxis(mat["frames"], -1, 0), arr)
              and np.allclose(mx, [float(r["x_mid_px"]) for r in rows])
              and float(mat["params"]["speed"][0, 0][0, 0]) == pw.speed
              and not (run_dir / ".frames_mat.tmp").exists())
        check("frames.mat (H x W x N, truth, params)", ok,
              f"frames {mat['frames'].shape} {mat['frames'].dtype}")

        per = [np.fromfile(run_dir / "raw" / f"frame_{k:05d}.raw", dtype="<u2")
               .reshape(80, 96) for k in range(6)]
        res2 = sw.run(OutputSpec(name="s", folder=tmp, npy=False, tiff=False,
                                 preview="none", truth=False, raw="single"), log=quiet_log)
        single = np.fromfile(Path(res2["folder"]) / "frames.raw", dtype="<u2").reshape(-1, 80, 96)
        ok = (np.array_equal(np.stack(per), arr) and np.array_equal(single, arr)
              and (run_dir / "raw_format.txt").exists())
        check("raw frames (per frame, single file)", ok,
              f"{len(per)} files + frames.raw {single.shape}")

        global MAT_PART_BYTES
        saved_part = MAT_PART_BYTES
        MAT_PART_BYTES = 2 * 80 * 96 * 2           # force 2 frames per .mat part
        try:
            res3 = sw.run(OutputSpec(name="m", folder=tmp, npy=False, tiff=False,
                                     preview="none", mat=True), log=quiet_log)
        finally:
            MAT_PART_BYTES = saved_part
        d3 = Path(res3["folder"])
        parts = sorted(d3.glob("frames_part*.mat"))
        loaded = [loadmat(f) for f in parts]
        joined = np.concatenate([m["frames"] for m in loaded], axis=-1)
        ok = (len(parts) == 3 and [float(m["first_frame"][0, 0]) for m in loaded] == [0, 2, 4]
              and np.array_equal(np.moveaxis(joined, -1, 0), arr))
        check("frames.mat split into parts", ok, f"{len(parts)} parts, first frames 0, 2, 4")
    print(f"selftest {'passed' if not fails else 'FAILED: ' + ', '.join(fails)} "
          f"({time.perf_counter() - t_start:.1f} s)")
    return 1 if fails else 0


# ------------------------------------------------------------------------ CLI
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Cassegrain video simulator (headless). Without a config "
                    "file the defaults are used.",
        epilog="examples: python cassegrain_video.py run.json --set speed=500 "
               "--set direction_deg=45 --name fast_pass | "
               "python cassegrain_video.py --set mat=true --set raw=per_frame "
               "(raw: none, per_frame, single)")
    ap.add_argument("config", nargs="?", help="JSON configuration (see --write-config)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override a simulation or output value (repeatable)")
    ap.add_argument("--name", help="run name = output sub-folder")
    ap.add_argument("--out", help="parent output folder (default: <project>/output)")
    ap.add_argument("--write-config", metavar="PATH",
                    help="write the configuration (defaults + --set) to PATH and exit")
    ap.add_argument("--selftest", action="store_true", help="run headless checks")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    try:
        params, out = load_config(args.config) if args.config else (VideoParams(), OutputSpec())
        apply_overrides(params, out, args.set)
        if args.name:
            out.name = args.name
        if args.out:
            out.folder = args.out
        if args.write_config:
            params.validate()
            save_config(args.write_config, params, out)
            print(f"wrote {args.write_config}")
            return 0
        t0 = time.perf_counter()
        sim = VideoSimulator(params)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    d = sim.derived()
    print(f"{params.frames} frames {params.width}x{params.height} px, pixel "
          f"{params.pixel_um:g} um = {d['pixel_arcsec']:.3f} arcsec, lambda/D = "
          f"{d['lambda_over_d_px']:.2f} px (model ready in {time.perf_counter() - t0:.1f} s)")
    print(f"speed {d['speed_px_s']:.1f} px/s = {d['speed_arcsec_s']:.1f} arcsec/s -> "
          f"{d['step_px_per_frame']:.2f} px/frame, blur {d['blur_px']:.2f} px per exposure")
    print(f"signal {d['signal_e_per_frame']:.4g} e-/frame; peak {d['static_peak_percent_fs']:.0f}% "
          f"of full scale static, {d['moving_peak_percent_fs']:.0f}% moving")

    def progress(done, total):
        if done == total or done % max(1, total // 10) == 0:
            print(f"  frame {done}/{total}")

    res = sim.run(out, progress=progress)
    print(f"saved {res['frames_written']} frames in {res['render_seconds']:.1f} s to "
          f"{res['folder']}: {', '.join(res['files'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
