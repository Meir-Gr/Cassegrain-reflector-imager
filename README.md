# Cassegrain reflector imager

Scalar Fourier-optics study of defocused ("donut") star images in a
Cassegrain telescope with a 3-vane spider: point and extended sources,
vibration, obstruction and vane-width effects, plus an interactive simulator
and a camera video simulator of a spot in linear motion.

## Layout

| Path | Contents |
| --- | --- |
| `studies/` | Batch simulation scripts; each one writes its figures to `docs/figures/` |
| `simulator/` | Interactive simulator: `cassegrain_gui.py` (Python/Tkinter) and its MATLAB port `cassegrain_simulator.m`; video simulator: `cassegrain_video.py` (engine, command line, Python API) and `cassegrain_video_gui.py` (GUI), with the MATLAB ports `cassegrain_video.m` and `cassegrain_video_gui.m` |
| `docs/` | `FINDINGS` (results) and `THEORY` (background) reports, Markdown sources + PDFs |
| `docs/figures/` | All generated figures (most are embedded in `FINDINGS.md`) |
| `tools/` | Report pipeline: Markdown-to-PDF builder, PDF table-layout check, GUI screenshot |
| `output/` | Video-simulator runs, one sub-folder per run (created on the first run; not tracked by git) |
| `.cache/` | Disposable math-render cache of the PDF builder (not tracked by git) |

## Study scripts

`cassegrain_defocus_sim.py` is the shared optics core (system parameters,
pupil, PSF) and is imported by every other study.

| Script | Topic (FINDINGS section) | Imported by |
| --- | --- | --- |
| `cassegrain_defocus_sim.py` | Point source through focus (§2) | all other studies |
| `cassegrain_extended_source_sim.py` | Extended sources (§3) | vibration, Arago, central-spot, SA scripts |
| `cassegrain_arago_vs_source.py` | Central (Arago) spot vs source size (§3.4) | — |
| `cassegrain_vibration_sim.py` | Vibration blur (§4) | — |
| `cassegrain_parameter_study.py` | Donut size, obstruction, vane width (§5) | Arago, SA scripts |
| `cassegrain_central_spot_cases.py` | Four causes of a central spot (side study, not in FINDINGS) | — |
| `cassegrain_sa_central_spot.py` | Spherical aberration filling the hole (side study, not in FINDINGS) | — |

## Running (from this folder)

```
python -m pip install -r requirements.txt
python studies/cassegrain_defocus_sim.py        # or any other study
python simulator/cassegrain_gui.py              # interactive simulator
python simulator/cassegrain_gui.py --selftest   # headless check
python simulator/cassegrain_video_gui.py        # video simulator (GUI)
python simulator/cassegrain_video.py run.json   # video simulator (script)
python simulator/cassegrain_video.py --selftest # headless check
python tools/make_gui_screenshot.py             # refresh the GUI figure
python tools/build_report.py                    # docs/FINDINGS.md -> PDF
python tools/build_report.py THEORY.md          # docs/THEORY.md -> PDF
python tools/verify_pdf_layout.py               # table pagination check
```

MATLAB: add `simulator/` to the path, then run `cassegrain_simulator`
(or `cassegrain_simulator selftest`). The video simulator has MATLAB ports
too (R2020a or newer, no toolboxes):

```
cassegrain_video_gui                       % video simulator (GUI)
cassegrain_video run.json --set speed=500  % video simulator (script)
cassegrain_video --selftest                % headless check
```

## Wavelength bands, including MWIR

In either still-image GUI enter **Lower wavelength [µm]** and **Upper
wavelength [µm]**, then **Apply system parameters**. The video GUIs have the
same fields on the **Optics** tab. Use `3.5` and `5.5` for the MWIR example,
or enter the same value twice for monochromatic light (`0.55` is the original
550 nm default). These are passband edges, not centre and half-width.

Video JSON and command-line settings retain their original units: `lambda_nm`
is the arithmetic midpoint in **nm**; `bandwidth` is full width divided by
that midpoint. For 3.5–5.5 µm use `lambda_nm=4500` and
`bandwidth=0.4444444444444444`. Wider bands are supported; an explicit error
is reported if the band and defocus need too large a grid or too many spectral
samples. Wide, strongly defocused bands take longer to render.

To write a video configuration for this band without rendering frames, run
the following from the repository root. Set the camera pixel pitch and noise
levels in `mwir.json` to match your detector before running the simulation.

```sh
python -m simulator.cassegrain_video --set lambda_nm=4500 --set bandwidth=0.4444444444444444 --write-config mwir.json
```

**What transfers to IR:** the on-axis diffraction pattern of an ideal
all-mirror telescope. The angular scale is λ/D, and the focal-plane scale
is λ × F-number. At f/10 the 3.5 µm edge requires a pitch of **17.5 µm or
finer** to Nyquist-sample the full band. This is the still simulator's
computational pitch, not a requirement for running the video model: enter
the camera's actual pixel pitch there to see the effect of undersampling.
When comparing bands at the same physical focus offset, keep Δz fixed,
**not** the number of waves at the midpoint.

**What does not transfer automatically:** real MWIR spectral response or
sensitivity. Equal wavelength intervals contribute equal integrated signal
here. For a thermal source and a real photon detector the weights depend on
source temperature/emissivity, atmospheric and optical transmission, and
detector quantum efficiency. Warm optics and sky also *add* background.
The video model accepts measured signal, background and dark-current rates;
it does not predict them from temperature or from the band edges. See
[the polychromatic model and its limits](docs/THEORY.md#73-polychromatic-light-implemented-in-the-interactive-simulators).

### Python still-image preview performance

The speedup is automatic; it does not drop wavelengths or coarsen the grid.
The simulator builds each pupil only within its bounding box. For a point
source with no jitter or motion, the GUI computes only the displayed Fourier
bins. If the brightest pixel might be outside the crop, it uses the full
field to preserve peak normalization. For an extended source or any blur,
it blurs the full field **before** cropping so light from outside the view
is included. `OpticalModel.image()` still returns the entire grid.

For example, at **3.5–5.5 µm, f/10, and 8 waves of defocus**, the unblurred
point-source preview went from **1.80 s to 0.24 s** on the benchmark machine
(7.4×; median of five fresh calculations, excluding startup and screen
drawing). Both versions used 14 wavelengths on the same 1344 × 1344 grid;
the largest displayed pixel difference was under 0.000001 on a 0–1 scale.
This timing is for the Python still-image GUI, not MATLAB or the video simulator.

## Video simulator

Renders a frame sequence of the spot as a real camera records it, with the
same optics as the interactive simulator (defocus, bandwidth, source size,
in-frame jitter) plus:

- **Sensor:** image size H × W, pixel pitch, frame rate, exposure time. Each
  pixel integrates the light over its area, with sub-pixel spot positions;
  undersampled sensors alias exactly as in a real camera.
- **Motion:** constant velocity (speed + direction) from a start point.
  Blur per frame = speed × exposure; step between frames = speed ÷ frame
  rate. Speed can be given in px/s, µm/s, arcsec/s or deg/s.
- **Noise and read-out:** shot noise, dark current, background, read noise,
  full-well saturation, gain, black level and bit depth. Options: rolling
  shutter (line time) and random frame-to-frame pointing jitter.
- **Intensity:** either the photoelectron rate of the whole spot [e⁻/s], or
  the peak a *static* spot would reach [% of full scale]. Both are converted
  to the same rate, so motion blur lowers the peak as in a real camera.

The three ways to run it share one JSON configuration format:

- **GUI:** `python simulator/cassegrain_video_gui.py`. Besides the full frame
  and a zoom on the spot, it has an exact pixel view: each camera pixel is
  drawn as a uniform square block, never resampled (1:1 up to 32×, with a
  pixel grid and a µm scale bar). Use the mouse wheel to zoom and drag to
  pan; hovering shows the pixel value in DN. The *Save config* button
  writes the file that the command line reads.
- **Command line:** write the defaults with `python simulator/cassegrain_video.py
  --write-config run.json`, edit them, then run `python simulator/cassegrain_video.py
  run.json`. Single values can be changed with `--set key=value`, for example
  `--set speed=500 --set mat=true --set raw=per_frame`.
- **Python:** see the example below.

```
import sys; sys.path.insert(0, "simulator")
from cassegrain_video import VideoParams, OutputSpec, VideoSimulator
sim = VideoSimulator(VideoParams(width=1280, height=1024, pixel_um=5.86,
                                 frames=200, fps=100, exposure_ms=5,
                                 speed=2000, direction_deg=30))
frame_dn, truth = sim.render_frame(0)               # one frame + ground truth
sim.run(OutputSpec(name="fast_pass", mat=True, raw="single"))
```

Each run is saved to `output/<name>/`. The output options are:

- `frames.npy` and `frames.tif`: lossless camera counts (DN).
- `frames.mat`: optional. Holds MATLAB `frames` (H × W × N) plus `truth`,
  `params` and `derived`; files over 1 GiB are split into numbered parts.
- Raw binary frames: optional. Headerless little-endian, either one file per
  frame in `raw/` or all frames in `frames.raw`. `raw_format.txt` describes
  the layout and gives MATLAB, Python and ImageJ readers.
- `preview.mp4`: needs `imageio-ffmpeg`; without it a GIF is written instead.
- `truth.csv` and `params.json`.

`truth.csv` gives the true spot centre per frame: at the start, middle and
end of the exposure, and where the recorded spot appears (this differs with
a rolling shutter).

Conventions:

- Pixel coordinates are x = column and y = row, with (0, 0) at the centre of
  the top-left pixel.
- The start point is the spot centre at the middle of the first exposure.
- Direction is counter-clockwise on screen: 0° = right, 90° = up.
- Noise is seeded per frame, so any frame can be reproduced alone.

Model limits:

- The PSF is the same over the whole sensor (on-axis model, no field-dependent
  coma).
- Motion is a straight line at constant speed.
- Frame-to-frame jitter is uncorrelated between frames.
- Noise-free image centroids agree with the true centre within 0.01 px.

### MATLAB version

`cassegrain_video.m` and `cassegrain_video_gui.m` are MATLAB ports with the
same settings, JSON configuration files, output files and command-line
options, so a configuration saved by one version runs in the other.

- Noise-free frames and the ground truth match the Python version (frames
  can differ by 1 DN in a few pixels that sit on a rounding boundary).
- Noise has the same statistics but uses MATLAB's random generator, so a noisy
  frame is not bit-identical to the Python frame with the same seed.
- MP4 previews need VideoWriter's MPEG-4 profile (Windows, macOS). Elsewhere a
  GIF preview is written instead.
- The GUI renders in MATLAB's own thread. During a run it stays responsive
  between frames, so Cancel works. An edited field takes effect when you press
  Enter or leave the field.

From MATLAB code:

```
cv = cassegrain_video('api');
sim = cv.VideoSimulator(cv.VideoParams('width', 1280, 'height', 1024, ...
    'pixel_um', 5.86, 'frames', 200, 'fps', 100, 'exposure_ms', 5, ...
    'speed', 2000, 'direction_deg', 30));
[frame_dn, truth] = sim.render_frame(0);   % one frame + ground truth
sim.run(cv.OutputSpec('name', 'fast_pass', 'mat', true, 'raw', 'single'));
```
