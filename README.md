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
