function varargout = cassegrain_video(varargin)
%CASSEGRAIN_VIDEO Cassegrain video simulator: a camera recording of a spot in linear motion.
%
% MATLAB port of cassegrain_video.py: same model, settings, JSON
% configuration, output files and command line. Every frame is the spot
% image of the interactive simulator (scalar diffraction: annular pupil
% with a 3-vane spider, defocus W020, flat-spectrum bandwidth, uniform-disk
% source, Gaussian in-frame jitter) as a real sensor records it:
%
%   * constant-velocity linear motion: during one exposure the spot travels
%     L = speed x exposure (an exact box blur, MTF sinc(L q_par) along the
%     motion); between frames it advances speed / fps;
%   * pixel integration (box MTF of the pixel pitch) and sampling on the
%     sensor grid, with the aliasing of undersampled sensors included
%     exactly (spectrum folding); sub-pixel spot positions by Fourier shift;
%   * photon shot noise, dark current, uniform background, full-well
%     saturation, read noise, gain / offset / bit-depth quantization;
%   * optional rolling shutter (row r exposed r x line time later) and
%     random frame-to-frame pointing jitter.
%
% Conventions: x = column, y = row in pixels (origin at the centre of the
% top-left pixel, y down; MATLAB row = y + 1, column = x + 1). (start_x,
% start_y) is the spot centre at the middle of the first frame's exposure
% (rolling shutter: of row 0). The motion direction is counter-clockwise
% from +x as seen on screen (0 = right, 90 = up). Frame k (0-based) is
% exposed from k/fps to k/fps + exposure. The PSF is the same everywhere on
% the sensor (on-axis model).
%
% Intensity: either the photoelectron rate of the whole spot [e-/s] (signal
% per frame = rate x exposure, so motion blur lowers the peak), or the peak
% a static spot would reach [% of full scale], converted to the same rate.
%
% Outputs (output/<name>/): frames.npy, frames.tif, frames.mat (frames
% H x W x N plus truth / params), raw frames (headerless little-endian
% binary, one file per frame in raw/ or all frames in frames.raw, layout in
% raw_format.txt), preview.mp4 / .gif, truth.csv (true position per frame)
% and params.json (all settings). Noise-free frames match the Python
% version; noise has the same statistics but comes from MATLAB's generator
% (seeded per frame, so any frame can be reproduced on its own).
%
% Command line (the options of cassegrain_video.py):
%   cassegrain_video [run.json] [--set key=value ...] [--name NAME] [--out FOLDER]
%   cassegrain_video --write-config run.json    (defaults + --set, to edit)
%   cassegrain_video --selftest
%   status = cassegrain_video(...)              0 = ok, 1 = selftest failed, 2 = error
%   cassegrain_video_gui                        (interactive)
%
% MATLAB API (the counterpart of the Python classes):
%   cv  = cassegrain_video('api');
%   sim = cv.VideoSimulator(cv.VideoParams('width', 1280, 'height', 1024, ...
%           'pixel_um', 5.86, 'frames', 200, 'fps', 100, 'exposure_ms', 5, ...
%           'speed', 2000, 'direction_deg', 30));
%   [frame_dn, truth] = sim.render_frame(0);       % one frame + ground truth
%   sim.run(cv.OutputSpec('name', 'fast_pass', 'mat', true, 'raw', 'single'));
%
% Requirements: MATLAB R2020a or newer, no toolboxes. MP4 previews use
% VideoWriter's MPEG-4 profile (Windows, macOS); elsewhere preview.gif is
% written instead.

if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'api')
    varargout{1} = api();
    return
end
status = main(varargin);
if nargout > 0
    varargout{1} = status;
end
end

% ================================================================ interface
function a = api()
%API Function handles of the engine (used by cassegrain_video_gui).
C = consts();
a.ROOT = C.ROOT;
a.DEFAULT_OUTPUT = C.DEFAULT_OUTPUT;
a.SPEED_UNITS = C.SPEED_UNITS;
a.INTENSITY_MODES = C.INTENSITY_MODES;
a.PREVIEW_KINDS = C.PREVIEW_KINDS;
a.RAW_KINDS = C.RAW_KINDS;
spec = param_spec();
a.PARAM_FIELDS = spec(:, 1)';
spec = output_spec();
a.OUTPUT_FIELDS = spec(:, 1)';
a.VideoParams = @VideoParams;           % p = VideoParams('speed', 500, ...) or (struct)
a.OutputSpec = @OutputSpec;
a.validate = @validate_params;
a.to_dict = @params_to_dict;            % typed struct as saved to JSON / .mat
a.output_dict = @output_to_dict;
a.scales = @scales;                     % derived scales (lf_um, pixel_arcsec, speed_px_s ...)
a.load_config = @load_config;           % [p, out] = load_config(path)
a.save_config = @save_config;           % save_config(path, p, out)
a.apply_overrides = @apply_overrides;   % [p, out] = apply_overrides(p, out, {'key=value'})
a.VideoSimulator = @VideoSimulator;
a.SpotModel = @spot_model;
a.get_model = @get_model;
a.crop_global = @crop_global;
a.crop_rolling = @crop_rolling;
a.rolling_apparent = @rolling_apparent;
a.n_lambda = @n_lambda;
a.even_fast_len = @even_fast_len;
a.to_display = @to_display;
a.read_npy = @read_npy;                 % H x W x N array of a frames.npy
a.npy_open = @npy_open;                 % header info of a frames.npy
a.npy_frame = @npy_frame;               % one frame (0-based) of a frames.npy
a.json_dumps = @json_dumps;
a.error_text = @error_text;             % 'ValueError: message' for an MException
a.selftest = @selftest;
a.main = @main;
end

function C = consts()
persistent K
if isempty(K)
    here = fileparts(mfilename('fullpath'));
    K.ROOT = fileparts(here);
    K.DEFAULT_OUTPUT = fullfile(K.ROOT, 'output');
    K.RAD_TO_ARCSEC = 180.0 / pi * 3600.0;
    K.SPEED_UNITS = {'px/s', 'um/s', 'arcsec/s', 'deg/s'};
    K.INTENSITY_MODES = {'static_peak', 'e_per_s'};
    K.PREVIEW_KINDS = {'mp4', 'gif', 'none'};
    K.RAW_KINDS = {'none', 'per_frame', 'single'};
    K.MAX_GRID = 8192;           % largest fine-grid FFT (complex single N^2 = 512 MB)
    K.MIN_PUPIL_RADIUS = 96;     % pupil samples: keeps grid-induced centroid bias < 0.001 px
    K.GIF_MAX_FRAMES = 300;
    K.PREVIEW_MAX_FPS = 60.0;
    K.MAT_INFO = ['frames(:,:,k) is frame first_frame+k-1 in camera counts (DN). truth has ' ...
        'one value per frame (truth.frame is the 0-based frame index). Pixel ' ...
        'coordinates are 0-based with (0,0) at the centre of the top-left pixel: ' ...
        'MATLAB row = y + 1, column = x + 1. params and derived hold all settings.'];
end
C = K;
end

function v = mat_part_bytes(newv)
% frames per .mat file stay below this (MAT v5: < 2 GB per variable)
persistent V
if isempty(V)
    V = 2 ^ 30;
end
if nargin
    V = newv;
end
v = V;
end

% =============================================================== parameters
function spec = param_spec()
% name, default, kind (f float, i int, b bool, t three floats, s text)
spec = { ...
    'd_mm',               200.0,          'f'   % primary diameter
    'obstruction',        0.35,           'f'   % secondary / primary diameter
    'vane_mm',            4.0,            'f'   % spider vane width
    'vane_angles_deg',    [90 210 330],   't'
    'lambda_nm',          550.0,          'f'   % central wavelength
    'fnum',               10.0,           'f'   % focal ratio
    'defocus_waves',      8.0,            'f'   % W020 at the central wavelength
    'bandwidth',          0.0,            'f'   % d(lambda)/lambda, flat spectrum
    'source_lod',         0.0,            'f'   % source diameter [lambda/D]
    'jitter_lod',         0.0,            'f'   % in-frame jitter RMS [lambda/D]
    'width',              640,            'i'   % W [px]
    'height',             512,            'i'   % H [px]
    'pixel_um',           3.45,           'f'   % pixel pitch
    'bit_depth',          12,             'i'
    'gain_e_per_dn',      2.5,            'f'
    'offset_dn',          64.0,           'f'   % black level
    'full_well_e',        10000.0,        'f'
    'read_noise_e',       2.5,            'f'   % RMS
    'dark_e_per_s',       5.0,            'f'   % per pixel
    'background_e_per_s', 0.0,            'f'   % per pixel (sky / scene)
    'shot_noise',         true,           'b'
    'rolling_shutter',    false,          'b'
    'line_time_us',       10.0,           'f'   % row-to-row delay (rolling shutter)
    'frames',             60,             'i'
    'fps',                30.0,           'f'
    'exposure_ms',        20.0,           'f'
    'start_x',            120.0,          'f'   % spot centre in frame 0 [px]
    'start_y',            300.0,          'f'
    'speed',              200.0,          'f'   % in speed_unit
    'speed_unit',         'px/s',         's'   % px/s, um/s, arcsec/s, deg/s
    'direction_deg',      20.0,           'f'   % counter-clockwise from +x on screen
    'frame_jitter_px',    0.0,            'f'   % frame-to-frame pointing RMS per axis
    'intensity_mode',     'static_peak',  's'   % 'static_peak' [% full scale] or 'e_per_s'
    'intensity',          50.0,           'f'
    'seed',               1,              'i'};
end

function spec = output_spec()
% what a run writes, and where: <folder>/<name>/
spec = { ...
    'name',    '',     's'   % '' -> video_YYYYmmdd_HHMMSS
    'folder',  '',     's'   % '' -> <project>/output
    'npy',     true,   'b'   % frames.npy: lossless DN, shape (frames, H, W)
    'tiff',    true,   'b'   % frames.tif: multi-page 8/16-bit stack
    'mat',     false,  'b'   % frames.mat: frames (H x W x N) + truth + params
    'raw',     'none', 's'   % raw binary: 'per_frame' (raw/*.raw), 'single' (frames.raw)
    'preview', 'mp4',  's'   % preview.mp4 / preview.gif / 'none'
    'truth',   true,   'b'}; % truth.csv: true spot position per frame
end

function p = VideoParams(varargin)
%VIDEOPARAMS All simulation inputs: defaults, overridden by name-value
% pairs or a struct (values are converted as when read from JSON).
p = make_struct(param_spec(), 'VideoParams', args_to_struct(varargin));
end

function o = OutputSpec(varargin)
o = make_struct(output_spec(), 'OutputSpec', args_to_struct(varargin));
C = consts();
if ~any(strcmp(o.preview, C.PREVIEW_KINDS))
    value_error('preview must be one of %s', py_tuple(C.PREVIEW_KINDS));
end
if ~any(strcmp(o.raw, C.RAW_KINDS))
    value_error('raw must be one of %s', py_tuple(C.RAW_KINDS));
end
end

function data = args_to_struct(args)
if isempty(args)
    data = struct();
elseif numel(args) == 1 && isstruct(args{1})
    data = args{1};
elseif mod(numel(args), 2) == 0 && all(cellfun(@ischar, args(1:2:end)))
    data = struct();
    for i = 1:2:numel(args)
        data.(args{i}) = args{i + 1};
    end
else
    error('cassegrain_video:TypeError', '%s', ...
          'expected name-value pairs or a struct of settings');
end
end

function s = make_struct(spec, cls, data)
s = struct();
for i = 1:size(spec, 1)
    s.(spec{i, 1}) = spec{i, 2};
end
keys = fieldnames(data);
unknown = sort(setdiff(keys, spec(:, 1)));
if ~isempty(unknown)
    value_error('unknown %s key(s) %s; valid keys: %s', cls, py_list(unknown), ...
                py_list(sort(spec(:, 1))));
end
for i = 1:numel(keys)
    j = find(strcmp(spec(:, 1), keys{i}), 1);
    s.(keys{i}) = coerce(keys{i}, spec{j, 3}, data.(keys{i}));
end
end

function v = coerce(name, kind, value)
% Convert a JSON / command-line / GUI value to the type of a field.
switch kind
    case 'b'
        if ischar(value)
            t = lower(strtrim(value));
            if any(strcmp(t, {'1', 'true', 'yes', 'on'}))
                v = true;
            elseif any(strcmp(t, {'0', 'false', 'no', 'off'}))
                v = false;
            else
                value_error('%s: expected true/false, got %s', name, py_repr(value));
            end
        else
            v = py_bool(value);
        end
    case 'i'
        x = py_float(value);
        if ~(isfinite(x) && x == fix(x))
            value_error('%s: expected an integer, got %s', name, py_repr(value));
        end
        v = x;
    case 'f'
        v = py_float(value);
    case 't'
        if ischar(value)
            value = strsplit(value, ',', 'CollapseDelimiters', false);
        end
        if iscell(value)
            v = cellfun(@py_float, value);
        else
            v = arrayfun(@py_float, double(value));
        end
        v = reshape(v, 1, []);
    otherwise
        v = py_str(value);
end
end

function validate_params(p)
need(p.d_mm > 0 && p.lambda_nm > 0 && p.fnum > 0, ...
     'D, wavelength and F# must be positive');
need(0.0 <= p.obstruction && p.obstruction <= 0.95, 'obstruction must be in 0 .. 0.95');
need(p.vane_mm >= 0, 'vane width cannot be negative');
need(numel(p.vane_angles_deg) == 3, 'three vane angles are required');
need(0.0 <= p.bandwidth && p.bandwidth <= 0.5, 'bandwidth must be in 0 .. 0.5');
need(p.source_lod >= 0 && p.jitter_lod >= 0, 'source size and jitter cannot be negative');
need(8 <= p.width && p.width <= 16384 && 8 <= p.height && p.height <= 16384, ...
     'image width and height must be 8 .. 16384 px');
need(p.pixel_um > 0, 'pixel size must be positive');
need(1 <= p.bit_depth && p.bit_depth <= 16, 'bit depth must be 1 .. 16');
need(p.gain_e_per_dn > 0 && p.full_well_e > 0, 'gain and full well must be positive');
need(0 <= p.offset_dn && p.offset_dn < dn_max(p), 'offset must lie inside the ADC range');
need(py_min(py_min(p.read_noise_e, p.dark_e_per_s), p.background_e_per_s) >= 0, ...
     'noise and background levels cannot be negative');
need(p.frames >= 1, 'the sequence needs at least one frame');
need(p.fps > 0 && p.exposure_ms > 0, 'frame rate and exposure must be positive');
need(exposure_s(p) <= (1.0 + 1e-9) / p.fps, ...
     sprintf('exposure %s ms is longer than the frame period %s ms', ...
             py_fmt('%g', p.exposure_ms), py_fmt('%.4g', 1e3 / p.fps)));
need(p.speed >= 0, 'speed cannot be negative (set the direction instead)');
C = consts();
need(any(strcmp(p.speed_unit, C.SPEED_UNITS)), ...
     ['speed unit must be one of ' py_tuple(C.SPEED_UNITS)]);
need(p.frame_jitter_px >= 0, 'frame jitter cannot be negative');
need(any(strcmp(p.intensity_mode, C.INTENSITY_MODES)), ...
     ['intensity mode must be one of ' py_tuple(C.INTENSITY_MODES)]);
need(p.intensity >= 0, 'intensity cannot be negative');
if p.rolling_shutter
    need(p.line_time_us > 0, 'line time must be positive');
    readout = p.height * p.line_time_us * 1e-6;
    need(readout <= (1.0 + 1e-9) / p.fps, ...
         sprintf(['rolling-shutter readout (%s ms for %d rows) is longer than ' ...
                  'the frame period'], py_fmt('%.4g', readout * 1e3), p.height));
    [~, vy] = velocity(p);
    need(abs(vy * p.line_time_us * 1e-6) < 0.5, ...
         'the spot moves vertically faster than half a row per line time');
end
end

function need(ok, msg)
if ~ok
    error('cassegrain_video:ValueError', '%s', msg);
end
end

function d = params_to_dict(p)
% Typed copy for JSON and .mat: integers as int64, angles as a column.
spec = param_spec();
d = struct();
for i = 1:size(spec, 1)
    v = p.(spec{i, 1});
    switch spec{i, 3}
        case 'i'
            v = int64(v);
        case 't'
            v = double(v(:));
        case 'b'
            v = logical(v);
        case 'f'
            v = double(v);
    end
    d.(spec{i, 1}) = v;
end
end

function d = output_to_dict(o)
spec = output_spec();
d = struct();
for i = 1:size(spec, 1)
    d.(spec{i, 1}) = o.(spec{i, 1});
end
end

% ---- derived scales (the VideoParams properties of the Python version) ----
function v = lf_um(p)
% lambda F# [um]: the focal-plane size of lambda/D
v = p.lambda_nm * 1e-3 * p.fnum;
end

function v = focal_mm(p)
v = p.d_mm * p.fnum;
end

function v = pixel_arcsec(p)
C = consts();
v = p.pixel_um * 1e-3 / focal_mm(p) * C.RAD_TO_ARCSEC;
end

function v = speed_px_s(p)
switch p.speed_unit
    case 'px/s'
        f = 1.0;
    case 'um/s'
        f = 1.0 / p.pixel_um;
    case 'arcsec/s'
        f = 1.0 / pixel_arcsec(p);
    case 'deg/s'
        f = 3600.0 / pixel_arcsec(p);
    otherwise
        error('cassegrain_video:KeyError', '%s', py_repr(p.speed_unit));
end
v = p.speed * f;
end

function [vx, vy] = velocity(p)
% (vx, vy) [px/s] with y pointing down the rows
a = deg2rad(p.direction_deg);
s = speed_px_s(p);
vx = s * cos(a);
vy = -s * sin(a);
end

function v = exposure_s(p)
v = p.exposure_ms * 1e-3;
end

function v = blur_px(p)
v = speed_px_s(p) * exposure_s(p);
end

function v = dn_max(p)
v = 2 ^ p.bit_depth - 1;
end

function v = full_scale_e(p)
% signal that saturates a pixel: the full well or the ADC range
v = py_min(p.full_well_e, (dn_max(p) - p.offset_dn) * p.gain_e_per_dn);
end

function s = scales(p)
%SCALES The derived quantities of a parameter set in one struct.
s.lf_um = lf_um(p);
s.focal_mm = focal_mm(p);
s.pixel_arcsec = pixel_arcsec(p);
s.speed_px_s = speed_px_s(p);
[s.vx, s.vy] = velocity(p);
s.exposure_s = exposure_s(p);
s.blur_px = blur_px(p);
s.dn_max = dn_max(p);
s.full_scale_e = full_scale_e(p);
end

% ================================================================== config
function save_config(path, params, output)
if nargin < 3 || isempty(output)
    output = OutputSpec();
end
data = struct('simulation', params_to_dict(params), 'output', output_to_dict(output));
write_text(path, json_dumps(data));
end

function [params, output] = load_config(path)
txt = read_text(path);
try
    data = jsondecode(txt);
catch err
    value_error('%s', err.message);
end
if ~isstruct(data) || ~isscalar(data)
    value_error('%s: expected a JSON object', path);
end
extra = sort(setdiff(fieldnames(data), {'simulation', 'output'}));
if ~isempty(extra)
    value_error('%s: unknown section(s) %s; expected ''simulation'' and ''output''', ...
                path, py_list(extra));
end
sim = struct();
out = struct();
if isfield(data, 'simulation')
    sim = data.simulation;
end
if isfield(data, 'output')
    out = data.output;
end
params = VideoParams(sim);
output = OutputSpec(out);
end

function [params, output] = apply_overrides(params, output, pairs)
%APPLY_OVERRIDES Apply 'key=value' strings to a VideoParams / OutputSpec pair.
ps = param_spec();
os = output_spec();
for i = 1:numel(pairs)
    pair = pairs{i};
    eq = find(pair == '=', 1);
    if isempty(eq)
        value_error('--set expects key=value, got %s', py_repr(pair));
    end
    key = strtrim(pair(1:eq - 1));
    value = pair(eq + 1:end);
    j = find(strcmp(ps(:, 1), key), 1);
    if ~isempty(j)
        params.(key) = coerce(key, ps{j, 3}, value);
        continue
    end
    j = find(strcmp(os(:, 1), key), 1);
    if ~isempty(j)
        output.(key) = coerce(key, os{j, 3}, value);
    else
        value_error('unknown key %s', py_repr(key));
    end
end
end

% ===================================================================== CLI
function status = main(argv)
[args, perr] = parse_args(argv);
if strcmp(perr, 'help')
    fprintf('%s', help_text());
    status = 0;
    return
elseif ~isempty(perr)
    fprintf(2, '%s', usage_text());
    fprintf(2, 'cassegrain_video: error: %s\n', perr);
    status = 2;
    return
end
if args.selftest
    status = selftest();
    return
end
try
    if ~isempty(args.config)
        [params, out] = load_config(args.config);
    else
        params = VideoParams();
        out = OutputSpec();
    end
    [params, out] = apply_overrides(params, out, args.set);
    if ~isempty(args.name)
        out.name = args.name;
    end
    if ~isempty(args.out)
        out.folder = args.out;
    end
    if ~isempty(args.write_config)
        validate_params(params);
        save_config(args.write_config, params, out);
        fprintf('wrote %s\n', args.write_config);
        status = 0;
        return
    end
    t0 = tic;
    sim = VideoSimulator(params);
catch err
    if any(strcmp(err.identifier, {'cassegrain_video:ValueError', 'cassegrain_video:OSError'}))
        fprintf(2, 'error: %s\n', err.message);
        status = 2;
        return
    end
    rethrow(err);
end
d = sim.derived();
fprintf(['%d frames %dx%d px, pixel %s um = %.3f arcsec, lambda/D = %.2f px ' ...
         '(model ready in %.1f s)\n'], params.frames, params.width, params.height, ...
        py_fmt('%g', params.pixel_um), d.pixel_arcsec, d.lambda_over_d_px, toc(t0));
fprintf('speed %.1f px/s = %.1f arcsec/s -> %.2f px/frame, blur %.2f px per exposure\n', ...
        d.speed_px_s, d.speed_arcsec_s, d.step_px_per_frame, d.blur_px);
fprintf('signal %.4g e-/frame; peak %.0f%% of full scale static, %.0f%% moving\n', ...
        d.signal_e_per_frame, d.static_peak_percent_fs, d.moving_peak_percent_fs);
res = sim.run(out, @cli_progress);
fprintf('saved %d frames in %.1f s to %s: %s\n', res.frames_written, ...
        res.render_seconds, res.folder, strjoin(res.files, ', '));
status = 0;
end

function cli_progress(done, total)
if done == total || mod(done, max(1, floor(total / 10))) == 0
    fprintf('  frame %d/%d\n', done, total);
end
end

function [a, err] = parse_args(argv)
% argparse rules: options may be abbreviated, '--opt=value' is accepted,
% '--' ends the options, unknown arguments are an error.
a = struct('config', '', 'set', {{}}, 'name', '', 'out', '', ...
           'write_config', '', 'selftest', false);
err = '';
opts = {'--help', '--set', '--name', '--out', '--write-config', '--selftest'};
takes = [false, true, true, true, true, false];
extras = {};
npos = 0;
only_pos = false;
i = 1;
while i <= numel(argv)
    t = argv{i};
    if ~ischar(t)
        t = num2str(t);
    end
    if ~only_pos && strcmp(t, '--')
        only_pos = true;
    elseif ~only_pos && is_flag(t)
        if strcmp(t, '-h')
            err = 'help';
            return
        end
        eq = find(t == '=', 1);
        if strncmp(t, '--', 2) && ~isempty(eq)
            name = t(1:eq - 1);
            val = t(eq + 1:end);
            hasval = true;
        else
            name = t;
            val = '';
            hasval = false;
        end
        j = find(strcmp(opts, name));
        if isempty(j) && strncmp(name, '--', 2)
            j = find(strncmp(opts, name, numel(name)));
            if numel(j) > 1
                err = sprintf('ambiguous option: %s could match %s', name, ...
                              strjoin(opts(j), ', '));
                return
            end
        end
        if isempty(j)
            extras{end + 1} = t; %#ok<AGROW>
        elseif strcmp(opts{j}, '--help')
            err = 'help';
            return
        elseif takes(j)
            if ~hasval
                if i == numel(argv) || is_flag(argv{i + 1})
                    err = sprintf('argument %s: expected one argument', opts{j});
                    return
                end
                i = i + 1;
                val = argv{i};
            end
            switch opts{j}
                case '--set'
                    a.set{end + 1} = val;
                case '--name'
                    a.name = val;
                case '--out'
                    a.out = val;
                otherwise
                    a.write_config = val;
            end
        elseif hasval
            err = sprintf('argument %s: ignored explicit argument %s', opts{j}, py_repr(val));
            return
        else
            a.selftest = true;
        end
    else
        npos = npos + 1;
        if npos == 1
            a.config = t;
        else
            extras{end + 1} = t; %#ok<AGROW>
        end
    end
    i = i + 1;
end
if ~isempty(extras)
    err = ['unrecognized arguments: ' strjoin(extras, ' ')];
end
end

function tf = is_flag(t)
tf = ischar(t) && numel(t) > 1 && t(1) == '-' && isnan(str2double(t));
end

function t = usage_text()
t = sprintf(['usage: cassegrain_video [-h] [--set KEY=VALUE] [--name NAME] [--out OUT]\n' ...
             '                        [--write-config PATH] [--selftest]\n' ...
             '                        [config]\n']);
end

function t = help_text()
t = [usage_text() sprintf(['\n' ...
    'Cassegrain video simulator (headless). Without a config file the defaults are\n' ...
    'used.\n\n' ...
    'positional arguments:\n' ...
    '  config               JSON configuration (see --write-config)\n\n' ...
    'options:\n' ...
    '  -h, --help           show this help message and exit\n' ...
    '  --set KEY=VALUE      override a simulation or output value (repeatable)\n' ...
    '  --name NAME          run name = output sub-folder\n' ...
    '  --out OUT            parent output folder (default: <project>/output)\n' ...
    '  --write-config PATH  write the configuration (defaults + --set) to PATH and\n' ...
    '                       exit\n' ...
    '  --selftest           run headless checks\n\n' ...
    'examples: cassegrain_video run.json --set speed=500 --set direction_deg=45\n' ...
    '--name fast_pass | cassegrain_video --set mat=true --set raw=per_frame (raw:\n' ...
    'none, per_frame, single)\n'])];
end

% ============================================================== spot model
function nl = n_lambda(w020, bandwidth)
% spectral samples for a flat band (same rule as cassegrain_gui.py)
if bandwidth <= 0
    nl = 1;
else
    nl = floor(min(max(3 + 20.0 * bandwidth + 0.7 * bandwidth * abs(w020), 3), 15));
end
end

function k = even_fast_len(n)
% smallest even FFT length >= n with prime factors <= 11 (scipy next_fast_len)
if isnan(n)
    value_error('cannot convert float NaN to integer');
elseif isinf(n)
    error('cassegrain_video:OverflowError', '%s', 'cannot convert float infinity to integer');
end
k = next_fast_len(max(2, ceil(n)));
while mod(k, 2)
    k = next_fast_len(k + 1);
end
end

function best = next_fast_len(t)
% smallest 11-smooth integer >= t
best = 2 ^ ceil(log2(t));
if best / 2 >= t
    best = best / 2;
end
p11 = 1;
while p11 < best
    p7 = p11;
    while p7 < best
        p5 = p7;
        while p5 < best
            p3 = p5;
            while p3 < best
                q = p3;
                while q < t
                    q = q * 2;
                end
                if q < best
                    best = q;
                end
                p3 = p3 * 3;
            end
            p5 = p5 * 5;
        end
        p7 = p7 * 7;
    end
    p11 = p11 * 11;
end
end

function f = fftfreq(n)
% numpy.fft.fftfreq(n): sample frequencies [cycles / sample], a row vector
f = [0:ceil(n / 2) - 1, -floor(n / 2):-1] * (1.0 / n);
end

function s = sinc_f(x)
% numpy.sinc in single precision: sin(pi x) / (pi x)
y = single(pi) * x;
s = sin(y) ./ y;
s(y == 0) = 1;
end

function c = cover(d, a, b)
% Area of a unit pixel lying on the side u <= d of a straight edge (u =
% distance from the pixel centre along the edge normal, whose components
% have magnitudes a, b): the CDF of U(-a/2, a/2) + U(-b/2, b/2), exact at
% any edge angle.
aa = max(a, b);
bb = min(a, b);
h1 = 0.5 * (aa + bb);
h2 = 0.5 * (aa - bb);
ab2 = max(2.0 * aa .* bb, 1e-9);
v2 = (d + h1) .^ 2 ./ ab2;
v3 = (d + 0.5 * aa) ./ aa;
v4 = 1.0 - (h1 - d) .^ 2 ./ ab2;
c = ones(size(d), class(d));
m = d < h1;
c(m) = v4(m);
m = d <= h2;
c(m) = v3(m);
m = d < -h2;
c(m) = v2(m);
c(d <= -h1) = 0;
end

function acc = psf_fine(p, pad0, n)
% Polychromatic PSF (sum 1, centre n/2 + 1) sampled at pad0 samples per
% lambda0 F#. The pupil Y axis is flipped so that row 1 is the top of the
% image. Pupil samples transmit power equal to their exact open area
% (anti-aliased edges), so thin vanes keep their true area and the spider
% geometry sets the centroid, not the sampling grid.
persistent KEY VAL
key = [p.obstruction, p.vane_mm / p.d_mm, reshape(p.vane_angles_deg, 1, []), ...
       p.defocus_waves, p.bandwidth, round(pad0 * 1e12) / 1e12, n];
if ~isempty(KEY) && isequal(KEY, key)
    acc = VAL;
    return
end
c = single((0:n - 1) - floor(n / 2));
[X, Y] = meshgrid(c, -c);
rr2 = X .* X + Y .* Y;
r = sqrt(rr2);
pos = r > 0;
rs = r;
rs(~pos) = 1;
na = abs(X) ./ rs;                       % radial edge normal
na(~pos) = 1;
nb = abs(Y) ./ rs;
ang = deg2rad(p.vane_angles_deg);
trig = single([cos(ang(:)), sin(ang(:))]);
nl = n_lambda(p.defocus_waves, p.bandwidth);
acc = zeros(n, n, 'single');
for k = 0:nl - 1
    if nl == 1
        rel = 1.0;
    else
        rel = 1.0 - p.bandwidth / 2.0 + p.bandwidth * (k + 0.5) / nl;
    end
    R = n / (2.0 * pad0 * rel);          % pupil radius [samples], ~ 1/lambda
    T = cover(single(R) - r, na, nb);
    if p.obstruction > 0
        T = T .* (1.0 - cover(single(p.obstruction * R) - r, na, nb));
    end
    hw = R * p.vane_mm / p.d_mm;
    if hw > 0
        for v = 1:size(trig, 1)
            ca = trig(v, 1);
            sa = trig(v, 2);
            across = Y * ca - X * sa;
            strip = cover(single(hw) - across, abs(sa), abs(ca)) ...
                  - cover(single(-hw) - across, abs(sa), abs(ca));
            T = T .* (1.0 - strip .* single(X * ca + Y * sa >= 0));
        end
    end
    % amplitude = sqrt(open area): each sample transmits power equal to its
    % open area, so energy and centroid follow the true pupil shape
    T = sqrt(T);
    % achromatic mirrors: same physical defocus, W020 ~ 1/lambda
    phase = single(2.0 * pi * p.defocus_waves / rel / (R * R)) * rr2;
    field = fft2(ifftshift(T .* exp(1i * phase)));
    psf = abs(fftshift(field)) .^ 2;
    acc = acc + psf / single(sum(T(:) .* T(:), 'double'));   % equal energy per wavelength
end
% fade the outermost halo smoothly to zero (circular cosine taper over the
% outer quarter of the field): otherwise the field edge shows as a faint
% square step, and the taper also keeps the PSF point-symmetric on the even
% grid (its unpaired edge row / column falls outside the taper)
r_edge = n / 2.0 - 1.0;
fade = min(max((single(r_edge) - r) / single(0.25 * r_edge), 0), 1);
acc = acc .* (0.5 - 0.5 * cos(single(pi) * fade));
acc = acc / single(sum(acc(:), 'double'));
KEY = key;
VAL = acc;
end

function M = spot_model(p)
%SPOT_MODEL Spot image on the sensor: PSF x source x jitter x motion x pixel.
% The intensity is sampled on a fine grid of pitch delta = pixel/m, with m
% chosen so that delta <= lambda_min F#/2 (alias-free). The whole PSF field
% (>= 192 lambda/D wide, so the faint diffraction halo is kept the same for
% any motion) sits on an N x N canvas (N = m n_s) centred at N/2 + 1,
% enlarged by the blur, source and jitter extents so nothing wraps around.
% A crop is n_s x n_s sensor pixels holding the fraction of the spot energy
% collected by each pixel.
C = consts();
lf = lf_um(p);
lam_min_f = lf * (1.0 - p.bandwidth / 2.0);
m = max(1, ceil(p.pixel_um / (0.5 * lam_min_f) - 1e-9));
delta = p.pixel_um / m;
pad0 = lf / delta;                       % fine samples per lambda0 F#
% defocus phase step <= pi/2 per pupil sample: n_pupil >= 16 |W020|
n_psf = even_fast_len(py_max(16.0 * abs(p.defocus_waves), ...
                             2.0 * C.MIN_PUPIL_RADIUS) * pad0);
[vx, vy] = velocity(p);
speed = hypot(vx, vy);
kernel_px = blur_px(p) + (p.source_lod + 8.0 * p.jitter_lod) * lf / p.pixel_um;
if p.rolling_shutter
    kernel_px = kernel_px + speed * p.line_time_us * 1e-6 * n_psf / m;
end
n_s = even_fast_len(n_psf / m + 1.15 * kernel_px + 8);
N = m * n_s;
if max(N, n_psf) > C.MAX_GRID
    value_error(['the spot plus its motion blur needs a %d^2 grid (limit %d); ' ...
                 'reduce defocus, speed x exposure or bandwidth, or use larger ' ...
                 'pixels'], max(N, n_psf), C.MAX_GRID);
end
M = struct('m', m, 'delta', delta, 'pad0', pad0, 'n_psf', n_psf, 'n_s', n_s, 'N', N);

psf = psf_fine(p, pad0, n_psf);
canvas = zeros(N, N, 'single');
o = N / 2 - n_psf / 2;
canvas(o + 1:o + n_psf, o + 1:o + n_psf) = psf;
F = fft2(canvas);
clear canvas psf
f = single(fftfreq(N));                  % cycles / fine sample (row vector)
pix = sinc_f(single(m) * f);             % pixel = m x m samples
F = F .* pix(:);
F = F .* pix;
src = single(p.source_lod * pad0);       % diameter [samples]
sig = single(p.jitter_lod * pad0);       % RMS [samples]
if src > 0 || sig > 0
    for j = 1:512:N                      % bounded temporaries
        cols = j:min(j + 511, N);
        q = hypot(f(:), f(cols));
        blk = F(:, cols);
        if src > 0                       % uniform disk 2 J1(x) / x
            x = single(pi) * src * q;
            S = ones(size(x), 'single');
            nz = x > 1e-12;
            S(nz) = 2 * single(besselj(1, double(x(nz)))) ./ x(nz);
            blk = blk .* S;
        end
        if sig > 0                       % Gaussian jitter
            blk = blk .* exp(single(-2.0 * pi ^ 2) * sig * sig * q .* q);
        end
        F(:, cols) = blk;
    end
end
M.kc = fftfreq(N) * N;                   % centred frequencies
M.F = F;
clear F
crop = crop_global(M, 0.0, 0.0);
M.static_peak = double(max(crop(:)));
L = single(blur_px(p) * m);              % blur length [samples]
if L > 0                                 % box blur along v
    ux = single(vx / speed);
    uy = single(vy / speed);
    for j = 1:512:N
        cols = j:min(j + 511, N);
        M.F(:, cols) = M.F(:, cols) .* sinc_f(L * (f(cols) * ux + f(:) * uy));
    end
end
crop = crop_global(M, 0.0, 0.0);
M.moving_peak = double(max(crop(:)));
end

function [crop, r0, c0] = crop_global(M, xc, yc)
%CROP_GLOBAL Energy fraction per pixel of a spot centred at (xc, yc) [px];
% crop(1, 1) is the frame pixel at 0-based row r0, column c0.
m = M.m;
n_s = M.n_s;
N = M.N;
c0 = int_round(xc) - n_s / 2;
r0 = int_round(yc) - n_s / 2;
% shift the canvas so that every m-th fine sample is a pixel centre
ry = exp(-2i * pi * M.kc * (m * (yc - r0 - n_s / 2) / N));
rx = exp(-2i * pi * M.kc * (m * (xc - c0 - n_s / 2) / N));
G = M.F .* single(ry(:));
G = G .* single(rx);
G = reshape(sum(sum(reshape(G, n_s, m, n_s, m), 2), 4), n_s, n_s);   % fold aliases
crop = real(ifft2(G));
end

function [crop, r0, c0] = crop_rolling(M, xc, yc, vx, vy, line_s)
%CROP_ROLLING Rolling shutter: row r is exposed r*line_s later than row 0
% and sees the spot centred at (xc + vx r line_s, yc + vy r line_s).
m = M.m;
n_s = M.n_s;
N = M.N;
ax_row = vx * line_s;                    % px of motion per row
ay_row = vy * line_s;
[x_app, y_app] = rolling_apparent(xc, yc, vx, vy, line_s);
c0 = int_round(x_app) - n_s / 2;
r0 = int_round(y_app) - n_s / 2;
rows = r0 + (0:n_s - 1)';
ny = N / 2 + m * (rows - (yc + ay_row * rows));     % canvas row of each crop row
nx0 = N / 2 + m * (c0 - (xc + ax_row * rows));      % canvas column of crop column 1
Ey = single(exp(2i * pi * (ny .* M.kc) / N));
Hm = Ey * M.F;                                      % evaluate rows: n_s x N
Hm = Hm .* single(exp(2i * pi * (nx0 .* M.kc) / N));
Hm = sum(reshape(Hm, n_s, n_s, m), 3);              % fold column aliases
crop = real(ifft(Hm, [], 2)) / n_s;
end

function [x_app, y_app] = rolling_apparent(xc, yc, vx, vy, line_s)
% Centre of the spot recorded with a rolling shutter: the row r* that sees
% the spot centre at its own exposure time, r* = yc + vy line_s r*.
y_app = yc / (1.0 - vy * line_s);
x_app = xc + vx * line_s * y_app;
end

function M = get_model(p)
% SpotModel for p, reusing the previous one when only non-optical settings
% (noise, intensity, position, frame count ...) changed.
persistent KEY MODEL
key = [p.d_mm, p.obstruction, p.vane_mm, reshape(p.vane_angles_deg, 1, []), ...
       p.lambda_nm, p.fnum, p.defocus_waves, p.bandwidth, p.source_lod, ...
       p.jitter_lod, p.pixel_um, p.exposure_ms, double(p.rolling_shutter), ...
       p.line_time_us, speed_px_s(p), p.direction_deg];
if isempty(KEY) || ~isequal(KEY, key)
    KEY = [];
    MODEL = [];
    MODEL = spot_model(p);
    KEY = key;
end
M = MODEL;
end

% =============================================================== simulator
function sim = VideoSimulator(params)
%VIDEOSIMULATOR Renders frames (DN) and their ground truth for one parameter set.
% Returns a struct with the settings (p, model, vx, vy, rate, dtype) and
% the methods position(k), truth_geometry(k), signal_frame(k),
% render_frame(k), display_range(), derived(), frames() and
% run(out, progress, cancel, log). Frame indices k are 0-based.
validate_params(params);
S.p = params;
S.model = get_model(params);
[S.vx, S.vy] = velocity(params);
if strcmp(params.intensity_mode, 'e_per_s')
    S.rate = params.intensity;
else
    S.rate = params.intensity / 100.0 * full_scale_e(params) ...
             / (S.model.static_peak * exposure_s(params));
end
if params.bit_depth <= 8
    S.dtype = 'uint8';
else
    S.dtype = 'uint16';
end
sim = S;
sim.position = @(k) sim_position(S, k);
sim.truth_geometry = @(k) sim_truth_geometry(S, k);
sim.signal_frame = @(k) sim_signal_frame(S, k);
sim.render_frame = @(k) sim_render_frame(S, k);
sim.display_range = @() sim_display_range(S);
sim.derived = @() sim_derived(S);
sim.frames = @() sim_frames(S);
sim.run = @(varargin) sim_run(S, varargin{:});
end

function [x, y, jx, jy] = sim_position(S, k)
% true spot centre at the middle of frame k's exposure (row 0), including
% the frame's pointing jitter (jx, jy)
p = S.p;
jx = 0.0;
jy = 0.0;
if p.frame_jitter_px > 0
    st = rng_stream(p.seed, k, 2);
    z = st.randn([1, 2]);
    jx = p.frame_jitter_px * z(1);
    jy = p.frame_jitter_px * z(2);
end
t = k / p.fps;
x = p.start_x + S.vx * t + jx;
y = p.start_y + S.vy * t + jy;
end

function g = sim_truth_geometry(S, k)
% Ground-truth geometry of frame k (no rendering): spot centre at exposure
% start / middle / end, where the recorded spot appears (differs only with
% a rolling shutter), jitter offset and blur.
p = S.p;
[xc, yc, jx, jy] = sim_position(S, k);
xa = xc;
ya = yc;
if p.rolling_shutter
    [xa, ya] = rolling_apparent(xc, yc, S.vx, S.vy, p.line_time_us * 1e-6);
end
half = 0.5 * exposure_s(p);
g = struct('frame', k, 't_start_s', k / p.fps, 't_mid_s', k / p.fps + half, ...
           'x_mid_px', xc, 'y_mid_px', yc, ...
           'x_start_px', xc - S.vx * half, 'y_start_px', yc - S.vy * half, ...
           'x_end_px', xc + S.vx * half, 'y_end_px', yc + S.vy * half, ...
           'x_apparent_px', xa, 'y_apparent_px', ya, ...
           'jitter_dx_px', jx, 'jitter_dy_px', jy, 'blur_px', blur_px(p));
end

function [spot_e, crop, geo] = sim_signal_frame(S, k)
% expected spot e- per pixel, crop and truth geometry of frame k
p = S.p;
geo = sim_truth_geometry(S, k);
if p.rolling_shutter
    [crop, r0, c0] = crop_rolling(S.model, geo.x_mid_px, geo.y_mid_px, S.vx, S.vy, ...
                                  p.line_time_us * 1e-6);
else
    [crop, r0, c0] = crop_global(S.model, geo.x_mid_px, geo.y_mid_px);
end
sig = paste(zeros(p.height, p.width), crop, r0, c0);
spot_e = sig * (S.rate * exposure_s(p));
end

function [dn, truth] = sim_render_frame(S, k)
% Frame k in DN and its ground truth. Noise is seeded per frame, so a frame
% is identical whether rendered alone or within a run.
p = S.p;
[spot_e, crop, truth] = sim_signal_frame(S, k);
e_spot = S.rate * exposure_s(p);
mean_e = spot_e + (p.dark_e_per_s + p.background_e_per_s) * exposure_s(p);
st = rng_stream(p.seed, k, 1);
if p.shot_noise
    e = poisson_sample(st, mean_e);
else
    e = mean_e;
end
full = e >= p.full_well_e;
e = min(e, p.full_well_e);
if p.read_noise_e > 0
    e = e + p.read_noise_e * st.randn(size(e));
end
top = dn_max(p);
dn = cast(min(max(rint(e / p.gain_e_per_dn + p.offset_dn), 0), top), S.dtype);
total = sum(double(crop(:))) * e_spot;
signal = sum(spot_e(:));
truth.signal_e = signal;
if total > 0
    truth.in_frame_fraction = signal / total;
else
    truth.in_frame_fraction = 0.0;
end
truth.peak_e = max(spot_e(:));
truth.saturated_px = nnz(full | (dn >= top));
end

function [lo, hi] = sim_display_range(S)
% DN range for previews: black level .. expected moving peak
p = S.p;
bg = (p.dark_e_per_s + p.background_e_per_s) * exposure_s(p);
peak = S.model.moving_peak * S.rate * exposure_s(p);
top = py_min(peak + bg, p.full_well_e) + 3.0 * sqrt(peak + bg + p.read_noise_e ^ 2);
lo = p.offset_dn;
hi = py_max(lo + 1.0, py_min(dn_max(p), lo + top / p.gain_e_per_dn));
end

function d = sim_derived(S)
% scales and photometry of this configuration (for readouts)
p = S.p;
mdl = S.model;
s = speed_px_s(p);
tau = exposure_s(p);
[vx, vy] = velocity(p);
e_frame = S.rate * tau;
fs = full_scale_e(p);
lf = lf_um(p);
pa = pixel_arcsec(p);
d = struct();
d.focal_length_mm = focal_mm(p);
d.pixel_arcsec = pa;
d.lambda_over_d_px = lf / p.pixel_um;
d.lambda_over_d_arcsec = lf / p.pixel_um * pa;
d.nyquist_ratio = 0.5 * lf / p.pixel_um;
d.donut_diameter_px = 8.0 * abs(p.defocus_waves) * lf / p.pixel_um;
d.speed_px_s = s;
d.speed_um_s = s * p.pixel_um;
d.speed_arcsec_s = s * pa;
d.speed_deg_s = s * pa / 3600.0;
d.vx_px_s = vx;
d.vy_px_s = vy;
d.step_px_per_frame = s / p.fps;
d.blur_px = blur_px(p);
d.duty_cycle = tau * p.fps;
d.signal_e_per_s = S.rate;
d.signal_e_per_frame = e_frame;
d.static_peak_e = mdl.static_peak * e_frame;
d.static_peak_percent_fs = 100.0 * mdl.static_peak * e_frame / fs;
d.moving_peak_e = mdl.moving_peak * e_frame;
d.moving_peak_percent_fs = 100.0 * mdl.moving_peak * e_frame / fs;
d.full_scale_e = fs;
d.dn_max = int64(dn_max(p));
d.background_e_per_frame = (p.dark_e_per_s + p.background_e_per_s) * tau;
d.oversampling = int64(mdl.m);
d.fine_grid = int64(mdl.N);
d.psf_grid = int64(mdl.n_psf);
d.crop_px = int64(mdl.n_s);
d.spectral_samples = int64(n_lambda(p.defocus_waves, p.bandwidth));
if p.rolling_shutter
    d.rolling_readout_ms = p.height * p.line_time_us * 1e-3;
end
end

function [dn, truth] = sim_frames(S)
% the whole sequence: frames H x W x N (DN) and the truth struct array
p = S.p;
dn = zeros(p.height, p.width, p.frames, S.dtype);
rows = cell(1, p.frames);
for k = 0:p.frames - 1
    [dn(:, :, k + 1), rows{k + 1}] = sim_render_frame(S, k);
end
truth = [rows{:}];
end

function frame = paste(frame, crop, r0, c0)
% add crop into frame at 0-based (r0, c0), clipped to the frame
[H, W] = size(frame);
[n_r, n_c] = size(crop);
ra = max(r0, 0);
rb = min(r0 + n_r, H);
ca = max(c0, 0);
cb = min(c0 + n_c, W);
if ra < rb && ca < cb
    frame(ra + 1:rb, ca + 1:cb) = frame(ra + 1:rb, ca + 1:cb) ...
        + double(crop(ra - r0 + 1:rb - r0, ca - c0 + 1:cb - c0));
end
end

function img = to_display(frame, lo, hi, gamma)
%TO_DISPLAY DN frame -> uint8 display image (gamma stretch between lo and hi).
if nargin < 4
    gamma = 0.5;
end
v = (single(frame) - single(lo)) / single(max(hi - lo, 1e-9));
v = min(max(v, 0), 1);
if gamma == 0.5
    v = sqrt(v);
else
    v = v .^ single(gamma);
end
img = uint8(floor(single(255.0) * v + single(0.5)));
end

% ------------------------------------------------------------ random numbers
function st = rng_stream(seed, k, tag)
% Independent, reproducible stream for frame k: tag 1 = pixel noise,
% tag 2 = frame-to-frame jitter.
if seed < 0
    value_error('expected non-negative integer');
end
if exist('RandStream', 'class') == 8 || exist('RandStream', 'file') == 2
    s = RandStream('threefry4x64_20', 'Seed', mod(seed, 2 ^ 32));
    s.Substream = 2 * k + tag;
    st.rand = @(n) rand(s, n, 1);
    st.randn = @(sz) randn(s, sz);
else                                     % GNU Octave (no RandStream)
    state = [mod(seed, 2 ^ 32); k; tag];
    rand('state', state); %#ok<RAND>
    randn('state', state); %#ok<RAND>
    st.rand = @(n) rand(n, 1);
    st.randn = @(sz) randn(sz);
end
end

function k = poisson_sample(st, lam)
% Exact Poisson variates (vectorized): multiplication method below lambda
% 10, transformed rejection (PTRS, Hoermann 1993) above, as numpy does.
k = zeros(size(lam));
idx = find(lam > 0 & lam < 10);
if ~isempty(idx)
    enlam = exp(-lam(idx));
    produ = ones(size(idx));
    cnt = zeros(size(idx));
    act = (1:numel(idx))';
    while ~isempty(act)
        produ(act) = produ(act) .* st.rand(numel(act));
        go = produ(act) > enlam(act);
        cnt(act(go)) = cnt(act(go)) + 1;
        act = act(go);
    end
    k(idx) = cnt;
end
idx = find(lam >= 10);
if ~isempty(idx)
    lb = lam(idx);
    slam = sqrt(lb);
    loglam = log(lb);
    b = 0.931 + 2.53 * slam;
    a = -0.059 + 0.02483 * b;
    invalpha = 1.1239 + 1.1328 ./ (b - 3.4);
    vr = 0.9277 - 3.6224 ./ (b - 2);
    out = zeros(size(idx));
    act = (1:numel(idx))';
    while ~isempty(act)
        U = st.rand(numel(act)) - 0.5;
        V = st.rand(numel(act));
        us = 0.5 - abs(U);
        kk = floor((2 * a(act) ./ us + b(act)) .* U + lb(act) + 0.43);
        ok = (us >= 0.07) & (V <= vr(act));
        bad = ~ok & ((kk < 0) | ((us < 0.013) & (V > us)));
        test = ~ok & ~bad;
        if any(test)
            c = act(test);
            ok(test) = (log(V(test)) + log(invalpha(c)) - log(a(c) ./ (us(test) .^ 2) + b(c))) ...
                       <= (-lb(c) + kk(test) .* loglam(c) - gammaln(kk(test) + 1));
        end
        out(act(ok)) = kk(ok);
        act = act(~ok);
    end
    k(idx) = out;
end
end

% ============================================================== whole runs
function summary = sim_run(S, out, progress, cancel, logf)
% Render the sequence into <out folder>/<name>/ and return a summary.
% progress(done, total) is called after each frame; the sequence stops
% after the current frame once cancel() returns true; logf(message)
% reports a changed output (default: print).
if nargin < 2 || isempty(out)
    out = OutputSpec();
end
if nargin < 3
    progress = [];
end
if nargin < 4
    cancel = [];
end
if nargin < 5 || isempty(logf)
    logf = @(msg) fprintf('%s\n', msg);
end
p = S.p;
folder = unique_dir(run_dir(out));
[ok, msg] = mkdir(folder);
if ~ok
    error('cassegrain_video:OSError', '%s', sprintf('%s: %s', msg, folder));
end
[lo, hi] = sim_display_range(S);
writers = {};
if out.npy
    writers{end + 1} = npy_writer(fullfile(folder, 'frames.npy'), p, S.dtype);
end
if out.tiff
    writers{end + 1} = tiff_writer(fullfile(folder, 'frames.tif'));
end
preview = out.preview;
if strcmp(preview, 'mp4')
    w = mp4_writer(fullfile(folder, 'preview.mp4'), p, lo, hi);
    if isempty(w)
        logf('MPEG-4 video (VideoWriter) is not available: writing preview.gif instead');
        preview = 'gif';
    else
        w.gif = gif_writer(fullfile(folder, 'preview.gif'), p, lo, hi);
        w.log = logf;
        writers{end + 1} = w;
    end
end
if strcmp(preview, 'gif')
    writers{end + 1} = gif_writer(fullfile(folder, 'preview.gif'), p, lo, hi);
end
if ~strcmp(out.raw, 'none')
    writers{end + 1} = raw_writer(folder, out.raw, S.dtype, p.frames);
end
if out.mat
    writers{end + 1} = mat_writer(folder, p, S.dtype);
end
rows = {};
t0 = tic;
cancelled = false;
try
    for k = 0:p.frames - 1
        if ~isempty(cancel) && cancel()
            cancelled = true;
            break
        end
        [dn, info] = sim_render_frame(S, k);
        for i = 1:numel(writers)
            writers{i} = writer_add(writers{i}, dn);
        end
        rows{end + 1} = info; %#ok<AGROW>
        if ~isempty(progress)
            progress(k + 1, p.frames);
        end
    end
catch err
    for i = 1:numel(writers)
        writer_close(writers{i});
    end
    rethrow(err);
end
for i = 1:numel(writers)
    writers{i} = writer_close(writers{i});
end
seconds = toc(t0);
meta = struct('simulation', params_to_dict(p), 'output', output_to_dict(out), ...
              'derived', sim_derived(S), 'created', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
if out.truth
    write_truth(fullfile(folder, 'truth.csv'), rows);
end
if ~strcmp(out.raw, 'none')
    write_raw_info(folder, p, out.raw, numel(rows), S.dtype);
end
for i = 1:numel(writers)
    if strcmp(writers{i}.kind, 'mat')
        mat_finish(writers{i}, rows, meta);
    end
end
listing = dir(folder);
names = {listing.name};
is_dir = [listing.isdir];
keep = ~ismember(names, {'.', '..'});
names = names(keep);
is_dir = is_dir(keep);
[names, order] = sort(names);
is_dir = is_dir(order);
files = cell(1, numel(names) + 1);
for i = 1:numel(names)
    if is_dir(i)
        sub = dir(fullfile(folder, names{i}));
        files{i} = sprintf('%s/ (%d files)', names{i}, ...
                           sum(~ismember({sub.name}, {'.', '..'})));
    else
        files{i} = names{i};
    end
end
files{end} = 'params.json';
summary = struct('folder', folder, 'files', {files}, 'frames_written', int64(numel(rows)), ...
                 'cancelled', cancelled, 'render_seconds', round(seconds * 1000) / 1000, ...
                 'preview_display_range_dn', [lo, hi]);
keys = fieldnames(summary);
for i = 1:numel(keys)
    meta.(keys{i}) = summary.(keys{i});
end
write_text(fullfile(folder, 'params.json'), json_dumps(meta));
end

function path = run_dir(out)
C = consts();
if isempty(out.folder)
    base = C.DEFAULT_OUTPUT;
else
    base = out.folder;
end
name = out.name;
if isempty(name)
    name = ['video_' datestr(now, 'yyyymmdd_HHMMSS')];
end
path = fullfile(base, name);
end

function path = unique_dir(path)
% path, or path_2, path_3 ... if it already exists (never overwrite)
if ~(isfolder(path) || isfile(path))
    return
end
i = 2;
while isfolder(sprintf('%s_%d', path, i)) || isfile(sprintf('%s_%d', path, i))
    i = i + 1;
end
path = sprintf('%s_%d', path, i);
end

function write_truth(path, rows)
fid = open_file(path, 'w');
closer = onCleanup(@() fclose(fid));
if isempty(rows)
    return
end
keys = fieldnames(rows{1})';
ints = ismember(keys, {'frame', 'saturated_px'});
fmts = repmat({'%.6f'}, 1, numel(keys));
fmts(ints) = {'%d'};
fwrite(fid, [strjoin(keys, ',') char([13 10])], 'char');
vals = zeros(numel(keys), numel(rows));
for j = 1:numel(rows)
    vals(:, j) = cellfun(@(key) double(rows{j}.(key)), keys)';
end
fwrite(fid, sprintf([strjoin(fmts, ',') '\r\n'], vals), 'char');
end

function w = npy_writer(path, p, dtype)
% frames.npy: lossless DN, shape (frames, H, W); a cancelled run keeps the
% frames written (the header is rewritten with the final count)
w = struct('kind', 'npy', 'path', path, 'dtype', dtype, 'k', 0, 'n', p.frames, ...
           'H', p.height, 'W', p.width);
w.fid = open_file(path, 'w');
fwrite(w.fid, npy_header(dtype, [p.frames, p.height, p.width]), 'uint8');
end

function bytes = npy_header(dtype, shape)
% NPY format 1.0 header exactly as numpy writes it
if strcmp(dtype, 'uint8')
    descr = '|u1';
else
    descr = '<u2';
end
hdr = sprintf('{''descr'': ''%s'', ''fortran_order'': False, ''shape'': (%d, %d, %d), }', ...
              descr, shape);
hdr = [hdr, repmat(' ', 1, 21 - numel(sprintf('%d', shape(1))))];  % room to grow the count
hlen = numel(hdr) + 1;
pad = 64 - mod(8 + 2 + hlen, 64);
total = hlen + pad;
bytes = [uint8([147, double('NUMPY'), 1, 0, mod(total, 256), floor(total / 256)]), ...
         uint8(hdr), repmat(uint8(32), 1, pad), uint8(10)];
end

function w = tiff_writer(path)
% frames.tif: multi-page, uncompressed 8/16-bit grayscale, laid out as
% PIL's AppendingTiffWriter writes it (one strip per page, each page
% 16-byte aligned behind a stub TIFF header)
w = struct('kind', 'tiff', 'path', path, 'k', 0, 'next_ptr', 0);
w.fid = open_file(path, 'w');
end

function w = mp4_writer(path, p, lo, hi)
% preview.mp4 (H.264); [] when VideoWriter cannot write MPEG-4 here
w = [];
if exist('VideoWriter', 'class') ~= 8 && exist('VideoWriter', 'file') ~= 2
    return
end
try
    C = consts();
    v = VideoWriter(path, 'MPEG-4');
    v.FrameRate = min(p.fps, C.PREVIEW_MAX_FPS);
    v.Quality = 80;
    open(v);
catch
    return
end
w = struct('kind', 'mp4', 'path', path, 'v', v, 'lo', lo, 'hi', hi, 'k', 0, ...
           'shape', [p.height + mod(p.height, 2), p.width + mod(p.width, 2)]);
end

function w = gif_writer(path, p, lo, hi)
% preview.gif: at most GIF_MAX_FRAMES frames, played at <= 50 fps
C = consts();
step = max(1, ceil(p.frames / C.GIF_MAX_FRAMES));
ms = max(20, py_round(1000.0 * step / min(p.fps, 50.0)));
w = struct('kind', 'gif', 'path', path, 'lo', lo, 'hi', hi, 'k', 0, 'step', step, ...
           'delay', floor(ms / 10) / 100, 'written', 0);
end

function w = raw_writer(folder, kind, dtype, n_frames)
% headerless little-endian binary frames: one file per frame in raw/, or
% all frames back to back in frames.raw
w = struct('kind', 'raw', 'mode', kind, 'dtype', dtype, 'k', 0, ...
           'digits', max(5, numel(sprintf('%d', max(n_frames - 1, 0)))), 'fid', -1);
if strcmp(kind, 'per_frame')
    w.dir = fullfile(folder, 'raw');
    [ok, msg] = mkdir(w.dir);
    if ~ok
        error('cassegrain_video:OSError', '%s', msg);
    end
else
    w.fid = open_file(fullfile(folder, 'frames.raw'), 'w');
end
end

function w = mat_writer(folder, p, dtype)
% Streams frames to a temporary file; mat_finish writes frames.mat in
% MATLAB layout (H x W x N) with the truth table and the settings.
w = struct('kind', 'mat', 'folder', folder, 'H', p.height, 'W', p.width, ...
           'dtype', dtype, 'k', 0, 'tmp', fullfile(folder, '.frames_mat.tmp'));
w.fid = open_file(w.tmp, 'w');
end

function w = writer_add(w, frame)
switch w.kind
    case 'npy'
        fwrite(w.fid, frame.', w.dtype);
    case 'tiff'
        [H, W] = size(frame);
        bits = 8 * (1 + isa(frame, 'uint16'));
        ifd = ftell(w.fid) + 8;
        fwrite(w.fid, [73, 73, 42, 0, 8, 0, 0, 0], 'uint8');
        if w.k > 0                                % link the previous page here
            fseek(w.fid, w.next_ptr, 'bof');
            fwrite(w.fid, ifd, 'uint32');
            fseek(w.fid, 0, 'eof');
        end
        tags = [256, 4, W; 257, 4, H; 258, 3, bits; 259, 3, 1; 262, 3, 1; ...
                273, 4, ifd + 114; 278, 4, H; 279, 4, H * W * bits / 8; 284, 3, 1];
        fwrite(w.fid, size(tags, 1), 'uint16');
        for i = 1:size(tags, 1)
            fwrite(w.fid, tags(i, 1:2), 'uint16');
            fwrite(w.fid, [1, tags(i, 3)], 'uint32');
        end
        w.next_ptr = ftell(w.fid);
        fwrite(w.fid, 0, 'uint32');
        fwrite(w.fid, frame.', class(frame));
        fwrite(w.fid, zeros(1, mod(-ftell(w.fid), 16)), 'uint8');
    case 'mp4'
        img = to_display(frame, w.lo, w.hi);
        if ~isequal(size(img), w.shape)
            img(w.shape(1), w.shape(2)) = 0;      % pad to even (yuv420p)
        end
        try
            writeVideo(w.v, repmat(img, [1, 1, 3]));
        catch err
            if w.k > 0
                rethrow(err);
            end
            % the encoder rejects this frame size: switch to the GIF preview
            close(w.v);
            delete(w.path);
            w.log(sprintf('MPEG-4 video cannot be written (%s): writing preview.gif instead', ...
                          err.message));
            w = w.gif;
            w = writer_add(w, frame);
            return
        end
    case 'gif'
        if mod(w.k, w.step) == 0
            img = to_display(frame, w.lo, w.hi);
            if w.written == 0
                imwrite(img, gray(256), w.path, 'gif', 'DelayTime', w.delay, 'LoopCount', Inf);
            else
                imwrite(img, gray(256), w.path, 'gif', 'DelayTime', w.delay, ...
                        'WriteMode', 'append');
            end
            w.written = w.written + 1;
        end
    case 'raw'
        if strcmp(w.mode, 'per_frame')
            name = sprintf(sprintf('frame_%%0%dd.raw', w.digits), w.k);
            fid = open_file(fullfile(w.dir, name), 'w');
            fwrite(fid, frame.', w.dtype);
            fclose(fid);
        else
            fwrite(w.fid, frame.', w.dtype);
        end
    case 'mat'
        fwrite(w.fid, frame, w.dtype);
end
w.k = w.k + 1;
end

function w = writer_close(w)
switch w.kind
    case 'npy'
        if w.k < w.n                         % cancelled: keep written frames
            fseek(w.fid, 0, 'bof');
            fwrite(w.fid, npy_header(w.dtype, [w.k, w.H, w.W]), 'uint8');
        end
        fclose(w.fid);
    case 'mp4'
        close(w.v);
    case {'raw', 'mat', 'tiff'}
        if w.fid >= 0
            fclose(w.fid);
            w.fid = -1;
        end
end
end

function mat_finish(w, rows, meta)
% frames.mat, split over frames_part001.mat, ... above mat_part_bytes()
% (MAT v5 holds < 2 GB per variable); variables as written by scipy
C = consts();
isz = 1 + strcmp(w.dtype, 'uint16');
per_part = max(1, floor(mat_part_bytes() / (w.H * w.W * isz)));
common = struct('info', C.MAT_INFO, 'params', meta.simulation, 'derived', meta.derived);
if ~isempty(rows)
    keys = fieldnames(rows{1});
    truth = struct();
    for i = 1:numel(keys)
        col = cellfun(@(r) double(r.(keys{i})), rows(:));
        if any(strcmp(keys{i}, {'frame', 'saturated_px'}))
            col = int64(col);
        end
        truth.(keys{i}) = col;
    end
    common.truth = truth;
end
try
    if w.k > 0
        fid = open_file(w.tmp, 'r');
        try
            starts = 0:per_part:w.k - 1;
            for i = 1:numel(starts)
                a = starts(i);
                n = min(per_part, w.k - a);
                frames = fread(fid, w.H * w.W * n, [w.dtype '=>' w.dtype]);
                if numel(starts) == 1
                    name = 'frames.mat';
                else
                    name = sprintf('frames_part%03d.mat', i);
                end
                s = common;
                s.frames = reshape(frames, w.H, w.W, n);
                s.first_frame = double(a);
                save(fullfile(w.folder, name), '-struct', 's', '-v6');
            end
        catch err
            fclose(fid);
            rethrow(err);
        end
        fclose(fid);
    end
catch err
    delete_quietly(w.tmp);
    rethrow(err);
end
delete_quietly(w.tmp);
end

function delete_quietly(path)
if isfile(path)
    try
        delete(path);
    catch
    end
end
end

function write_raw_info(folder, p, kind, n, dtype)
% raw_format.txt: layout of the raw files and how to read them
sz = 1 + strcmp(dtype, 'uint16');
bits = 8 * sz;
W = p.width;
H = p.height;
mtype = sprintf('uint%d', bits);
ptype = sprintf('<u%d', sz);
if strcmp(kind, 'per_frame')
    d = max(5, numel(sprintf('%d', max(p.frames - 1, 0))));
    fmt = sprintf('raw/frame_%%0%dd.raw', d);
    first = sprintf(fmt, 0);
    where = sprintf('files  : %s ... %s (%d files, one frame of %d bytes each)', ...
                    first, sprintf(fmt, max(n - 1, 0)), n, W * H * sz);
    matlab = sprintf(['fid = fopen(''%s''); img = fread(fid, [%d %d], ''%s=>%s'')''; ' ...
                      'fclose(fid);'], first, W, H, mtype, mtype);
    python = sprintf('img = np.fromfile(''%s'', dtype=''%s'').reshape(%d, %d)', ...
                     first, ptype, H, W);
    imagej = sprintf(['File > Import > Raw...: %d-bit Unsigned, width %d, height %d, ' ...
                      'little-endian byte order'], bits, W, H);
else
    where = sprintf('file   : frames.raw (%d frames back to back, %d bytes each, no header)', ...
                    n, W * H * sz);
    matlab = sprintf(['fid = fopen(''frames.raw''); v = fread(fid, Inf, ''%s=>%s''); ' ...
                      'fclose(fid); frames = permute(reshape(v, %d, %d, []), [2 1 3]);'], ...
                     mtype, mtype, W, H);
    python = sprintf(['frames = np.fromfile(''frames.raw'', dtype=''%s'')' ...
                      '.reshape(-1, %d, %d)'], ptype, H, W);
    imagej = sprintf(['File > Import > Raw...: %d-bit Unsigned, width %d, height %d, ' ...
                      'number of images %d, little-endian byte order'], bits, W, H, n);
end
lines = { ...
    'Raw frames written by cassegrain_video.m'
    sprintf('format : headerless binary, unsigned %d-bit integers, little-endian,', bits)
    '         row-major (each row left to right, row 0 = top of the image)'
    sprintf('size   : width %d x height %d pixels', W, H)
    where
    sprintf('values : camera counts (DN), %d-bit right-aligned (0 .. %d), black level %s DN', ...
            p.bit_depth, dn_max(p), py_fmt('%g', p.offset_dn))
    ''
    'Read with'
    ['  MATLAB : ' matlab]
    ['  Python : ' python]
    ['  ImageJ : ' imagej]
    ''
    'Ground truth per frame: truth.csv; all settings: params.json.'};
write_text(fullfile(folder, 'raw_format.txt'), [strjoin(lines', char(10)) char(10)]);
end

% ---------------------------------------------------------------- npy files
function info = npy_open(path)
%NPY_OPEN Header of a C-order uint8 / uint16 .npy file: shape, type, offset.
fid = open_file(path, 'r');
closer = onCleanup(@() fclose(fid));
magic = fread(fid, 8, 'uint8=>double')';
if numel(magic) < 8 || ~isequal(magic(1:6), [147, double('NUMPY')])
    value_error('%s: not a .npy file', path);
end
if magic(7) == 1
    hlen = fread(fid, 1, 'uint16=>double', 0, 'ieee-le');
    offset = 10 + hlen;
else
    hlen = fread(fid, 1, 'uint32=>double', 0, 'ieee-le');
    offset = 12 + hlen;
end
hdr = fread(fid, hlen, 'uint8=>char')';
descr = regexp(hdr, '''descr'':\s*''([^'']*)''', 'tokens', 'once');
shape = regexp(hdr, '''shape'':\s*\(([^)]*)\)', 'tokens', 'once');
if isempty(descr) || isempty(shape) || ~isempty(strfind(hdr, 'True')) %#ok<STREMP>
    value_error('%s: unsupported .npy layout', path);
end
switch descr{1}
    case {'|u1', '<u1'}
        dtype = 'uint8';
    case '<u2'
        dtype = 'uint16';
    otherwise
        value_error('%s: unsupported .npy type %s', path, descr{1});
end
dims = str2double(regexp(shape{1}, '\d+', 'match'));
info = struct('path', path, 'dtype', dtype, 'shape', dims, 'offset', offset, ...
              'bytes', 1 + strcmp(dtype, 'uint16'));
end

function img = npy_frame(info, k)
%NPY_FRAME Frame k (0-based) of a (frames, H, W) .npy file as an H x W array.
H = info.shape(2);
W = info.shape(3);
fid = open_file(info.path, 'r');
closer = onCleanup(@() fclose(fid));
fseek(fid, info.offset + k * H * W * info.bytes, 'bof');
img = fread(fid, [W, H], [info.dtype '=>' info.dtype], 0, 'ieee-le')';
end

function arr = read_npy(path)
%READ_NPY Whole (frames, H, W) .npy file as an H x W x frames array.
info = npy_open(path);
fid = open_file(path, 'r');
closer = onCleanup(@() fclose(fid));
fseek(fid, info.offset, 'bof');
n = prod(info.shape);
v = fread(fid, n, [info.dtype '=>' info.dtype], 0, 'ieee-le');
arr = permute(reshape(v, info.shape(3), info.shape(2), info.shape(1)), [2, 1, 3]);
end

% ------------------------------------------------------------- text files
function fid = open_file(path, mode)
fid = fopen(path, [mode 'b'], 'ieee-le');
if fid < 0
    [folder, ~] = fileparts(path);
    missing = (strcmp(mode, 'r') && ~isfile(path)) || (~isempty(folder) && ~isfolder(folder));
    if missing
        error('cassegrain_video:OSError', '%s', ...
              sprintf('[Errno 2] No such file or directory: %s', py_repr(path)));
    end
    error('cassegrain_video:OSError', '%s', ...
          sprintf('[Errno 13] Permission denied: %s', py_repr(path)));
end
end

function write_text(path, txt)
% text with the platform's line endings, as Python's Path.write_text
if ispc
    txt = strrep(txt, char(10), char([13 10]));
end
fid = open_file(path, 'w');
closer = onCleanup(@() fclose(fid));
fwrite(fid, unicode2native(txt, 'UTF-8'), 'uint8');
end

function txt = read_text(path)
fid = open_file(path, 'r');
closer = onCleanup(@() fclose(fid));
txt = native2unicode(fread(fid, [1, Inf], 'uint8=>uint8'), 'UTF-8');
end

% ==================================================== JSON (Python format)
function s = json_dumps(v, ind)
%JSON_DUMPS Text identical to Python's json.dumps(v, indent=2): integer
% classes print as integers, doubles as Python floats ('200.0').
if nargin < 2
    ind = 0;
end
nl = char(10);
inner = repmat(' ', 1, ind + 2);
if isstruct(v)
    keys = fieldnames(v);
    if isempty(keys)
        s = '{}';
        return
    end
    parts = cell(1, numel(keys));
    for i = 1:numel(keys)
        parts{i} = [inner json_str(keys{i}) ': ' json_dumps(v.(keys{i}), ind + 2)];
    end
    s = ['{' nl strjoin(parts, [',' nl]) nl repmat(' ', 1, ind) '}'];
elseif ischar(v)
    s = json_str(v);
elseif iscell(v) || numel(v) ~= 1
    if iscell(v)
        items = reshape(v, 1, []);
    else
        items = num2cell(reshape(v, 1, []));
    end
    if isempty(items)
        s = '[]';
        return
    end
    parts = cellfun(@(x) [inner json_dumps(x, ind + 2)], items, 'UniformOutput', false);
    s = ['[' nl strjoin(parts, [',' nl]) nl repmat(' ', 1, ind) ']'];
elseif islogical(v)
    if v
        s = 'true';
    else
        s = 'false';
    end
elseif isinteger(v)
    s = sprintf('%d', v);
else
    s = py_float_repr(double(v), true);
end
end

function s = json_str(t)
% JSON string with ensure_ascii escapes
c = double(t);
out = cell(1, numel(c));
for i = 1:numel(c)
    switch c(i)
        case 34
            out{i} = '\"';
        case 92
            out{i} = '\\';
        case 10
            out{i} = '\n';
        case 13
            out{i} = '\r';
        case 9
            out{i} = '\t';
        case 8
            out{i} = '\b';
        case 12
            out{i} = '\f';
        otherwise
            if c(i) < 32 || c(i) > 126
                out{i} = sprintf('\\u%04x', c(i));
            else
                out{i} = char(c(i));
            end
    end
end
s = ['"' out{:} '"'];
end

function s = py_float_repr(x, json)
% Python repr(float): the shortest digits that read back exactly
if nargin < 2
    json = false;
end
if isnan(x)
    if json
        s = 'NaN';
    else
        s = 'nan';
    end
    return
elseif isinf(x)
    if json
        s = 'Infinity';
    else
        s = 'inf';
    end
    if x < 0
        s = ['-' s];
    end
    return
elseif x == 0
    if 1 / x < 0
        s = '-0.0';
    else
        s = '0.0';
    end
    return
end
for prec = 1:17
    t = sprintf(sprintf('%%.%de', prec - 1), x);
    if str2double(t) == x
        break
    end
end
neg = t(1) == '-';
if neg
    t = t(2:end);
end
e = find(t == 'e', 1);
decpt = str2double(t(e + 1:end)) + 1;
digits = regexprep(strrep(t(1:e - 1), '.', ''), '0+$', '');
n = numel(digits);
if decpt <= -4 || decpt > 16
    s = digits(1);
    if n > 1
        s = [s '.' digits(2:end)];
    end
    if decpt - 1 < 0
        sg = '-';
    else
        sg = '+';
    end
    s = sprintf('%se%s%02d', s, sg, abs(decpt - 1));
elseif decpt <= 0
    s = ['0.' repmat('0', 1, -decpt) digits];
elseif decpt >= n
    s = [digits repmat('0', 1, decpt - n) '.0'];
else
    s = [digits(1:decpt) '.' digits(decpt + 1:end)];
end
if neg
    s = ['-' s];
end
end

% ============================================== Python-compatible helpers
function value_error(fmt, varargin)
error('cassegrain_video:ValueError', '%s', sprintf(fmt, varargin{:}));
end

function t = error_text(err)
%ERROR_TEXT 'ValueError: message', as the Python GUI shows an exception.
id = err.identifier;
if strncmp(id, 'cassegrain_video:', 17)
    name = id(18:end);
elseif isempty(id)
    name = 'Error';
else
    name = id;
end
t = sprintf('%s: %s', name, err.message);
end

function r = py_min(a, b)
% Python min(a, b): a unless b < a (NaN handling included)
if b < a
    r = b;
else
    r = a;
end
end

function r = py_max(a, b)
if b > a
    r = b;
else
    r = a;
end
end

function r = py_round(x)
% round half to even (Python round, numpy rint)
r = round(x);
tie = abs(x - fix(x)) == 0.5;
r(tie) = 2 * round(x(tie) / 2);
end

function r = rint(x)
r = py_round(x);
end

function r = int_round(x)
% Python int(round(x))
if isnan(x)
    value_error('cannot convert float NaN to integer');
elseif isinf(x)
    error('cassegrain_video:OverflowError', '%s', 'cannot convert float infinity to integer');
end
r = py_round(x);
end

function x = py_float(value)
% Python float(): numbers, or text such as ' 1e3 ', '1_000', 'inf', 'nan'
if ischar(value)
    s = strtrim(value);
    d = '\d(_?\d)*';
    pat = ['^[+-]?((' d ')(\.(' d ')?)?|\.' d ')([eE][+-]?' d ')?$'];
    if ~isempty(regexp(s, pat, 'once'))
        x = str2double(strrep(s, '_', ''));
    elseif ~isempty(regexpi(s, '^[+-]?(inf|infinity|nan)$', 'once'))
        if ~isempty(regexpi(s, 'nan', 'once'))
            x = NaN;
        elseif s(1) == '-'
            x = -Inf;
        else
            x = Inf;
        end
    else
        value_error('could not convert string to float: %s', py_repr(value));
    end
elseif (isnumeric(value) || islogical(value)) && isscalar(value)
    x = double(value);
else
    error('cassegrain_video:TypeError', '%s', sprintf( ...
          'float() argument must be a string or a real number, not ''%s''', py_type(value)));
end
end

function b = py_bool(value)
% Python bool() of a JSON value
if (isnumeric(value) || islogical(value)) && isscalar(value)
    b = value ~= 0;
elseif isstruct(value)
    b = ~isempty(fieldnames(value));
else
    b = ~isempty(value);
end
end

function s = py_str(value)
% Python str() of a JSON value
if ischar(value)
    s = reshape(value, 1, []);
elseif islogical(value) && isscalar(value)
    s = py_repr(value);
elseif isnumeric(value) && isscalar(value)
    if value == fix(value) && abs(value) < 2 ^ 53
        s = sprintf('%d', value);
    else
        s = py_float_repr(double(value));
    end
elseif isempty(value)
    s = 'None';
else
    error('cassegrain_video:TypeError', '%s', ...
          sprintf('expected text, got a %s', py_type(value)));
end
end

function s = py_repr(v)
% Python repr() of a string, bool or number (for messages)
if ischar(v)
    if any(v == '''') && ~any(v == '"')
        q = '"';
    else
        q = '''';
    end
    t = strrep(reshape(v, 1, []), '\', '\\');
    if q == ''''
        t = strrep(t, '''', '\''');
    end
    t = strrep(strrep(strrep(t, char(10), '\n'), char(13), '\r'), char(9), '\t');
    s = [q t q];
elseif islogical(v) && isscalar(v)
    if v
        s = 'True';
    else
        s = 'False';
    end
elseif isinteger(v) && isscalar(v)
    s = sprintf('%d', v);
elseif isnumeric(v) && isscalar(v)
    s = py_float_repr(double(v));
else
    s = sprintf('<%s>', py_type(v));
end
end

function t = py_type(v)
if isstruct(v)
    t = 'dict';
elseif isempty(v) && isnumeric(v)
    t = 'NoneType';
else
    t = 'list';
end
end

function s = py_list(c)
c = reshape(c, 1, []);
s = ['[' strjoin(cellfun(@py_repr, c, 'UniformOutput', false), ', ') ']'];
end

function s = py_tuple(c)
c = reshape(c, 1, []);
s = ['(' strjoin(cellfun(@py_repr, c, 'UniformOutput', false), ', ') ')'];
end

function s = py_fmt(fmt, x)
% sprintf with Python's spelling of infinities and NaN
s = strrep(strrep(sprintf(fmt, x), 'Inf', 'inf'), 'NaN', 'nan');
end

% ================================================================ selftest
function status = selftest()
%SELFTEST Headless physics and I/O checks; returns 0 if all pass.
t_start = tic;
fails = {};
quiet = struct('shot_noise', false, 'read_noise_e', 0.0, 'dark_e_per_s', 0.0, ...
               'background_e_per_s', 0.0, 'width', 256, 'height', 256);

% energy and sub-pixel centroid of a static spot (centres as drawn by
% numpy's default_rng(7).uniform(100, 156) in the Python selftest)
sim = VideoSimulator(VideoParams(with(quiet, 'speed', 0.0, 'defocus_waves', 6.0)));
centres = [135.00534612986135, 150.24397285429623; 143.43839865373084, 112.61160263947315
           116.80931195502862, 148.91899294219067; 100.29485705567218, 145.9887914294349
           144.6358880101146, 126.20435735924836; 116.96981590188156, 115.5918342776433];
errs = zeros(1, 6);
sums = zeros(1, 6);
for i = 1:6
    [crop, r0, c0] = crop_global(sim.model, centres(i, 1), centres(i, 2));
    [cx, cy] = moments(crop);
    errs(i) = hypot(cx + c0 - centres(i, 1), cy + r0 - centres(i, 2));
    sums(i) = sum(double(crop(:)));
end
check('energy conservation', max(abs(sums - 1)) < 1e-4, ...
      sprintf('sum of pixel fractions = %.6f .. %.6f', min(sums), max(sums)));
check('sub-pixel centroid', max(errs) < 0.01, sprintf('max error %.4f px', max(errs)));

% motion blur adds L^2/12 along the motion direction to the covariance
base = VideoParams(with(quiet, 'speed', 0.0, 'defocus_waves', 6.0));
sim = VideoSimulator(base);
[~, ~, cov0] = moments(crop_global(sim.model, 0.3, -0.2));
for direction = [0.0, 90.0, 30.0]
    pm = VideoParams(with(base, 'speed', 2000.0, 'direction_deg', direction));
    L = blur_px(pm);
    [vx, vy] = velocity(pm);
    u = [vx; vy] / hypot(vx, vy);
    sim = VideoSimulator(pm);
    [~, ~, cv] = moments(crop_global(sim.model, 0.3, -0.2));
    expect = L * L / 12.0 * (u * u');
    err = max(max(abs(cv - cov0 - expect))) / (L * L / 12.0);
    check(sprintf('blur covariance, dir %g deg', direction), err < 0.01, ...
          sprintf('L = %.0f px, relative error %.4f', L, err));
end

% rolling shutter: zero line time equals the global shutter; otherwise the
% recorded spot is centred where the analytic model puts it
pr = VideoParams(with(base, 'speed', 3000.0, 'direction_deg', 60.0, ...
                      'rolling_shutter', true, 'line_time_us', 30.0));
mdl = spot_model(pr);
[vx, vy] = velocity(pr);
[g, r0, c0] = crop_global(mdl, 128.4, 120.7);
[r, r1, c1] = crop_rolling(mdl, 128.4, 120.7, vx, vy, 0.0);
diff_max = max(abs(double(g(:)) - double(r(:))));
check('rolling (0 line time) == global', r0 == r1 && c0 == c1 && diff_max < 1e-6 * max(g(:)), ...
      sprintf('max diff %.2e', diff_max));
[crop, r0, c0] = crop_rolling(mdl, 128.4, 120.7, vx, vy, 30e-6);
[xa, ya] = rolling_apparent(128.4, 120.7, vx, vy, 30e-6);
[cx, cy] = moments(crop);
err = hypot(cx + c0 - xa, cy + r0 - ya);
check('rolling-shutter apparent centre', err < 0.05, ...
      sprintf('shift (%+.2f, %+.2f) px, error %.4f px', xa - 128.4, ya - 120.7, err));

% static-peak intensity calibration and ADC
ps = VideoParams(with(quiet, 'speed', 0.0, 'intensity_mode', 'static_peak', ...
                      'intensity', 40.0, 'start_x', 128.0, 'start_y', 128.0));
sim = VideoSimulator(ps);
dn = sim.render_frame(0);
want = py_round(0.40 * full_scale_e(ps) / ps.gain_e_per_dn + ps.offset_dn);
check('static-peak calibration', abs(double(max(dn(:))) - want) <= 1, ...
      sprintf('peak %d DN, expected %d DN', max(dn(:)), want));

% noise statistics on a flat field (spot outside the frame)
pn = VideoParams('start_x', -5000.0, 'speed', 0.0, 'width', 200, 'height', 200, ...
                 'background_e_per_s', 20000.0, 'dark_e_per_s', 1000.0, ...
                 'read_noise_e', 5.0, 'gain_e_per_dn', 2.0, 'exposure_ms', 20.0);
sim = VideoSimulator(pn);
dn = double(sim.render_frame(3));
mean_e = (20000.0 + 1000.0) * exposure_s(pn);
m_th = mean_e / pn.gain_e_per_dn + pn.offset_dn;
v_th = (mean_e + pn.read_noise_e ^ 2) / pn.gain_e_per_dn ^ 2 + 1.0 / 12.0;
m = mean(dn(:));
v = var(dn(:), 1);
check('noise mean / variance', abs(m - m_th) < 0.2 && abs(v / v_th - 1) < 0.03, ...
      sprintf('mean %.2f (th %.2f), var %.2f (th %.2f) DN', m, m_th, v, v_th));

% writers and per-frame reproducibility
pw = VideoParams('width', 96, 'height', 80, 'frames', 6, 'speed', 400.0, 'start_x', 30.0, ...
                 'start_y', 40.0, 'defocus_waves', 3.0, 'frame_jitter_px', 0.5);
sw = VideoSimulator(pw);
tmp = tempname();
mkdir(tmp);
remover = onCleanup(@() rmdir(tmp, 's'));
res = sw.run(OutputSpec('name', 't', 'folder', tmp, 'preview', 'gif', 'mat', true, ...
                        'raw', 'per_frame'), [], [], @noop);
run_dir = res.folder;
info = npy_open(fullfile(run_dir, 'frames.npy'));
arr = read_npy(fullfile(run_dir, 'frames.npy'));
n_tif = numel(imfinfo(fullfile(run_dir, 'frames.tif')));
tif5 = imread(fullfile(run_dir, 'frames.tif'), 5);
lines = strsplit(strtrim(read_text(fullfile(run_dir, 'truth.csv'))), char(10));
header = strsplit(strtrim(lines{1}), ',');
col = find(strcmp(header, 'x_mid_px'));
csv_x = zeros(1, numel(lines) - 1);
for i = 2:numel(lines)
    cells = strsplit(strtrim(lines{i}), ',');
    csv_x(i - 1) = str2double(cells{col});
end
meta = jsondecode(read_text(fullfile(run_dir, 'params.json')));
again = sw.render_frame(4);
ok = isequal(info.shape, [6, 80, 96]) && strcmp(info.dtype, 'uint16') && n_tif == 6 ...
     && isequal(tif5, arr(:, :, 5)) && numel(csv_x) == 6 ...
     && isfile(fullfile(run_dir, 'preview.gif')) && meta.frames_written == 6 ...
     && isequal(again, arr(:, :, 5));
check('writers + reproducibility', ok, sprintf('npy (%d, %d, %d), tif %d pages, csv %d rows, gif', ...
                                               info.shape, n_tif, numel(csv_x)));

mat = load(fullfile(run_dir, 'frames.mat'));
mx = mat.truth.x_mid_px;
ok = isequal(size(mat.frames), [80, 96, 6]) && isa(mat.frames, 'uint16') ...
     && isequal(mat.frames, arr) && max(abs(mx(:)' - csv_x)) < 1e-5 ...
     && mat.params.speed == pw.speed && ~isfile(fullfile(run_dir, '.frames_mat.tmp'));
check('frames.mat (H x W x N, truth, params)', ok, sprintf('frames (%d, %d, %d) %s', ...
                                                          size(mat.frames), class(mat.frames)));

per = zeros(80, 96, 6, 'uint16');
for k = 0:5
    fid = open_file(fullfile(run_dir, 'raw', sprintf('frame_%05d.raw', k)), 'r');
    per(:, :, k + 1) = fread(fid, [96, 80], 'uint16=>uint16', 0, 'ieee-le')';
    fclose(fid);
end
res2 = sw.run(OutputSpec('name', 's', 'folder', tmp, 'npy', false, 'tiff', false, ...
                         'preview', 'none', 'truth', false, 'raw', 'single'), [], [], @noop);
fid = open_file(fullfile(res2.folder, 'frames.raw'), 'r');
v = fread(fid, Inf, 'uint16=>uint16', 0, 'ieee-le');
fclose(fid);
one = permute(reshape(v, 96, 80, []), [2, 1, 3]);
ok = isequal(per, arr) && isequal(one, arr) && isfile(fullfile(run_dir, 'raw_format.txt'));
check('raw frames (per frame, single file)', ok, sprintf('%d files + frames.raw (%d, %d, %d)', ...
                                                         size(per, 3), size(one, 3), 80, 96));

saved = mat_part_bytes();
mat_part_bytes(2 * 80 * 96 * 2);                 % force 2 frames per .mat part
try
    res3 = sw.run(OutputSpec('name', 'm', 'folder', tmp, 'npy', false, 'tiff', false, ...
                             'preview', 'none', 'mat', true), [], [], @noop);
catch err
    mat_part_bytes(saved);
    rethrow(err);
end
mat_part_bytes(saved);
parts = dir(fullfile(res3.folder, 'frames_part*.mat'));
names = sort({parts.name});
firsts = zeros(1, numel(names));
joined = zeros(80, 96, 0, 'uint16');
for i = 1:numel(names)
    part = load(fullfile(res3.folder, names{i}));
    firsts(i) = part.first_frame;
    joined = cat(3, joined, part.frames);
end
ok = numel(names) == 3 && isequal(firsts, [0, 2, 4]) && isequal(joined, arr);
check('frames.mat split into parts', ok, sprintf('%d parts, first frames 0, 2, 4', numel(names)));
if isempty(fails)
    verdict = 'passed';
else
    verdict = ['FAILED: ' strjoin(fails, ', ')];
end
fprintf('selftest %s (%.1f s)\n', verdict, toc(t_start));
status = double(~isempty(fails));

    function check(label, passed, detail)
        if passed
            tag = 'OK ';
        else
            tag = 'FAIL';
            fails{end + 1} = label;
        end
        fprintf('  %s %-34s %s\n', tag, label, detail);
    end
end

function s = with(s, varargin)
for i = 1:2:numel(varargin)
    s.(varargin{i}) = varargin{i + 1};
end
end

function [cx, cy, cv] = moments(img)
img = double(img);
[ny, nx] = size(img);
[x, y] = meshgrid(0:nx - 1, 0:ny - 1);
s = sum(img(:));
cx = sum(img(:) .* x(:)) / s;
cy = sum(img(:) .* y(:)) / s;
dx = x(:) - cx;
dy = y(:) - cy;
cxy = sum(img(:) .* dx .* dy);
cv = [sum(img(:) .* dx .^ 2), cxy; cxy, sum(img(:) .* dy .^ 2)] / s;
end

function noop(varargin) %#ok<INUSD>
end
