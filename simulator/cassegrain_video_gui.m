function cassegrain_video_gui(config)
%CASSEGRAIN_VIDEO_GUI Interactive front end of the Cassegrain video simulator (cassegrain_video.m).
%
% MATLAB port of cassegrain_video_gui.py.
%
% Left : parameters in tabs (Optics, Camera, Noise, Motion, Output) with
%        live readouts (scales, speed in all units, blur, photometry);
%        buttons to run the whole sequence (saved to output/<name>/ as
%        .npy / .tif / .mat / raw binary, chosen in the Output tab), cancel
%        it, and save / load the JSON configuration that the command line
%        also uses.
% Right: the full sensor frame with the true track and the current
%        exposure streak, a zoom on the spot, an exact pixel view (every
%        camera pixel drawn as a uniform k x k block of screen pixels, 1:1
%        up to 32x, with pixel grid, scale bar and the value of the pixel
%        under the mouse), and playback (frame slider, play, display
%        stretch).
%
% Every edit re-renders the frame under the slider. Noise is seeded per
% frame, so a preview is identical to the same frame of a run.
%
%   cassegrain_video_gui               opens the GUI with the defaults
%   cassegrain_video_gui('run.json')   ... with a saved configuration
%
% Rendering runs in MATLAB's thread: during a run the window stays
% responsive between frames (Cancel, slider, edits), edits are applied when
% a field is confirmed (Enter or leaving the field).
% Requirements: MATLAB R2020a or newer (UTF-8 source), no toolboxes, and
% cassegrain_video.m in the same folder.

if nargin < 1
    config = '';
end
cv = cassegrain_video('api');
if isempty(config)
    start_params = cv.VideoParams();
    start_output = cv.OutputSpec();
else
    [start_params, start_output] = cv.load_config(config);
end
is_octave = exist('OCTAVE_VERSION', 'builtin') > 0;
has_timer = ~is_octave;                 % Octave has no timer objects

% ------------------------------------------------------------- constants
% choice fields: {stored values, labels shown}
CHOICES = struct( ...
    'speed_unit', {{cv.SPEED_UNITS, {'px/s', 'µm/s', 'arcsec/s', 'deg/s'}}}, ...
    'intensity_mode', {{cv.INTENSITY_MODES, ...
                        {'Static peak [% of full scale]', 'Photoelectrons [e⁻/s]'}}}, ...
    'preview', {{cv.PREVIEW_KINDS, {'MP4 video', 'Animated GIF', 'none'}}}, ...
    'raw', {{cv.RAW_KINDS, {'none', 'one file per frame (raw/)', 'one file, all frames'}}});
% (field, label, kind): f float, i int, b bool, t three floats, c choice, s text
TABS = { ...
    'Optics', { ...
        'd_mm', 'Primary diameter D [mm]', 'f'
        'obstruction', 'Obstruction ratio d/D', 'f'
        'vane_mm', 'Vane width [mm]', 'f'
        'vane_angles_deg', 'Vane angles [deg]', 't'
        'lambda_nm', 'Wavelength [nm]', 'f'
        'fnum', 'Focal ratio F#', 'f'
        'defocus_waves', 'Defocus W020 [waves]', 'f'
        'bandwidth', 'Bandwidth Δλ/λ (flat)', 'f'
        'source_lod', 'Source diameter [λ/D]', 'f'
        'jitter_lod', 'In-frame jitter RMS [λ/D]', 'f'}
    'Camera', { ...
        'width', 'Image width W [px]', 'i'
        'height', 'Image height H [px]', 'i'
        'pixel_um', 'Pixel size [µm]', 'f'
        'frames', 'Sequence length [frames]', 'i'
        'fps', 'Frame rate [fps]', 'f'
        'exposure_ms', 'Exposure [ms]', 'f'
        'rolling_shutter', 'Rolling shutter', 'b'
        'line_time_us', 'Line time [µs] (rolling)', 'f'}
    'Noise', { ...
        'bit_depth', 'Bit depth', 'i'
        'gain_e_per_dn', 'Gain [e⁻/DN]', 'f'
        'offset_dn', 'Black level [DN]', 'f'
        'full_well_e', 'Full well [e⁻]', 'f'
        'read_noise_e', 'Read noise [e⁻ RMS]', 'f'
        'dark_e_per_s', 'Dark current [e⁻/s/px]', 'f'
        'background_e_per_s', 'Background [e⁻/s/px]', 'f'
        'shot_noise', 'Photon shot noise', 'b'
        'seed', 'Noise seed', 'i'}
    'Motion', { ...
        'start_x', 'Start x (column) [px]', 'f'
        'start_y', 'Start y (row) [px]', 'f'
        'speed', 'Speed', 'f'
        'speed_unit', 'Speed unit', 'c'
        'direction_deg', 'Direction [deg] (0 right, 90 up)', 'f'
        'frame_jitter_px', 'Frame-to-frame jitter [px RMS]', 'f'
        'intensity_mode', 'Intensity given as', 'c'
        'intensity', 'Intensity value', 'f'}
    'Output', { ...
        'name', 'Run name (blank: time stamp)', 's'
        'folder', 'Parent folder (blank: output/)', 's'
        'npy', 'frames.npy (lossless DN)', 'b'
        'tiff', 'frames.tif (multi-page)', 'b'
        'mat', 'frames.mat (MATLAB, H×W×N)', 'b'
        'raw', 'Raw binary frames (.raw)', 'c'
        'preview', 'Preview video', 'c'
        'truth', 'truth.csv (true positions)', 'b'}};
STRETCH = {'sqrt', 'linear', 'log'};
ZOOM_STEPS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32];         % pixel-view magnifications
PIXEL_ZOOMS = [{'fit spot', '1:1'}, ...
               arrayfun(@(z) sprintf('%d×', z), ZOOM_STEPS(2:end), 'UniformOutput', false)];
GRID_MIN_ZOOM = 6;                                       % pixel grid drawn from 6× up
GRID_RGB = reshape([70, 120, 190], 1, 1, 3);             % mixed 50 % into grid lines
OUTSIDE_RGB = uint8(reshape([38, 44, 66], 1, 1, 3));     % pixel view beyond the sensor
PV_HINT = 'wheel: zoom · drag: pan · double-click: re-centre';
% matplotlib's 'gray' and 'inferno' colour tables as 8-bit RGB
GRAY_LUT = repmat(hex_bytes([ ...
    '000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f2020222324242627' ...
    '28282a2b2c2c2e2f303032333434363738383a3b3c3c3e3f40414143444546474849494b4c4d4e4f' ...
    '50515153545556575859595b5c5d5e5f60616163646566676869696b6c6d6e6f7071717374757677' ...
    '7879797b7c7d7e7f808182838385868788898a8b8c8d8e8f909192939395969798999a9b9c9d9e9f' ...
    'a0a1a2a3a3a5a6a7a8a9aaabacadaeafb0b1b2b3b3b5b6b7b8b9babbbcbdbebfc0c1c2c3c3c5c6c7' ...
    'c8c9cacbcccdcecfd0d1d2d3d3d5d6d7d8d9dadbdcdddedfe0e1e2e3e3e5e6e7e8e9eaebecedeeef' ...
    'f0f1f2f3f3f5f6f7f8f9fafbfcfdfeff']), 1, 3);
INFERNO_LUT = reshape(hex_bytes([ ...
    '00000300000400000601000701010901010b02010e02021003021204031404031605041806041b07' ...
    '051d08061f0906210a07230b07260d08280e082a0f092d10092f120a32130a34140b36160b39170b' ...
    '3b190b3e1a0b401c0c431d0c451f0c47200c4a220b4c240b4e260b50270b52290b542b0a562d0a58' ...
    '2e0a5a300a5c32095d34095f3509603709613909623b09643c09653e0966400966410967430a6845' ...
    '0a69460a69480b6a4a0b6a4b0c6b4d0c6b4f0d6c500d6c520e6c530e6d550f6d570f6d58106d5a11' ...
    '6d5b116e5d126e5f126e60136e62146e63146e65156e66156e68166e6a176e6b176e6d186e6e186e' ...
    '70196e72196d731a6d751b6d761b6d781c6d7a1c6d7b1d6c7d1d6c7e1e6c801f6b811f6b83206b85' ...
    '206a86216a88216a8922698b22698d23698e24689024689125679325679526669626669827659928' ...
    '649b28649c29639e2963a02a62a12b61a32b61a42c60a62c5fa72d5fa92e5eab2e5dac2f5cae305b' ...
    'af315bb1315ab23259b43358b53357b73456b83556ba3655bb3754bd3753be3852bf3951c13a50c2' ...
    '3b4fc43c4ec53d4dc73e4cc83e4bc93f4acb4049cc4148cd4247cf4446d04544d14643d24742d448' ...
    '41d54940d64a3fd74b3ed94d3dda4e3bdb4f3adc5039dd5238de5337df5436e05634e25733e35832' ...
    'e45a31e55b30e65c2ee65e2de75f2ce8612be9622aea6428eb6527ec6726ed6825ed6a23ee6c22ef' ...
    '6d21f06f1ff0701ef1721df2741cf2751af37719f37918f47a16f57c15f57e14f68012f68111f783' ...
    '10f7850ef8870df8880cf88a0bf98c09f98e08f99008fa9107fa9306fa9506fa9706fb9906fb9b06' ...
    'fb9d06fb9e07fba007fba208fba40afba60bfba80dfbaa0efbac10fbae12fbb014fbb116fbb318fb' ...
    'b51afbb71cfbb91efabb21fabd23fabf25fac128f9c32af9c52cf9c72ff8c931f8cb34f8cd37f7cf' ...
    '3af7d13cf6d33ff6d542f5d745f5d948f4db4bf4dc4ff3de52f3e056f3e259f2e45df2e660f1e864' ...
    'f1e968f1eb6cf1ed70f1ee74f1f079f1f27df2f381f2f485f3f689f4f78df5f891f6fa95f7fb99f9' ...
    'fc9dfafda0fcfea4']), 3, 256)';
if is_octave
    MONO = 'Courier New';
else
    MONO = 'FixedWidth';
end
CYAN = [0, 1, 1];

% ---------------------------------------------------------------- state
params = [];            % VideoParams on screen ([] while a field is invalid)
output = [];
cur = [];               % struct(k, dn, sim) on screen
last_sim = [];          % simulator of the latest preview
run_info = [];          % last run: frames (npy header), params, sim, folder
playing = false;
gen = 0;                % bumped on every parameter change
want = [];              % [frame, generation] last requested
shown = [];             % struct(img, dn, sim, geo) on screen
track_key = [];
track_xy = zeros(0, 2);
pv_map = [];            % pixel view: struct(x0, y0, k, dn, p)
pv_offset = [0, 0];     % pixel view panning [px] from the spot
pv_drag = [];           % [x, y, offset] at the start of a drag
pv_mouse = [];          % mouse position over the pixel view
canvas = [320, 320];    % pixel view size [screen px]
busy = false;           % a frame render or a run is in progress
pending = [];           % frame index waiting to be rendered
run_request = {};       % {params, output} waiting to run
cancel_flag = false;
play_timer = [];
play_clock = tic;
play_due = 0;
last_draw = tic;

% ---------------------------------------------------------------- window
fig = figure('Name', 'Cassegrain video simulator — moving spot on a real camera', ...
             'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', 'WindowStyle', 'normal', ...
             'Units', 'pixels', 'Position', [60, 60, 1440, 780], 'Visible', 'off', ...
             'HandleVisibility', 'callback', 'IntegerHandle', 'off', ...
             'Color', get(0, 'defaultUicontrolBackgroundColor'), ...
             'CloseRequestFcn', @(~, ~) close_gui(), 'DeleteFcn', @(~, ~) stop_timer());
bg = get(fig, 'Color');
[tabgroup, pages, tab_buttons] = build_tabs();
tab_rows = cell(1, size(TABS, 1));
vars = struct();
for tab_i = 1:size(TABS, 1)
    tab_rows{tab_i} = build_tab(pages{tab_i}, TABS{tab_i, 2});
end
readout = uicontrol(fig, 'Style', 'text', 'HorizontalAlignment', 'left', ...
                    'FontName', MONO, 'FontSize', 9, 'BackgroundColor', bg);
b_run = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Run sequence', ...
                  'Callback', @(~, ~) start_run());
b_cancel = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Enable', 'off', ...
                     'Callback', @(~, ~) request_cancel());
prog_box = uipanel(fig, 'Units', 'pixels', 'BorderType', 'line', ...
                   'HighlightColor', [0.6, 0.6, 0.6], 'BackgroundColor', [0.93, 0.93, 0.93]);
prog_bar = uipanel(prog_box, 'Units', 'normalized', 'Position', [0, 0, 0.001, 1], ...
                   'BorderType', 'none', 'BackgroundColor', [0.02, 0.48, 0.84], 'Visible', 'off');
b_load = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Load config…', ...
                   'Callback', @(~, ~) load_cfg());
b_save = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save config…', ...
                   'Callback', @(~, ~) save_cfg());
err_lbl = uicontrol(fig, 'Style', 'text', 'HorizontalAlignment', 'left', ...
                    'ForegroundColor', [0.85, 0, 0], 'BackgroundColor', bg);

fig_panel = uipanel(fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', 'w');
ax = axes('Parent', fig_panel, 'Units', 'normalized', 'Position', [0.06, 0.07, 0.532, 0.86], ...
          'FontSize', 7, 'YDir', 'reverse', 'DataAspectRatio', [1, 1, 1], 'Box', 'on', ...
          'Layer', 'top', 'NextPlot', 'add', 'XLim', [0, 1], 'YLim', [0, 1]);
axz = axes('Parent', fig_panel, 'Units', 'normalized', 'Position', [0.676, 0.07, 0.304, 0.86], ...
           'FontSize', 7, 'YDir', 'reverse', 'DataAspectRatio', [1, 1, 1], 'Box', 'on', ...
           'Layer', 'top', 'NextPlot', 'add', 'XLim', [0, 1], 'YLim', [0, 1]);
im = image('Parent', ax, 'CData', zeros(2, 2, 3, 'uint8'), 'XData', [0, 1], 'YData', [0, 1]);
track_line = line('Parent', ax, 'XData', [], 'YData', [], 'Color', [0.839, 0.153, 0.157], ...
                  'LineWidth', 0.7);
streak = line('Parent', ax, 'XData', [], 'YData', [], 'Color', CYAN, 'LineWidth', 1.5);
imz = image('Parent', axz, 'CData', zeros(2, 2, 3, 'uint8'), 'XData', [0, 1], 'YData', [0, 1]);
streak_z = line('Parent', axz, 'XData', [], 'YData', [], 'Color', CYAN, 'LineWidth', 1.2);
mark_z = line('Parent', axz, 'XData', [], 'YData', [], 'Color', CYAN, 'LineStyle', 'none', ...
              'Marker', '+', 'MarkerSize', 10, 'LineWidth', 1.2);
msg_z = text(0.5, 0.5, '', 'Parent', axz, 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
title(ax, ' ', 'FontSize', 9);
title(axz, ' ', 'FontSize', 9);

% exact pixel view: an axes whose data units are screen pixels (y down),
% placed directly in the figure so that it gets the mouse events
pv_edges = [uipanel(fig, 'Units', 'pixels', 'BorderType', 'none'), ...
            uipanel(fig, 'Units', 'pixels', 'BorderType', 'none'), ...
            uipanel(fig, 'Units', 'pixels', 'BorderType', 'none'), ...
            uipanel(fig, 'Units', 'pixels', 'BorderType', 'none')];
set(pv_edges, 'BackgroundColor', [0.72, 0.72, 0.72]);
pv_title = uicontrol(fig, 'Style', 'text', 'String', 'Pixel view', ...
                     'HorizontalAlignment', 'left', 'BackgroundColor', bg);
pv_ax = axes('Parent', fig, 'Units', 'pixels', 'Position', [1, 1, 320, 320], ...
             'Visible', 'off', 'YDir', 'reverse', 'XLim', [0, 320], 'YLim', [0, 320], ...
             'NextPlot', 'add');
pv_img = image('Parent', pv_ax, 'CData', uint8(reshape([30, 30, 30], 1, 1, 3)), ...
               'XData', [0.5, 319.5], 'YData', [0.5, 319.5]);
pv_streak = line('Parent', pv_ax, 'XData', [], 'YData', [], 'Color', CYAN, 'LineWidth', 2);
pv_cross_h = line('Parent', pv_ax, 'XData', [], 'YData', [], 'Color', CYAN);
pv_cross_v = line('Parent', pv_ax, 'XData', [], 'YData', [], 'Color', CYAN);
pv_back = patch('Parent', pv_ax, 'XData', [0; 0; 0; 0], 'YData', [0; 0; 0; 0], ...
                'FaceColor', 'k', 'EdgeColor', 'none');
pv_bar = line('Parent', pv_ax, 'XData', [], 'YData', [], 'Color', 'w', 'LineWidth', 3);
pv_label = text(10, 300, '', 'Parent', pv_ax, 'Color', 'w', 'FontSize', 9, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
pv_msg = text(160, 160, '', 'Parent', pv_ax, 'Color', 'w', 'FontSize', 9, ...
              'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
pv_zoom_lbl = uicontrol(fig, 'Style', 'text', 'String', 'zoom', ...
                        'HorizontalAlignment', 'left', 'BackgroundColor', bg);
pv_zoom = uicontrol(fig, 'Style', 'popupmenu', 'String', PIXEL_ZOOMS, 'Value', 1, ...
                    'Callback', @(~, ~) draw_pixels());
pv_grid = uicontrol(fig, 'Style', 'checkbox', 'Value', 1, 'BackgroundColor', bg, ...
                    'String', sprintf('pixel grid (from %d×)', GRID_MIN_ZOOM), ...
                    'Callback', @(~, ~) draw_pixels());
pv_info = uicontrol(fig, 'Style', 'text', 'HorizontalAlignment', 'left', 'BackgroundColor', bg);
pv_hint = uicontrol(fig, 'Style', 'text', 'String', PV_HINT, 'HorizontalAlignment', 'left', ...
                    'ForegroundColor', [0.4, 0.4, 0.4], 'BackgroundColor', bg);

% playback controls and status bar
b_play = uicontrol(fig, 'Style', 'pushbutton', 'String', '▶ Play', ...
                   'Callback', @(~, ~) toggle_play());
frame_lbl = uicontrol(fig, 'Style', 'text', 'String', 'frame  0', ...
                      'HorizontalAlignment', 'left', 'BackgroundColor', bg);
slider = uicontrol(fig, 'Style', 'slider', 'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(~, ~) on_slider());
try                                     % live updates while dragging (MATLAB)
    addlistener(slider, 'ContinuousValueChange', @(~, ~) on_slider_drag());
catch
end
stretch_lbl = uicontrol(fig, 'Style', 'text', 'String', 'stretch', ...
                        'HorizontalAlignment', 'right', 'BackgroundColor', bg);
stretch_pop = uicontrol(fig, 'Style', 'popupmenu', 'String', STRETCH, 'Value', 1, ...
                        'Callback', @(~, ~) redraw());
truth_chk = uicontrol(fig, 'Style', 'checkbox', 'String', 'true track', 'Value', 1, ...
                      'BackgroundColor', bg, 'Callback', @(~, ~) redraw());
status_lbl = uicontrol(fig, 'Style', 'text', 'HorizontalAlignment', 'left', 'BackgroundColor', bg);

set(fig, 'WindowButtonMotionFcn', @(~, ~) on_motion(), ...
         'WindowButtonDownFcn', @(~, ~) on_press(), ...
         'WindowButtonUpFcn', @(~, ~) on_release(), ...
         'WindowScrollWheelFcn', @(~, evt) on_wheel(evt), ...
         'SizeChangedFcn', @(~, ~) layout());

try
    set_values(start_params, start_output);
catch start_err
    delete(fig);
    rethrow(start_err);
end
scr = get(0, 'ScreenSize');
fig_pos = get(fig, 'Position');
if fig_pos(3) > 0.95 * scr(3) || fig_pos(4) > 0.9 * scr(4)
    try
        set(fig, 'WindowState', 'maximized');   % larger than the screen: start maximized
    catch
        set(fig, 'Position', [1, 41, scr(3), scr(4) - 80]);
    end
else
    set(fig, 'Position', [max(1, round((scr(3) - fig_pos(3)) / 2)), ...
                          max(41, round((scr(4) - fig_pos(4)) / 2)), fig_pos(3:4)]);
end
layout();
set(fig, 'Visible', 'on');
drawnow;
layout();
edited();

% ================================================== parameter widgets
    function [group, tab_pages, buttons] = build_tabs()
        titles = TABS(:, 1)';
        tab_pages = cell(1, numel(titles));
        buttons = [];
        if ~is_octave
            group = uitabgroup(fig, 'Units', 'pixels');
            for n = 1:numel(titles)
                tab_pages{n} = uitab(group, 'Title', titles{n}, 'BackgroundColor', bg);
            end
        else                                % GNU Octave has no uitabgroup
            group = uipanel(fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', bg);
            buttons = zeros(1, numel(titles));
            for n = 1:numel(titles)
                buttons(n) = uicontrol(group, 'Style', 'togglebutton', 'String', titles{n}, ...
                                       'Value', n == 1, 'Callback', @(~, ~) show_tab(n));
                tab_pages{n} = uipanel(group, 'Units', 'pixels', 'BackgroundColor', bg, ...
                                       'Visible', onoff(n == 1));
            end
        end
    end

    function show_tab(n)
        for m = 1:numel(pages)
            set(pages{m}, 'Visible', onoff(m == n));
            set(tab_buttons(m), 'Value', m == n);
        end
    end

    function rows = build_tab(page, specs)
        rows = cell(size(specs, 1), 1);
        for n = 1:size(specs, 1)
            key = specs{n, 1};
            kind = specs{n, 3};
            lab = uicontrol(page, 'Style', 'text', 'String', specs{n, 2}, ...
                            'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            switch kind
                case 'b'
                    w = uicontrol(page, 'Style', 'checkbox', 'BackgroundColor', bg, ...
                                  'Callback', @(~, ~) on_change());
                case 'c'
                    w = uicontrol(page, 'Style', 'popupmenu', 'String', CHOICES.(key){2}, ...
                                  'Callback', @(~, ~) on_change());
                otherwise
                    w = uicontrol(page, 'Style', 'edit', 'HorizontalAlignment', 'right', ...
                                  'BackgroundColor', 'w', 'Callback', @(~, ~) on_change());
            end
            extra = [];
            if strcmp(key, 'folder')
                extra = uicontrol(page, 'Style', 'pushbutton', 'String', '…', ...
                                  'Callback', @(~, ~) browse());
            end
            rows{n} = struct('label', lab, 'widget', w, 'extra', extra, 'kind', kind);
            vars.(key) = struct('h', w, 'kind', kind);
        end
    end

    function layout()
        if ~ishghandle(fig)
            return
        end
        fp = getpixelposition(fig);
        fw = fp(3);
        fh = fp(4);
        lw = 372;                                    % left column
        x = 6;
        top = fh - 6;
        set(status_lbl, 'Position', [8, 2, max(fw - 16, 10), 18]);
        set(tabgroup, 'Position', [x, top - 300, lw, 300]);
        layout_tabs();
        y = top - 304;
        set(readout, 'Position', [x + 2, y - 214, lw - 4, 214]);
        y = y - 218;
        set(b_run, 'Position', [x, y - 28, lw - 90, 28]);
        set(b_cancel, 'Position', [x + lw - 86, y - 28, 86, 28]);
        y = y - 36;
        set(prog_box, 'Position', [x, y - 12, lw, 12]);
        y = y - 20;
        half = floor((lw - 4) / 2);
        set(b_load, 'Position', [x, y - 28, half, 28]);
        set(b_save, 'Position', [x + half + 4, y - 28, lw - half - 4, 28]);
        y = y - 34;
        set(err_lbl, 'Position', [x, max(y - 64, 24), lw, max(min(64, y - 24), 1)]);
        xr = x + lw + 10;                            % right side
        rw = max(fw - xr - 6, 300);
        pw = min(max(round(0.3 * rw), 300), 560);
        figw = max(rw - pw - 6, 100);
        yb = 24;
        main_y = yb + 48;
        mh = max(top - main_y, 120);
        set(fig_panel, 'Position', [xr, main_y, figw, mh]);
        set(b_play, 'Position', [xr, yb + 8, 92, 28]);
        sw = max(rw - 98 - 260, 80);
        set(frame_lbl, 'Position', [xr + 100, yb + 26, 160, 16]);
        set(slider, 'Position', [xr + 100, yb + 8, sw, 18]);
        cx = xr + 100 + sw + 8;
        set(stretch_lbl, 'Position', [cx, yb + 10, 50, 18]);
        set(stretch_pop, 'Position', [cx + 54, yb + 10, 80, 22]);
        set(truth_chk, 'Position', [cx + 142, yb + 10, 110, 22]);
        layout_pixel_view(xr + figw + 6, main_y, pw, mh);
    end

    function layout_tabs()
        tp = getpixelposition(tabgroup);
        if is_octave
            bw = floor(tp(3) / numel(pages));
            for n = 1:numel(pages)
                set(tab_buttons(n), 'Position', [1 + (n - 1) * bw, tp(4) - 25, bw - 2, 24]);
                set(pages{n}, 'Position', [1, 1, tp(3) - 2, tp(4) - 28]);
            end
        end
        for n = 1:numel(pages)
            pp = getpixelposition(pages{n});
            ph = pp(4);
            pw_ = pp(3);
            rows = tab_rows{n};
            for r = 1:numel(rows)
                y = ph - 6 - r * 25;
                set(rows{r}.label, 'Position', [6, y + 2, 198, 18]);
                wx = 206;
                ww = pw_ - wx - 8;
                if ~isempty(rows{r}.extra)
                    ww = ww - 30;
                    set(rows{r}.extra, 'Position', [wx + ww + 4, y, 26, 23]);
                end
                if strcmp(rows{r}.kind, 'b')
                    set(rows{r}.widget, 'Position', [wx, y + 1, 22, 22]);
                else
                    set(rows{r}.widget, 'Position', [wx, y, max(ww, 40), 23]);
                end
            end
        end
    end

    function layout_pixel_view(x, y, w, h)
        % framed box like a labelled frame: title on the top edge, canvas,
        % zoom controls, pixel readout and hint below
        top = y + h - 9;
        set(pv_edges(1), 'Position', [x, top, w, 1]);
        set(pv_edges(2), 'Position', [x, y, w, 1]);
        set(pv_edges(3), 'Position', [x, y, 1, top - y]);
        set(pv_edges(4), 'Position', [x + w - 1, y, 1, top - y]);
        set(pv_title, 'Position', [x + 8, top - 8, min(w - 16, 300), 17]);
        set(pv_hint, 'Position', [x + 8, y + 4, w - 16, 16]);
        set(pv_info, 'Position', [x + 8, y + 22, w - 16, 18]);
        set(pv_zoom_lbl, 'Position', [x + 8, y + 46, 36, 18]);
        set(pv_zoom, 'Position', [x + 46, y + 44, 84, 22]);
        set(pv_grid, 'Position', [x + 138, y + 44, max(w - 146, 40), 22]);
        cw = max(round(w - 12), 16);
        ch = max(round(top - 12 - (y + 72)), 16);
        set(pv_ax, 'Position', [x + 6, y + 72, cw, ch], 'XLim', [0, cw], 'YLim', [0, ch]);
        if ~isequal(canvas, [cw, ch])
            canvas = [cw, ch];
            if isempty(shown)               % dark canvas until the first frame
                set(pv_img, 'XData', [0.5, cw - 0.5], 'YData', [0.5, ch - 0.5]);
            end
            draw_pixels();
        end
    end

    function set_values(p, o)
        keys = fieldnames(vars);
        for n = 1:numel(keys)
            key = keys{n};
            v = vars.(key);
            if any(strcmp(key, cv.PARAM_FIELDS))
                val = p.(key);
            else
                val = o.(key);
            end
            switch v.kind
                case 'b'
                    set(v.h, 'Value', double(logical(val)));
                case 'c'
                    idx = find(strcmp(CHOICES.(key){1}, val), 1);
                    if isempty(idx)
                        error('cassegrain_video:ValueError', '%s', ...
                              sprintf('%s: %s is not one of %s', key, val, ...
                                      strjoin(CHOICES.(key){1}, ', ')));
                    end
                    set(v.h, 'Value', idx);
                case 't'
                    set(v.h, 'String', strjoin(arrayfun(@(a) sprintf('%g', a), val, ...
                                                        'UniformOutput', false), ', '));
                case 'f'
                    set(v.h, 'String', sprintf('%g', val));
                case 'i'
                    set(v.h, 'String', sprintf('%d', val));
                otherwise
                    set(v.h, 'String', val);
            end
        end
    end

    function [p, o] = read_fields()
        sim_d = struct();
        out_d = struct();
        keys = fieldnames(vars);
        for n = 1:numel(keys)
            key = keys{n};
            v = vars.(key);
            switch v.kind
                case 'b'
                    val = logical(get(v.h, 'Value'));
                case 'c'
                    val = CHOICES.(key){1}{get(v.h, 'Value')};
                otherwise
                    val = get(v.h, 'String');
            end
            if any(strcmp(key, cv.PARAM_FIELDS))
                sim_d.(key) = val;
            else
                out_d.(key) = val;
            end
        end
        p = cv.VideoParams(sim_d);
        cv.validate(p);
        o = cv.OutputSpec(out_d);
    end

    function browse()
        path = uigetdir(cv.ROOT);
        if ischar(path) && ~isempty(path)
            set(vars.folder.h, 'String', path);
            on_change();
        end
    end

    function on_change()
        edited();
    end

    function edited()
        try
            [p, o] = read_fields();
        catch ME
            params = [];
            set(err_lbl, 'String', ME.message);
            readouts();
            return
        end
        set(err_lbl, 'String', '');
        same_sim = ~isempty(params) && isequal(cv.to_dict(p), cv.to_dict(params));
        params = p;
        output = o;
        if same_sim                         % only output settings changed
            return
        end
        gen = gen + 1;
        last = max(p.frames - 1, 0);
        k = min(slider_k(), last);
        set(slider, 'Max', max(last, 1), 'Value', k, 'Enable', onoff(last > 0), ...
            'SliderStep', min([1, 10] / max(last, 1), 1));
        set(frame_lbl, 'String', sprintf('frame  %d', k));
        readouts();
        go_to(k);
    end

% ============================================================ readouts
    function readouts()
        p = params;
        if isempty(p)
            set(readout, 'String', '(fix the parameter above)');
            return
        end
        s = cv.scales(p);
        lod = s.lf_um / p.pixel_um;
        v = s.speed_px_s;
        nyq = 0.5 * lod;
        bytes_px = 1 + (p.bit_depth > 8);
        if nyq >= 1
            sampling = 'Nyquist-sampled';
        else
            sampling = sprintf('%.2f × Nyquist', nyq);
        end
        lines = { ...
            sprintf('λ/D      %.2f px = %.3f″', lod, lod * s.pixel_arcsec)
            sprintf('pixel    %.3f″  (%s)', s.pixel_arcsec, sampling)
            sprintf('donut ⌀  %.0f px', 8 * abs(p.defocus_waves) * lod)
            sprintf('speed    %.1f px/s = %.0f µm/s', v, v * p.pixel_um)
            sprintf('         = %.2f″/s = %.4f°/s', v * s.pixel_arcsec, v * s.pixel_arcsec / 3600)
            sprintf('step     %.2f px/frame', v / p.fps)
            sprintf('blur     %.2f px per exposure (duty %.0f%%)', s.blur_px, 100 * s.exposure_s * p.fps)
            sprintf('length   %.2f s, travel %.0f px', p.frames / p.fps, v * p.frames / p.fps)
            sprintf('raw data %.0f MB per file', p.frames * p.width * p.height * bytes_px / 1e6)};
        sim = last_sim;
        if ~isempty(sim) && isequal(cv.to_dict(sim.p), cv.to_dict(p))
            d = sim.derived();
            lines = [lines; { ...
                sprintf('signal   %.3g e⁻/frame', d.signal_e_per_frame)
                sprintf('         = %.3g e⁻/s', d.signal_e_per_s)
                sprintf('peak     %.0f%% FS static → %.0f%% moving', d.static_peak_percent_fs, ...
                        d.moving_peak_percent_fs)
                sprintf('backgr.  %.1f e⁻/px/frame', d.background_e_per_frame)
                sprintf('grid     %d² (pixel/%d)', d.fine_grid, d.oversampling)}];
        else
            lines = [lines; {'signal   …'; ''; 'peak     …'; 'backgr.  …'; 'grid     …'}];
        end
        set(readout, 'String', lines);
    end

% ================================== rendering / background work
    function go_to(k)
        % Show frame k unless it is already requested for these parameters.
        if isequal(want, [k, gen])
            return
        end
        want = [k, gen];
        request_frame(k);
    end

    function request_frame(k)
        % Show frame k: read it from the last run when that run used the
        % current parameters, otherwise render it.
        r = run_info;
        if ~isempty(r) && ~isempty(r.frames) && ~isempty(params) ...
                && isequal(r.params, cv.to_dict(params)) && k < r.frames.shape(1)
            try
                display_frame(k, cv.npy_frame(r.frames, k), r.sim);
                return
            catch
            end
        end
        pending = k;
        kick();
    end

    function kick()
        % The single worker: a requested run first, then the latest requested
        % frame; requests made meanwhile (callbacks during drawnow) wait.
        if busy
            return
        end
        busy = true;
        try
            while ishghandle(fig)
                if ~isempty(run_request)
                    req = run_request;
                    run_request = {};
                    do_run(req{1}, req{2});
                elseif ~isempty(pending) && ~isempty(params)
                    k = pending;
                    pending = [];
                    do_frame(params, k);
                else
                    break
                end
            end
        catch ME
            busy = false;
            rethrow(ME);
        end
        busy = false;
    end

    function do_frame(p, k)
        try
            t0 = tic;
            sim = cv.VideoSimulator(p);
            k = min(k, p.frames - 1);
            dn = sim.render_frame(k);
            dt = toc(t0);
            last_sim = sim;
            display_frame(k, dn, sim);
            readouts();
            set(status_lbl, 'String', sprintf('frame %d rendered in %.0f ms (fine grid %d²)', ...
                                              k, dt * 1e3, sim.model.N));
        catch ME                            % report any failure in the GUI
            if ishghandle(fig)
                set(err_lbl, 'String', cv.error_text(ME));
            end
        end
    end

    function do_run(p, o)
        cancel_flag = false;
        set_busy(true);
        last_draw = tic;
        try
            sim = cv.VideoSimulator(p);
            res = sim.run(o, @run_progress, @is_cancelled, @run_log);
        catch ME
            if ishghandle(fig)
                set_busy(false);
                set(err_lbl, 'String', cv.error_text(ME));
            end
            return
        end
        run_finished(res, sim);
    end

    function run_progress(done, total)
        if ishghandle(fig)
            set_progress(done / total);
            set(status_lbl, 'String', sprintf('rendering frame %d / %d …', done, total));
            if done == total || toc(last_draw) > 0.04
                drawnow;                    % Cancel, slider and edits work between frames
                last_draw = tic;
            end
        end
    end

    function run_log(msg)
        if ishghandle(fig)
            set(status_lbl, 'String', msg);
        end
    end

    function tf = is_cancelled()
        tf = cancel_flag || ~ishghandle(fig);
    end

    function request_cancel()
        cancel_flag = true;
    end

% ===================================================== run / playback
    function start_run()
        if isempty(params)
            return
        end
        stop_play();
        release_run();
        run_request = {params, output};
        set(status_lbl, 'String', 'starting the run …');
        kick();
    end

    function set_busy(on)
        set(b_run, 'Enable', onoff(~on));
        set(b_cancel, 'Enable', onoff(on));
        if on
            set_progress(0);
        end
    end

    function set_progress(frac)
        if frac <= 0
            set(prog_bar, 'Visible', 'off');
        else
            set(prog_bar, 'Position', [0, 0, min(frac, 1), 1], 'Visible', 'on');
        end
    end

    function run_finished(res, sim)
        if ~ishghandle(fig)
            return
        end
        set_busy(false);
        folder = res.folder;
        npy = fullfile(folder, 'frames.npy');
        frames_info = [];
        if isfile(npy)
            try
                frames_info = cv.npy_open(npy);
            catch
            end
        end
        run_info = struct('frames', frames_info, 'params', cv.to_dict(sim.p), 'sim', sim, ...
                          'folder', folder);
        last_sim = sim;
        readouts();
        if res.cancelled
            what = 'cancelled after';
        else
            what = 'saved';
        end
        set(status_lbl, 'String', sprintf('%s %d frames in %.1f s → %s', what, ...
                                          res.frames_written, res.render_seconds, folder));
        k = slider_k();
        want = [k, gen];
        request_frame(k);                   % now served from the saved frames
    end

    function release_run()
        run_info = [];
    end

    function toggle_play()
        playing = ~playing;
        set_play_label();
        if ~playing
            return
        end
        play_due = 0;
        if has_timer
            if isempty(play_timer) || ~isvalid(play_timer)
                play_timer = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.01, ...
                                   'Name', 'cassegrain_video_gui playback', ...
                                   'TimerFcn', @(~, ~) timer_fire());
            end
            play_clock = tic;
            if strcmp(get(play_timer, 'Running'), 'off')
                start(play_timer);
            end
        else                                % GNU Octave: no timer objects
            while playing && ishghandle(fig) && ~isempty(params)
                delay = play_tick();
                pause(delay / 1e3);
            end
        end
    end

    function set_play_label()
        if playing
            set(b_play, 'String', '❚❚ Pause');
        else
            set(b_play, 'String', '▶ Play');
        end
    end

    function stop_play()
        playing = false;
        set_play_label();
    end

    function timer_fire()
        try
            if ~playing || ~ishghandle(fig) || isempty(params)
                stop(play_timer);           % Play again restarts the ticks
                return
            end
            if toc(play_clock) * 1e3 >= play_due
                delay = play_tick();
                play_due = toc(play_clock) * 1e3 + delay;
            end
        catch ME
            if ishghandle(fig)
                set(err_lbl, 'String', cv.error_text(ME));
            end
        end
    end

    function delay = play_tick()
        % Advance one frame once the current one is on screen, then leave the
        % event loop some slack.
        t0 = tic;
        k = slider_k();
        if ~isempty(cur) && cur.k == k
            k = mod(k + 1, params.frames);
            set_slider(k);
            go_to(k);
        end
        spent = toc(t0) * 1e3;
        interval = 1e3 / min(params.fps, 25.0);
        delay = max(10, fix(interval - spent));
    end

    function stop_timer()
        if ~isempty(play_timer)
            try
                stop(play_timer);
                delete(play_timer);
            catch
            end
            play_timer = [];
        end
    end

    function on_slider()
        k = slider_k();
        set_slider(k);                      % snap to whole frames
        go_to(k);
    end

    function on_slider_drag()
        k = slider_k();
        set(frame_lbl, 'String', sprintf('frame  %d', k));
        go_to(k);
    end

    function k = slider_k()
        k = round(get(slider, 'Value'));
    end

    function set_slider(k)
        set(slider, 'Value', k);
        set(frame_lbl, 'String', sprintf('frame  %d', k));
    end

% ============================================================= drawing
    function display_frame(k, dn, sim)
        cur = struct('k', k, 'dn', dn, 'sim', sim);
        redraw();
    end

    function img = stretch_img(dn, lo, hi)
        img = (single(dn) - single(lo)) / single(max(hi - lo, 1e-9));
        img = min(max(img, 0), 1);
        switch STRETCH{get(stretch_pop, 'Value')}
            case 'sqrt'
                img = sqrt(img);
            case 'log'
                img = min(max(1 + log10(max(img, single(1e-4))) / 4, 0), 1);
        end
    end

    function tr = track_of(sim)
        p = sim.p;
        key = [p.start_x, p.start_y, sim.vx, sim.vy, p.fps, p.frame_jitter_px, p.seed, ...
               p.frames, double(p.rolling_shutter), p.line_time_us];
        if ~isequal(key, track_key)
            xy = zeros(p.frames, 2);
            for n = 0:p.frames - 1
                g = sim.truth_geometry(n);
                xy(n + 1, :) = [g.x_apparent_px, g.y_apparent_px];
            end
            track_key = key;
            track_xy = xy;
        end
        tr = track_xy;
    end

    function redraw()
        if isempty(cur) || ~ishghandle(fig)
            return
        end
        k = cur.k;
        dn = cur.dn;
        sim = cur.sim;
        p = sim.p;
        [lo, hi] = sim.display_range();
        img = stretch_img(dn, lo, hi);
        geo = sim.truth_geometry(k);
        xa = geo.x_apparent_px;
        ya = geo.y_apparent_px;
        truth = get(truth_chk, 'Value') == 1;
        sx = [geo.x_start_px, geo.x_end_px];
        sy = [geo.y_start_px, geo.y_end_px];
        set(im, 'CData', to_rgb(img, GRAY_LUT), 'XData', [0, p.width - 1], 'YData', [0, p.height - 1]);
        if truth
            tr = track_of(sim);
            set(track_line, 'XData', tr(:, 1), 'YData', tr(:, 2));
            set(streak, 'XData', sx, 'YData', sy);
        else
            set(track_line, 'XData', [], 'YData', []);
            set(streak, 'XData', [], 'YData', []);
        end
        set(ax, 'XLim', [-0.5, p.width - 0.5], 'YLim', [-0.5, p.height - 0.5]);
        set(get(ax, 'Title'), 'String', sprintf('frame %d / %d    t = %.4f s    %d×%d px', ...
                                                k, p.frames - 1, geo.t_mid_s, p.width, p.height));
        half = fix(0.5 * spot_extent_px(cv, p));
        x0 = max(fix(xa) - half, 0);
        x1 = min(fix(xa) + half + 1, p.width);
        y0 = max(fix(ya) - half, 0);
        y1 = min(fix(ya) + half + 1, p.height);
        inside = x0 < x1 && y0 < y1;
        set(imz, 'Visible', onoff(inside));
        if inside
            set(msg_z, 'String', '');
        else
            set(msg_z, 'String', 'spot outside the frame');
        end
        if truth && inside
            set(streak_z, 'XData', sx, 'YData', sy);
            set(mark_z, 'XData', xa, 'YData', ya);
        else
            set(streak_z, 'XData', [], 'YData', []);
            set(mark_z, 'XData', [], 'YData', []);
        end
        if inside
            set(imz, 'CData', to_rgb(img(y0 + 1:y1, x0 + 1:x1), INFERNO_LUT), ...
                'XData', [x0, x1 - 1], 'YData', [y0, y1 - 1]);
            set(axz, 'XLim', [x0 - 0.5, x1 - 0.5], 'YLim', [y0 - 0.5, y1 - 0.5]);
            sub = dn(y0 + 1:y1, x0 + 1:x1);
            set(get(axz, 'Title'), 'String', sprintf('spot zoom    peak %d DN of %d', ...
                                                     max(sub(:)), 2 ^ p.bit_depth - 1));
        else
            set(get(axz, 'Title'), 'String', 'spot zoom');
        end
        shown = struct('img', img, 'dn', dn, 'sim', sim, 'geo', geo);
        draw_pixels();
        drawnow;
    end

% ==================================================== exact pixel view
    function k = pv_zoom_factor(p, cw, ch)
        % screen pixels per camera pixel (an integer, so blocks are exact)
        idx = get(pv_zoom, 'Value');
        if idx > 1
            k = ZOOM_STEPS(idx - 1);
            return
        end
        fit = floor(min(cw, ch) / max(spot_extent_px(cv, p), 1.0));   % largest zoom that fits
        k = min(max(fit, 1), ZOOM_STEPS(end));
    end

    function draw_pixels()
        % Every camera pixel becomes a uniform k x k block of screen pixels
        % (integer magnification, never resampled), centred on the spot plus
        % the panning offset; the same display stretch as the other views.
        pv_map = [];
        if isempty(shown) || ~ishghandle(fig)
            pv_readout();
            return
        end
        img = shown.img;
        geo = shown.geo;
        p = shown.sim.p;
        cw = canvas(1);
        ch = canvas(2);
        k = pv_zoom_factor(p, cw, ch);
        nx = ceil(cw / k);                          % blocks covering the canvas
        ny = ceil(ch / k);
        x0 = py_round(geo.x_apparent_px + pv_offset(1)) - floor(nx / 2);
        y0 = py_round(geo.y_apparent_px + pv_offset(2)) - floor(ny / 2);
        rgb = repmat(OUTSIDE_RGB, ny, nx);
        xs0 = max(x0, 0);
        xs1 = min(x0 + nx, p.width);
        ys0 = max(y0, 0);
        ys1 = min(y0 + ny, p.height);
        inside = xs0 < xs1 && ys0 < ys1;
        if inside
            v = uint8(floor(single(255) * img(ys0 + 1:ys1, xs0 + 1:xs1) + single(0.5)));
            rgb(ys0 - y0 + 1:ys1 - y0, xs0 - x0 + 1:xs1 - x0, :) = repmat(v, [1, 1, 3]);
        end
        big = rgb(ceil((1:ch) / k), ceil((1:cw) / k), :);
        if get(pv_grid, 'Value') && k >= GRID_MIN_ZOOM  % first row / column of each block
            big(1:k:end, :, :) = uint8(floor((double(big(1:k:end, :, :)) + GRID_RGB) / 2));
            big(:, 1:k:end, :) = uint8(floor((double(big(:, 1:k:end, :)) + GRID_RGB) / 2));
        end
        set(pv_img, 'CData', big, 'XData', [0.5, cw - 0.5], 'YData', [0.5, ch - 0.5], ...
            'Visible', 'on');
        if k == 1
            ttl = 'Pixel view  ·  1:1 (actual pixels)';
        else
            ttl = sprintf('Pixel view  ·  1 px = %d×%d screen px', k, k);
        end
        if ~strcmp(get(pv_title, 'String'), ttl)
            set(pv_title, 'String', ttl);
        end
        at = @(x, y) [(x - x0 + 0.5) * k, (y - y0 + 0.5) * k];   % sensor px -> canvas
        set(pv_msg, 'String', '');
        for h = [pv_streak, pv_cross_h, pv_cross_v]
            set(h, 'XData', [], 'YData', []);
        end
        if ~inside
            if isequal(pv_offset, [0, 0])
                msg = 'spot outside the sensor';
            else
                msg = sprintf('outside the sensor\n(double-click to re-centre)');
            end
            set(pv_msg, 'Position', [cw / 2, ch / 2, 0], 'String', msg);
        elseif get(truth_chk, 'Value')
            a = at(geo.x_start_px, geo.y_start_px);
            b = at(geo.x_end_px, geo.y_end_px);
            c = at(geo.x_apparent_px, geo.y_apparent_px);
            set(pv_streak, 'XData', [a(1), b(1)], 'YData', [a(2), b(2)]);
            set(pv_cross_h, 'XData', [c(1) - 7, c(1) + 8], 'YData', [c(2), c(2)]);
            set(pv_cross_v, 'XData', [c(1), c(1)], 'YData', [c(2) - 7, c(2) + 8]);
        end
        % scale bar: n camera pixels with their size on the sensor and the sky
        steps = [1, 2, 5, 10, 20, 50, 100, 200, 500];
        steps = steps(steps * k <= 0.4 * cw);
        n = 1;
        if ~isempty(steps)
            n = steps(end);
        end
        y_bar = ch - 10;
        s = cv.scales(p);
        set(pv_bar, 'XData', [10, 10 + n * k], 'YData', [y_bar, y_bar]);
        set(pv_label, 'Position', [10, y_bar - 5, 0], 'String', ...
            sprintf('%d px = %s µm = %s″', n, sprintf('%.4g', n * p.pixel_um), ...
                    sprintf('%.3g', n * s.pixel_arcsec)));
        ext = get(pv_label, 'Extent');
        bx1 = max(10 + ext(3), 10 + n * k) + 4;
        by0 = y_bar - 5 - abs(ext(4)) - 2;
        set(pv_back, 'XData', [6, bx1, bx1, 6], 'YData', [by0, by0, ch - 3, ch - 3]);
        pv_map = struct('x0', x0, 'y0', y0, 'k', k, 'dn', shown.dn, 'p', p);
        pv_readout();
    end

    function pv_readout()
        % value of the camera pixel under the mouse
        if isempty(pv_map) || isempty(pv_mouse)
            set(pv_info, 'String', 'hover over a pixel to read its value');
            return
        end
        m = pv_map;
        ix = m.x0 + floor(floor(pv_mouse(1)) / m.k);
        iy = m.y0 + floor(floor(pv_mouse(2)) / m.k);
        if ix >= 0 && ix < m.p.width && iy >= 0 && iy < m.p.height
            v = double(m.dn(iy + 1, ix + 1));
            txt = sprintf('pixel (%d, %d):  %d DN  ≈ %.0f e⁻ above black', ix, iy, v, ...
                          (v - m.p.offset_dn) * m.p.gain_e_per_dn);
        else
            txt = sprintf('(%d, %d): outside the sensor', ix, iy);
        end
        set(pv_info, 'String', txt);
    end

    function [x, y, over] = pv_point()
        cp = get(pv_ax, 'CurrentPoint');
        x = cp(1, 1);
        y = cp(1, 2);
        over = x >= 0 && x < canvas(1) && y >= 0 && y < canvas(2);
    end

    function on_motion()
        [x, y, over] = pv_point();
        if ~isempty(pv_drag)
            pv_mouse = [x, y];
            if isempty(pv_map)
                return
            end
            k = pv_map.k;
            offset = [pv_drag(3) - py_round((floor(x) - pv_drag(1)) / k), ...
                      pv_drag(4) - py_round((floor(y) - pv_drag(2)) / k)];
            if ~isequal(offset, pv_offset)      % whole pixels: blocks stay aligned
                pv_offset = offset;
                draw_pixels();
            else
                pv_readout();
            end
        elseif over
            pv_mouse = [x, y];
            set(fig, 'Pointer', 'crosshair');
            pv_readout();
        elseif ~isempty(pv_mouse)
            pv_mouse = [];
            set(fig, 'Pointer', 'arrow');
            pv_readout();
        end
    end

    function on_press()
        [x, y, over] = pv_point();
        if ~over
            return
        end
        switch get(fig, 'SelectionType')
            case 'normal'
                pv_drag = [floor(x), floor(y), pv_offset];
            case 'open'                         % double-click
                pv_recentre();
        end
    end

    function on_release()
        pv_drag = [];
    end

    function on_wheel(evt)
        [~, ~, over] = pv_point();
        if ~over || isempty(pv_map)
            return
        end
        k = pv_map.k;
        if evt.VerticalScrollCount < 0
            steps = ZOOM_STEPS(ZOOM_STEPS > k);
            pick = 1;
        else
            steps = ZOOM_STEPS(ZOOM_STEPS < k);
            pick = numel(steps);
        end
        if ~isempty(steps)
            set(pv_zoom, 'Value', find(ZOOM_STEPS == steps(pick)) + 1);
            draw_pixels();
        end
    end

    function pv_recentre()
        pv_offset = [0, 0];
        draw_pixels();
    end

% ================================================ configuration files
    function save_cfg()
        try
            [p, o] = read_fields();
        catch ME
            set(err_lbl, 'String', ME.message);
            return
        end
        [f, d] = uiputfile({'*.json', 'JSON'}, 'Save configuration', [cv.ROOT, filesep]);
        if ischar(f)
            path = fullfile(d, f);
            try
                cv.save_config(path, p, o);
            catch ME
                set(err_lbl, 'String', sprintf('%s: %s', f, ME.message));
                return
            end
            set(status_lbl, 'String', ['saved ' path]);
        end
    end

    function load_cfg()
        [f, d] = uigetfile({'*.json', 'JSON'}, 'Load configuration', [cv.ROOT, filesep]);
        if ~ischar(f)
            return
        end
        path = fullfile(d, f);
        try
            [p, o] = cv.load_config(path);
            set_values(p, o);
        catch ME
            set(err_lbl, 'String', sprintf('%s: %s', f, ME.message));
            return
        end
        edited();
        set(status_lbl, 'String', ['loaded ' path]);
    end

    function close_gui()
        cancel_flag = true;
        playing = false;
        stop_timer();
        release_run();
        delete(fig);
    end
end

% ================================================================ helpers
function b = hex_bytes(h)
b = uint8(sscanf(h, '%2x'));
end

function s = onoff(tf)
if tf
    s = 'on';
else
    s = 'off';
end
end

function r = py_round(x)
% round half to even, as Python's round()
r = round(x);
tie = abs(x - fix(x)) == 0.5;
r(tie) = 2 * round(x(tie) / 2);
end

function rgb = to_rgb(v, lut)
% values 0..1 -> RGB through a 256-entry colour table (as matplotlib maps them)
idx = min(floor(v * 256), 255) + 1;
rgb = reshape(lut(idx(:), :), [size(v, 1), size(v, 2), 3]);
end

function e = spot_extent_px(cv, p)
% Side [px] of a box holding the spot: donut + motion blur + source, plus a
% margin of 3 lambda/D + 8 px on each side.
s = cv.scales(p);
lod = s.lf_um / p.pixel_um;
e = 8 * abs(p.defocus_waves) * lod + s.blur_px + p.source_lod * lod + 6 * lod + 16;
end