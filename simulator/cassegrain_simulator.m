function cassegrain_simulator(mode)
%CASSEGRAIN_SIMULATOR Interactive defocused-image simulator.
%
%   cassegrain_simulator            opens the GUI
%   cassegrain_simulator selftest   headless compute check (no GUI)
%
% MATLAB port of cassegrain_gui.py, feature-identical:
%   * System parameters: primary diameter, obstruction ratio, vane width,
%     3 independent vane angles, wavelength bounds [um], focal ratio.
%   * Live sliders: defocus W020 (+-64 waves), source diameter [lambda/D], jitter RMS
%     [lambda/D], linear motion blur [um on sensor] + direction.
%
% Physics (scalar Fourier optics, all-reflective =>  achromatic focus):
%   image = SUM_lambda |FFT(pupil_lambda * exp(i 2 pi W020(lambda) rho^2))|^2
%           then multiplied in the frequency domain by
%           source disk 2J1(x)/x, jitter exp(-2 pi^2 sigma^2 q^2),
%           drift sinc(L q_par).
%   Sub-wavelengths are realized by rescaling the pupil on a fixed grid,
%   which puts every wavelength on the same angular pixel grid exactly;
%   W020(lambda) = W020(lambda0) * lambda0/lambda (same physical defocus).
%   Pupil sampling obeys n >= 16|W020| (phase step <= pi/2 per sample),
%   via grid tiers 512/1024 for |W| <= 32/64.
%
% Requirements: MATLAB R2019b+ (uifigure/uigridlayout). No toolboxes
% (besselj is base MATLAB; sinc is inlined).

if nargin > 0 && (strcmpi(mode, 'selftest') || strcmpi(mode, '--selftest'))
    selftest();
    return
end
buildGui();
end

% ====================== state & physics (pure) ==========================

function S = initState()
S.pad  = 2;                       % samples per central lambda F# (mono initially)
S.tierTable = [32 512; 64 1024];  % |W020| limit -> pupil samples
S.wMax = 64;
S = setParams(S, 200, 0.35, 4, [90 210 330], 550, 10);
end

function S = setParams(S, d_mm, obs, vane_mm, angles_deg, lambda_nm, fnum)
if d_mm <= 0 || lambda_nm <= 0 || fnum <= 0
    error('cassegrain:params', 'D, wavelength and F# must be positive');
end
if obs < 0 || obs > 0.95
    error('cassegrain:params', 'obstruction must be in 0 .. 0.95');
end
if vane_mm < 0
    error('cassegrain:params', 'vane width cannot be negative');
end
S.D    = d_mm * 1e-3;
S.obs  = obs;
S.wv   = vane_mm * 1e-3;
S.ang  = angles_deg(:).';
S.lam  = lambda_nm * 1e-9;
S.fnum = fnum;
S.lodArcsec  = rad2deg(S.lam / S.D) * 3600;
S.dzPerWaveMm = 8 * S.lam * S.fnum^2 * 1e3;
S.tiers   = containers.Map('KeyType', 'double', 'ValueType', 'any');
S.specKey = [];
S.specF   = [];
S.facKey  = '';
S.fac     = [];
end

function npup = nPupilFor(S, w)
npup = S.tierTable(end, 2);
for k = 1:size(S.tierTable, 1)
    if abs(w) <= S.tierTable(k, 1)
        npup = S.tierTable(k, 2);
        return
    end
end
end

function [S, t] = getTier(S, npup, pad)
if isKey(S.tiers, npup)
    t = S.tiers(npup);
    if t.pad == pad, return; end
end
n = 64 * ceil(npup * pad / 64);
if n > 3072
    error('cassegrain:sampling', ...
        'Band and defocus need a %d^2 grid (limit 3072); narrow the band or reduce defocus', n);
end
dx = S.D * pad / n;
c  = single(((0:n-1) - n/2) * dx);
[X, Y] = meshgrid(c, c);
f  = single(fftfreqv(n));
t  = struct('n', n, 'pad', pad, 'X', X, 'Y', Y, 'rr2', X.^2 + Y.^2, ...
            'fy', f(:), 'fx', f(:).');
t.q  = sqrt(t.fy.^2 + t.fx.^2);          % cycles/px, n x n
t.q2 = t.q.^2;
S.tiers(npup) = t;
end

function [mask, rho2] = maskRho2(S, t, scale)
% Aperture mask and rho^2 for a pupil scaled by lambda0/lambda_i:
% enlarging the pupil in samples keeps the ANGULAR pixel identical for
% every sub-wavelength (no image interpolation needed).
R = single(S.D / 2 * scale);
mask = (t.rr2 <= R^2) & (t.rr2 >= (single(S.obs) * R)^2);
if S.wv > 0
    hw = single(S.wv / 2 * scale);
    for a = S.ang
        ca = single(cosd(a));  sa = single(sind(a));
        along  = t.X * ca + t.Y * sa;
        across = t.Y * ca - t.X * sa;
        mask = mask & ~((along >= 0) & (abs(across) <= hw));
    end
end
mask = single(mask);
rho2 = t.rr2 * single(1 / R^2);
end

function [S, F, t] = psfSpectrum(S, w, b)
npup = nPupilFor(S, w);
key = [w, b, npup];
nl = cassegrain_spectrum.n_lambda(w, b);
S.pad = 2 / cassegrain_spectrum.short_rel(b);
[S, t] = getTier(S, npup, S.pad);
if ~isequal(S.specKey, key)
    acc = single(0);
    for k = 1:nl
        if nl == 1, rel = 1; else, rel = 1 - b/2 + b*(k - 0.5)/nl; end
        [mask, rho2] = maskRho2(S, t, 1/rel);
        P = mask .* exp(1i * single(2*pi*w/rel) * rho2);
        field = fft2(P);              % |.|^2 is shift-invariant: no shifts
        acc = acc + abs(field).^2 / sum(mask(:));  % equal energy per lambda
    end
    S.specF   = fft2(acc);            % PSF kept in FFT-native layout
    S.specKey = key;
end
F = S.specF;
end

function [S, img] = computeImage(S, w, b, srcLod, jitLod, motLod, motDeg)
[S, F, t] = psfSpectrum(S, w, b);
[S, fac] = mtfFactor(S, t, srcLod, jitLod, motLod, motDeg);
if ~isempty(fac), F = F .* fac; end
img = fftshift(real(ifft2(F)));       % single shift, at the very end
img = max(img, 0);
m = max(img(:));
if m > 0, img = img / m; end
end

function [S, fac] = mtfFactor(S, t, srcLod, jitLod, motLod, motDeg)
% Combined source/jitter/drift factor, cached across defocus changes.
key = sprintf('%.4f_%.4f_%.6f_%.1f_%d_%.17g', srcLod, jitLod, motLod, motDeg, t.n, t.pad);
if strcmp(S.facKey, key)
    fac = S.fac;
    return
end
if srcLod <= 0 && jitLod <= 0 && motLod <= 0
    fac = [];
else
    fac = ones(t.n, t.n, 'single');
    if srcLod > 0
        fac = fac .* diskFactor(single(pi * srcLod * S.pad) * t.q);
    end
    if jitLod > 0
        fac = fac .* exp(single(-2 * pi^2 * (jitLod * S.pad)^2) * t.q2);
    end
    if motLod > 0
        qpar = t.fx * single(cosd(motDeg)) + t.fy * single(sind(motDeg));
        fac = fac .* sinc1(single(motLod * S.pad) * qpar);
    end
end
S.facKey = key;
S.fac = fac;
end

function v = diskFactor(x)
% 2*J1(x)/x via lookup table: besselj on megapixel grids is slow
persistent vals stepInv nv
if isempty(vals)
    nv = 16384;  xmax = 256;               % src<=16, pad<=3072/512, |q|<=sqrt(1/2)
    xs = linspace(0, xmax, nv);
    v0 = 2 * besselj(1, xs) ./ xs;  v0(1) = 1;
    vals = single(v0);  stepInv = single((nv - 1) / xmax);
end
idx = max(1, min(nv, 1 + round(x * stepInv)));
v = vals(idx);
end

function s = sinc1(x)
% normalized sinc without the Signal Processing Toolbox
s = ones(size(x), 'single');
nz = x ~= 0;
s(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end

function half = cropHalf(S, w, srcLod, jitLod, motLod, n)
geom  = 4 * abs(w) * S.pad;
extra = S.pad * (srcLod/2 + 3*jitLod + motLod/2);
half  = floor(min(n/2 - 2, max(80, 1.5*geom + extra + 20)));
end

function f = fftfreqv(n)
% numpy.fft.fftfreq equivalent (n even or odd), in cycles/sample
f = [0:ceil(n/2)-1, -floor(n/2):-1] / n;
end

function cm = infernoMap()
anchors = [  0   0   4;  40  11  84; 101  21 110; 159  42  99; ...
           212  72  66; 245 125  21; 250 193  39; 252 255 164] / 255;
cm = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 256));
end

% ============================== GUI ====================================

function buildGui()
S = initState();
b = 0;
imP = gobjects(0);          % persistent image objects: CData updates are
imH = gobjects(0);          % far cheaper than recreating them each frame
lastRender = tic;

fig = uifigure('Name', 'Cassegrain defocus simulator - 3-vane spider', ...
               'Position', [60 60 1220 620]);
gl = uigridlayout(fig, [2 2]);
gl.RowHeight   = {'1x', 22};
gl.ColumnWidth = {300, '1x'};

% ---- left control column -------------------------------------------
lp = uipanel(gl, 'Title', 'Controls');
lp.Layout.Row = 1;  lp.Layout.Column = 1;
lg = uigridlayout(lp, [27 2]);
lg.RowHeight = [repmat({22}, 1, 9), {28}, {26}, ...
                repmat({18, 40, 16}, 1, 5), {22}];
lg.ColumnWidth = {'1.6x', '1x'};
lg.RowSpacing = 2;  lg.Padding = [6 6 6 6];

paramDefs = {'Primary diameter D [mm]', 200; 'Obstruction ratio d/D', 0.35; ...
             'Vane width [mm]', 4; 'Vane angle 1 [deg]', 90; ...
             'Vane angle 2 [deg]', 210; 'Vane angle 3 [deg]', 330; ...
             'Lower wavelength [um]', 0.55; 'Upper wavelength [um]', 0.55; ...
             'Focal ratio F#', 10};
edits = gobjects(1, size(paramDefs, 1));
for k = 1:numel(edits)
    lab = uilabel(lg, 'Text', paramDefs{k, 1});
    lab.Layout.Row = k;  lab.Layout.Column = 1;
    edits(k) = uieditfield(lg, 'numeric', 'Value', paramDefs{k, 2}, ...
                           'HorizontalAlignment', 'right');
    edits(k).Layout.Row = k;  edits(k).Layout.Column = 2;
end
btn = uibutton(lg, 'Text', 'Apply system parameters', ...
               'ButtonPushedFcn', @(~,~) applyParamsCb());
btn.Layout.Row = 10;  btn.Layout.Column = [1 2];
errLab = uilabel(lg, 'Text', '', 'FontColor', [0.8 0 0]);
errLab.Layout.Row = 11;  errLab.Layout.Column = [1 2];

sliderDefs = { ...
    'defocus', 'Defocus W020 [waves]',                  -S.wMax, S.wMax, 8; ...
    'source',  'Source diameter [lambda/D]',             0, 16,   0; ...
    'jitter',  'Jitter RMS sigma [lambda/D]',            0, 4,    0; ...
    'motion',  'Motion blur per integration [um]',       0, 300,  0; ...
    'motdir',  'Motion direction [deg]',                 0, 180,  0};
sld = struct();  rdo = struct();
row = 12;
for k = 1:size(sliderDefs, 1)
    key = sliderDefs{k, 1};
    lab = uilabel(lg, 'Text', sliderDefs{k, 2});
    lab.Layout.Row = row;  lab.Layout.Column = [1 2];
    s = uislider(lg, 'Limits', [sliderDefs{k, 3} sliderDefs{k, 4}], ...
                 'Value', sliderDefs{k, 5}, 'MajorTicks', [], ...
                 'MinorTicks', [], ...
                 'ValueChangedFcn',  @(~,~) renderCb(), ...
                 'ValueChangingFcn', @(~,~) liveCb());
    s.Interruptible = 'off';  s.BusyAction = 'cancel';
    s.Layout.Row = row + 1;  s.Layout.Column = [1 2];
    r = uilabel(lg, 'Text', '', 'FontSize', 10, 'FontColor', [0.3 0.3 0.3]);
    r.Layout.Row = row + 2;  r.Layout.Column = [1 2];
    sld.(key) = s;  rdo.(key) = r;
    row = row + 3;
end
chkLog = uicheckbox(lg, 'Text', 'Log stretch (in-focus detail)', ...
                    'ValueChangedFcn', @(~,~) renderCb());
chkLog.Layout.Row = row;  chkLog.Layout.Column = [1 2];

% ---- right display column ------------------------------------------
rp = uigridlayout(gl, [1 2]);
rp.Layout.Row = 1;  rp.Layout.Column = 2;
axPupil = uiaxes(rp);  axImg = uiaxes(rp);
colormap(axPupil, gray);
colormap(axImg, infernoMap());
status = uilabel(gl, 'Text', '');
status.Layout.Row = 2;  status.Layout.Column = [1 2];

drawPupil();
renderCb();

% ---------------- nested callbacks (share S and handles) -------------
    function applyParamsCb()
        v = arrayfun(@(e) e.Value, edits);
        try
            [lambda_nm, newBand] = cassegrain_spectrum.from_um(v(7), v(8));
            S = setParams(S, v(1), v(2), v(3), v(4:6), lambda_nm, v(9));
            b = newBand;
            errLab.Text = '';
        catch err
            errLab.Text = err.message;
            return
        end
        drawPupil();
        renderCb();
    end

    function drawPupil()
        [S, t] = getTier(S, S.tierTable(1, 2), 2);
        mask = maskRho2(S, t, 1);
        npupBase = S.tierTable(1, 2);
        half = npupBase/2 + 24;
        cc = t.n/2 + 1;
        crop = mask(cc-half:cc+half-1, cc-half:cc+half-1);
        rmm = half * (S.D / npupBase) * 1e3;
        if isempty(imP) || ~isvalid(imP)
            imP = imagesc(axPupil, [-rmm rmm], [-rmm rmm], crop);
            axis(axPupil, 'image');  axPupil.YDir = 'normal';
            xlabel(axPupil, 'pupil [mm]');
        else
            set(imP, 'CData', crop, 'XData', [-rmm rmm], 'YData', [-rmm rmm]);
            axPupil.XLim = [-rmm rmm];  axPupil.YLim = [-rmm rmm];
        end
        title(axPupil, sprintf('Pupil   (\\lambda/D = %.2f arcsec)', ...
              S.lodArcsec), 'FontSize', 10);
    end

    function liveCb()
        % throttle drag events: render at most ~8 fps while sliding
        if toc(lastRender) > 0.12
            renderCb();
        end
    end

    function renderCb()
        lastRender = tic;
        w   = round(sld.defocus.Value * 2) / 2;
        src = round(sld.source.Value * 4) / 4;
        jit = round(sld.jitter.Value, 1);
        motUm = round(sld.motion.Value);
        mdir  = round(sld.motdir.Value / 5) * 5;
        motLod = motUm / (S.lam * S.fnum * 1e6);   % um on sensor -> lambda/D

        tic;
        try
            [S, img] = computeImage(S, w, b, src, jit, motLod, mdir);
        catch err
            if ~strcmp(err.identifier, 'cassegrain:sampling'), rethrow(err); end
            if ~isempty(imH), imH.Visible = 'off'; end
            title(axImg, 'No image for these settings');
            status.Text = err.message;
            return
        end
        n = size(img, 1);
        half = cropHalf(S, w, src, jit, motLod, n);
        cc = n/2 + 1;
        crop = img(cc-half:cc+half-1, cc-half:cc+half-1);
        dt = toc;

        % decimate for display: uifigure CData upload dominates large crops
        stp = max(1, ceil(size(crop, 1) / 640));
        if stp > 1, crop = crop(1:stp:end, 1:stp:end); end

        pixUm = S.lam * S.fnum / S.pad * 1e6;
        ext = double(half) * pixUm;
        if chkLog.Value
            data = log10(double(crop) + 1e-8);  cl = [-6 0];
        else
            data = sqrt(double(crop));          cl = [0 1];
        end
        if isempty(imH) || ~isvalid(imH)
            imH = imagesc(axImg, [-ext ext], [-ext ext], data);
            axis(axImg, 'image');  axImg.YDir = 'normal';
            xlabel(axImg, 'focal plane [um]');
        else
            set(imH, 'CData', data, 'XData', [-ext ext], 'YData', [-ext ext], 'Visible', 'on');
            axImg.XLim = [-ext ext];  axImg.YLim = [-ext ext];
        end
        caxis(axImg, cl);                    % clim() needs R2022a+
        title(axImg, sprintf('Image   (donut diameter ~ %.0f um)', ...
              8 * abs(w) * S.lam * S.fnum * 1e6), 'FontSize', 10);

        rdo.defocus.Text = sprintf('%+.1f waves   (dz ~ %+.1f mm)', ...
                                   w, w * S.dzPerWaveMm);
        rdo.source.Text = sprintf('%.2f lambda/D   (= %.2f arcsec)', ...
                                  src, src * S.lodArcsec);
        rdo.jitter.Text = sprintf('%.2f lambda/D   (= %.2f arcsec RMS)', ...
                                  jit, jit * S.lodArcsec);
        rdo.motion.Text = sprintf('%d um smear   (= %.2f lambda/D)', ...
                                  motUm, motLod);
        if motUm > 0
            rdo.motdir.Text = sprintf('%d deg from +x', mdir);
        else
            rdo.motdir.Text = '(inactive)';
        end
        status.Text = sprintf(['lambda0 %g um | bandwidth %.2f%% | grid %d^2 | N_lambda %d | ' ...
            'render %.0f ms | image pixel %.2f um'], ...
            S.lam*1e6, 100*b, n, cassegrain_spectrum.n_lambda(w, b), dt*1e3, pixUm);
        drawnow limitrate
    end
end

% ============================ selftest ==================================

function selftest()
S = initState();
cases = [  0 0   0   0  0   0;   8 0   0   0  0   0; ...
          32 0   4   1  0   0; -16 0   2 0.5  0   0; ...
          16 0   0   0  8  45;  16 0.2 0   0  0   0; ...
          64 0.4 4   1  4  90];
for k = 1:size(cases, 1)
    w = cases(k,1); b = cases(k,2); src = cases(k,3);
    jit = cases(k,4); mot = cases(k,5); mdir = cases(k,6);
    tic;
    [S, img] = computeImage(S, w, b, src, jit, mot, mdir);
    dt = toc * 1e3;
    fprintf(['W020=%+6.1f  b=%4.2f  src=%4.1f  jit=%3.1f  ' ...
             'mot=%4.1f@%3.0f  Nl=%2d  grid=%4d^2  max=%.2f  ' ...
             'sum=%.3e  %5.0f ms\n'], w, b, src, jit, mot, mdir, ...
            cassegrain_spectrum.n_lambda(w, b), size(img,1), max(img(:)), sum(img(:)), dt);
    assert(abs(max(img(:)) - 1) < 1e-6, 'normalization failed');
    assert(all(img(:) >= 0), 'negative intensities');
end
S = setParams(S, 300, 0.2, 1, [0 120 240], 650, 8);
[~, img] = computeImage(S, 10, 0.1, 1, 0.2, 2, 30);
assert(abs(max(img(:)) - 1) < 1e-6);
fprintf('param change OK, image max=%.2f\n', max(img(:)));

[lam, b] = cassegrain_spectrum.from_um(3.5, 5.5);
S = setParams(S, 200, 0, 0, [90 210 330], lam, 10);
[S, img] = computeImage(S, 0, b, 0, 0, 0, 0);
pitch = S.lam * S.fnum / S.pad * 1e6;
assert(abs(pitch - 17.5) < 1e-12);
wavelengths = 3.5 + ((1:2048) - 0.5) * 2/2048;
u = pi * (1:8)' * pitch ./ (10*wavelengths);
weights = wavelengths.^-2;              % peak of a unit-flux clear-aperture PSF
airy = (2*besselj(1,u)./u).^2 * weights' / sum(weights);
cc = size(img,1)/2 + 1;
err = max(abs(double(img(cc,cc+(1:8)))' - airy));
assert(err < 0.001, 'MWIR Airy profile error');
fprintf('3.5-5.5 um: pitch %g um, Airy profile max error %.2e of peak\nselftest passed\n', pitch, err);
end
