%% processSpectra.m
% Background-corrects a hyperspectral dark-field scan and extracts
% single-particle scattering spectra for fitSpectraLorentzian.m.
%
% Inputs (.mat files made by convertTDMStoMAT.m):
%   1. White count  - lamp/white-light reference (system response)
%   2. Dark count   - detector dark/offset reference
%   3. Sample       - hyperspectral scan of the particles
%
% Processing:
%   1. Dark-subtract the white reference, keep it per CCD row (corrects
%      the illumination fall-off along the slit) and bin it to the
%      sample's wavelength grid.
%   2. Remove cosmic-ray spikes from the sample.
%   3. Subtract the substrate background: median over strips for each CCD
%      row and wavelength, then median over rows for each strip, computed
%      with the particles masked out (two passes). This also removes the
%      sample's dark offset.
%   4. Normalize by the white reference, corrected for exposure time and
%      gain:  I = (S - S_bg)/(t_S*g_S) ./ ((W - D)/(t_W*g_W))
%   5. Build a particle image (mean intensity over imageRange) and choose
%      particles, either
%        'auto'   - found automatically, numbered brightest first, or
%        'manual' - clicked on the image with a crosshair, numbered in
%                   the order they are clicked.
%   6. Sum each particle's spectrum over a small box around it (its total
%      scattering) with a local background (surrounding ring) subtracted.
%   7. Classify each particle as a single particle or an aggregate, using
%      its spectrum (number of peaks, single-Lorentzian R^2, line width)
%      and its brightness compared with all the particles in the scan.
%      Particles with a neighbour nearby or an unusually large spot (likely
%      two particles too close to separate) are flagged as crowded.
%
% Outputs (next to the sample file):
%   <sample>_particles.xlsx  - sheet 1: column A wavelength, then one column
%                              per particle (input for fitSpectraLorentzian.m)
%                              'Particle Info': position, class, reasons
%                              'Aggregate Spectra': spectra left out of sheet 1
%   <sample>_processed.mat   - spectra, positions, classes, image, cube
%   <sample>_particles.png   - particle image with numbered particles
%
% Base MATLAB only; no toolboxes required.

clear; close all; clc;

%% ---------------- User settings ----------------
selectionMode   = 'manual';      % 'auto' = find particles automatically
                               % 'manual' = click particles on the image
snapToPeak      = true;        % manual: move each click to the nearest brightness peak

imageRange      = [650 950];   % nm, wavelengths averaged into the particle image
exportRange     = [500 1000];  % nm, wavelengths written to the Excel file
detectThreshold = 5;           % particles must exceed background by this many noise SDs
roiHalfSize     = [1 3];       % [rows strips] half-size of box around each particle
spectrumMode    = 'sum';       % 'sum' (total scattering) or 'mean' over the box
peakSeparation  = [3 5];       % [rows strips] minimum distance between particles
localBackground = true;        % subtract median of a ring around each particle
removeSpikes    = true;        % remove cosmic-ray spikes from the sample
spikeThreshold  = 8;           % spike = this many noise SDs above local median
gainCorrection  = true;        % divide by metadata 'gain' (assumes linear gain value)
stripAspect     = 'auto';      % strip step / CCD row size for display ('auto' or number)

% Single particle vs. aggregate. A particle is called an aggregate if ANY
% test fails. Brightness and width are compared with the median of all
% particles found automatically in the scan.
excludeAggregates = true;      % auto mode: leave aggregates out of the Excel spectra
fitWindow       = [650 950];   % nm, window for the quick Lorentzian fit
minSingleR2     = 0.9;         % single-Lorentzian fit must reach this R^2
peakProminence  = 0.1;         % a 2nd peak this prominent (fraction of max) = multiple peaks
edgeMargin      = 10;          % nm; peaks this close to a fitWindow edge skip the fit tests
maxWidthSD      = 4;           % line width more than this many robust SDs above median
maxBrightness   = 3;           % more than this many times the median brightness
maxSpotArea     = 2;           % crowded if half-max spot area > this many times the median

%% ---------------- Load files ----------------
white  = loadScan('Select the WHITE COUNT .mat file');
dark   = loadScan('Select the DARK COUNT .mat file', white.folder);
sample = loadScan('Select the SAMPLE (particle) .mat file', white.folder);

wl     = sample.wavelength;
nWl    = numel(wl);
nRows  = sample.scanSize(1);
nStrip = sample.scanSize(2);
fprintf('Sample: %d rows x %d strips, %d wavelengths\n', nRows, nStrip, nWl);

%% ---------------- 1. White reference (dark-subtracted, per CCD row) ----------------
darkSpec = mean(dark.intensities, 2);                        % detector dark spectrum
if numel(darkSpec) ~= numel(white.wavelength)
    error('White and dark files have different wavelength grids.');
end
wCube   = reshape(white.intensities, [], white.scanSize(1), white.scanSize(2));
whiteRow = mean(wCube, 3) - darkSpec;                        % nWl_W x nRows_W

% Match the sample's CCD rows
rowsS = sample.metadata.topPixel : sample.metadata.bottomPixel;
rowsW = white.metadata.topPixel  : white.metadata.bottomPixel;
[found, rowIdx] = ismember(rowsS, rowsW);
if all(found)
    whiteRow = whiteRow(:, rowIdx);
else
    warning('White reference does not cover all sample rows; using its mean spectrum.');
    whiteRow = repmat(mean(whiteRow, 2), 1, nRows);
end

% Match the sample's wavelength grid (sum binned detector pixels)
binFactor = numel(white.wavelength) / nWl;
if binFactor == round(binFactor) && ...
        max(abs(white.wavelength(binFactor:binFactor:end) - wl)) < 1e-3
    whiteRow = squeeze(sum(reshape(whiteRow, binFactor, nWl, nRows), 1));
    fprintf('Binned white reference by %d to match sample.\n', binFactor);
elseif binFactor ~= 1 || max(abs(white.wavelength - wl)) > 1e-3
    whiteRow = interp1(white.wavelength, whiteRow, wl, 'linear', NaN) * binFactor;
    fprintf('Interpolated white reference onto sample wavelengths.\n');
end

% Exposure / gain scaling
gW = 1; gS = 1;
if gainCorrection
    gW = white.metadata.gain;
    gS = sample.metadata.gain;
end
whiteRate = whiteRow / (white.metadata.exposure * gW);       % counts per s per gain
whiteRate(whiteRate < 0.01 * max(whiteRate(:))) = NaN;       % no lamp light -> undefined

%% ---------------- 2. Sample cube and spike removal ----------------
S = reshape(sample.intensities, nWl, nRows, nStrip);         % wavelength x row x strip

if removeSpikes
    med   = movmedian(S, 7, 1);
    resid = S - med;
    noise = 1.4826 * median(abs(resid), 1);                  % robust per-pixel noise
    spikes = resid > spikeThreshold * max(noise, 1);
    S(spikes) = med(spikes);
    fprintf('Removed %d cosmic-ray points.\n', nnz(spikes));
    clear med resid noise spikes
end

%% ---------------- 3. Substrate background subtraction ----------------
% Pass 1 uses all pixels; pass 2 repeats it with the particles masked out so
% crowded rows/strips are not over-subtracted.
inImage = wl >= imageRange(1) & wl <= imageRange(2);
mask = false(nRows, nStrip);
for pass = 1:2
    Sm = S;
    Sm(:, mask) = NaN;
    bgRow = median(Sm, 3, 'omitnan');          % per row & wavelength (slit illumination + dark)
    Sm = Sm - bgRow;
    bgStrip = median(Sm, 2, 'omitnan');        % per strip & wavelength (drift during scan)
    clear Sm
    if pass == 2                               % rows/strips that were fully masked
        bgRow   = fillmissing(bgRow, 'nearest', 2);
        bgStrip = fillmissing(bgStrip, 'nearest', 3);
    end
    if pass == 1
        % Mask everything noticeably brighter than the substrate
        im = squeeze(mean(S(inImage, :, :) - bgRow(inImage, :) - bgStrip(inImage, 1, :), 1));
        lvl = median(im(:));
        sd  = 1.4826 * median(abs(im(:) - lvl));
        mask = movmax(movmax(im > lvl + 3 * sd, 3, 1), 5, 2);   % grow mask a little
    end
end
S = S - bgRow - bgStrip;
fprintf('Background masked %.0f%% of pixels as particles.\n', 100 * mean(mask(:)));

%% ---------------- 4. White-light normalization ----------------
cube = S / (sample.metadata.exposure * gS) ./ whiteRate;    % normalized scattering
clear S

%% ---------------- 5. Particle image and automatic detection ----------------
img = squeeze(mean(cube(inImage, :, :), 1, 'omitnan'));      % rows x strips

% Light smoothing for detection only
[gy, gx] = ndgrid(-2:2, -3:3);
kernel = exp(-(gy.^2 / 2 + gx.^2 / 4.5));
imgSmooth = conv2(img, kernel / sum(kernel(:)), 'same');

bgPix   = imgSmooth(~mask);                                 % substrate only
bgLevel = median(bgPix);
bgNoise = 1.4826 * median(abs(bgPix - bgLevel));
isPeak  = imgSmooth == movmax(movmax(imgSmooth, 2*peakSeparation(1)+1, 1), ...
                              2*peakSeparation(2)+1, 2) ...
          & imgSmooth > bgLevel + detectThreshold * bgNoise;

% Drop particles whose box would run off the scan
edge = roiHalfSize;
isPeak([1:edge(1), end-edge(1)+1:end], :) = false;
isPeak(:, [1:edge(2), end-edge(2)+1:end]) = false;

[autoRow, autoStrip] = find(isPeak);
[~, order] = sort(imgSmooth(isPeak), 'descend');            % brightest first
autoRow = autoRow(order);  autoStrip = autoStrip(order);
fprintf('Automatic detection found %d particles.\n', numel(autoRow));

if ischar(stripAspect) && strcmp(stripAspect, 'auto')
    stripAspect = estimateAspect(img);
end
lo = prctile_(img(:), 1);
hi = prctile_(img(:), 99.7);

% Settings shared by the measurement functions
P = struct('roiHalfSize', roiHalfSize, 'ringPad', roiHalfSize + [2 3], ...
           'localBackground', localBackground, 'spectrumMode', spectrumMode, ...
           'inImage', inImage, 'imgSmooth', imgSmooth, 'bgLevel', bgLevel, ...
           'stripAspect', stripAspect, 'wl', wl, 'fitWindow', fitWindow, ...
           'peakProminence', peakProminence, 'edgeMargin', edgeMargin);

%% ---------------- 6. Choose particles ----------------
switch lower(selectionMode)
    case 'auto'
        pRow = autoRow;  pStrip = autoStrip;
    case 'manual'
        [pRow, pStrip] = selectParticles(img, imgSmooth, lo, hi, stripAspect, ...
                                         roiHalfSize, snapToPeak);
        if isempty(pRow)
            error('No particles were selected.');
        end
        fprintf('Selected %d particles manually.\n', numel(pRow));
    otherwise
        error('selectionMode must be ''auto'' or ''manual''.');
end
nP = numel(pRow);

%% ---------------- 7. Spectra and single/aggregate classification ----------------
fprintf('Measuring particles...\n');
feat = measureParticles(cube, pRow, pStrip, P);
if strcmpi(selectionMode, 'auto')
    ref = feat;                                              % same population
else
    ref = measureParticles(cube, autoRow, autoStrip, P);     % whole scan as reference
end
spectra = feat.spectra;

refBright = median(ref.brightness, 'omitnan');
refArea   = median(ref.area, 'omitnan');
okWidth   = ref.R2 >= minSingleR2 & ref.nPeaks == 1 & ~ref.atEdge;  % clean singles set the width scale
refWidth  = median(ref.fwhm(okWidth), 'omitnan');
sdWidth   = 1.4826 * median(abs(ref.fwhm(okWidth) - refWidth), 'omitnan');

testMulti  = feat.nPeaks >= 2;
% Fit-based tests are skipped when the peak is at the edge of the fit window
% (the width cannot be measured there).
testFit    = ~feat.atEdge & ~(feat.R2 >= minSingleR2);       % also catches failed fits
testBroad  = ~feat.atEdge & (feat.fwhm - refWidth) > maxWidthSD * sdWidth;
testBright = feat.brightness > maxBrightness * refBright;
testLarge  = feat.area > maxSpotArea * refArea;              % -> crowded, not aggregate
isAggregate = testMulti | testFit | testBroad | testBright;

fprintf(['Aggregate tests: >= 2 peaks, R^2 < %.2f, line width > %.0f meV, ' ...
         'brightness > %.1fx median\n'], minSingleR2, 1000 * (refWidth + maxWidthSD * sdWidth), ...
        maxBrightness);

reasons = strings(nP, 1);
for k = 1:nP
    r = strings(0);
    if testMulti(k),  r(end+1) = sprintf('%d peaks', feat.nPeaks(k)); end %#ok<SAGROW>
    if testFit(k),    r(end+1) = sprintf('R2 = %.2f', feat.R2(k)); end %#ok<SAGROW>
    if testBroad(k),  r(end+1) = sprintf('broad (%.0f meV)', 1000 * feat.fwhm(k)); end %#ok<SAGROW>
    if testBright(k), r(end+1) = sprintf('bright (%.1fx median)', feat.brightness(k) / refBright); end %#ok<SAGROW>
    if feat.atEdge(k) && ~isAggregate(k), r(end+1) = "peak at fit-window edge (width not tested)"; end %#ok<SAGROW>
    reasons(k) = strjoin(r, '; ');
end

% Crowded: another bright spot inside the background ring, or a spot so
% large it is probably two particles
anySpot = imgSmooth == movmax(movmax(imgSmooth, 3, 1), 5, 2) ...
          & imgSmooth > bgLevel + detectThreshold * bgNoise;
crowded = false(nP, 1);
for k = 1:nP
    rr = max(1, pRow(k) - P.ringPad(1))   : min(nRows,  pRow(k) + P.ringPad(1));
    cc = max(1, pStrip(k) - P.ringPad(2)) : min(nStrip, pStrip(k) + P.ringPad(2));
    crowded(k) = nnz(anySpot(rr, cc)) > 1 || testLarge(k);
end

% Numbering. Exported particles are numbered 1..N in the order they appear
% in the Excel file (= Spectrum # in the fitting script). In auto mode with
% excludeAggregates, aggregates are labelled A1, A2, ... instead.
exported = true(nP, 1);
if strcmpi(selectionMode, 'auto') && excludeAggregates
    exported = ~isAggregate;
end
labels = strings(nP, 1);
labels(exported)  = string(1:nnz(exported));
labels(~exported) = "A" + string(1:nnz(~exported));

fprintf('%d single particles, %d aggregates (%d crowded).\n', ...
        nnz(~isAggregate), nnz(isAggregate), nnz(crowded));
for k = find(isAggregate)'
    fprintf('  %-4s aggregate: %s\n', labels(k), reasons(k));
end

%% ---------------- 8. Save outputs ----------------
[~, baseName] = fileparts(sample.file);
outBase = fullfile(sample.folder, baseName);

% Excel file in the layout fitSpectraLorentzian.m expects
inExport    = wl >= exportRange(1) & wl <= exportRange(2);
wavelength  = wl(inExport);
intensities = spectra(inExport, exported);
xlsxFile = [outBase '_particles.xlsx'];
if exist(xlsxFile, 'file'), delete(xlsxFile); end            % avoid stale columns/sheets
writecell([{'Wavelength (nm)'}, cellstr("Particle " + labels(exported))'], xlsxFile, ...
          'Sheet', 1, 'Range', 'A1');
writematrix([wavelength, intensities], xlsxFile, 'Sheet', 1, 'Range', 'A2');

classNames = repmat("single", nP, 1);
classNames(isAggregate) = "aggregate";
info = table(labels, pRow, pStrip, classNames, crowded, feat.E0, 1000 * feat.fwhm, ...
             feat.R2, feat.nPeaks, feat.brightness / refBright, feat.area / refArea, reasons, ...
    'VariableNames', {'Particle', 'CCDRow', 'Strip', 'Class', 'Crowded', 'PeakEnergy_eV', ...
                      'FWHM_meV', 'R2', 'NumPeaks', 'BrightnessVsMedian', ...
                      'SpotAreaVsMedian', 'Notes'});
writetable(info, xlsxFile, 'Sheet', 'Particle Info');
if any(~exported)
    writecell([{'Wavelength (nm)'}, cellstr("Particle " + labels(~exported))'], xlsxFile, ...
              'Sheet', 'Aggregate Spectra', 'Range', 'A1');
    writematrix([wavelength, spectra(inExport, ~exported)], xlsxFile, ...
                'Sheet', 'Aggregate Spectra', 'Range', 'A2');
end

% Full results (all particles, including aggregates)
particleLabels    = labels;
particleRows      = pRow;
particleStrips    = pStrip;
particleCrowded   = crowded;
particleAggregate = isAggregate;
particleInfo      = info;
particleImage     = img;
correctedCube     = single(cube);                            % wavelength x row x strip
fullWavelength    = wl;
fullSpectra       = spectra;
save([outBase '_processed.mat'], 'wavelength', 'intensities', 'fullWavelength', ...
     'fullSpectra', 'particleLabels', 'particleRows', 'particleStrips', ...
     'particleCrowded', 'particleAggregate', 'particleInfo', 'particleImage', ...
     'correctedCube', 'imageRange', 'roiHalfSize', 'selectionMode');

%% ---------------- 9. Particle image ----------------
% Vertical orientation: CCD rows run vertically (row 1 at the top) and
% strips horizontally (strip 1 on the left).
fig = figure('Color', 'w', 'Position', [50 50 900 1300]);
if exist('theme', 'file'), theme(fig, 'light'); end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl);
showParticleImage(ax1, img, lo, hi, stripAspect);
title(ax1, {'Background-corrected image', sprintf('%d-%d nm', imageRange)});

ax2 = nexttile(tl);
showParticleImage(ax2, img, lo, hi, stripAspect);
hold(ax2, 'on');
colSingle    = [0 1 1];                                      % cyan
colCrowded   = [1 0.2 1];                                    % magenta
colAggregate = [0.4 1 0.2];                                  % lime green
for k = 1:nP
    c = colSingle;
    if crowded(k),     c = colCrowded;   end
    if isAggregate(k), c = colAggregate; end
    rectangle(ax2, 'Position', [pStrip(k)-roiHalfSize(2)-0.5, pRow(k)-roiHalfSize(1)-0.5, ...
              2*roiHalfSize(2)+1, 2*roiHalfSize(1)+1], 'EdgeColor', c, 'LineWidth', 0.8);
    if pRow(k) - roiHalfSize(1) > 5            % label above the box ...
        yLabel = pRow(k) - roiHalfSize(1) - 0.8;  vAlign = 'bottom';
    else                                       % ... or below it near the top edge
        yLabel = pRow(k) + roiHalfSize(1) + 0.8;  vAlign = 'top';
    end
    text(ax2, pStrip(k), yLabel, labels(k), 'Color', c, 'FontSize', 6, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', vAlign);
end
hold(ax2, 'off');
title(ax2, {sprintf('%d particles (%s selection)', nP, lower(selectionMode)), ...
            'cyan = single, magenta = crowded, green = aggregate'});

cb = colorbar(ax1, 'westoutside');
cb.Label.String = 'Normalized intensity (per pixel)';

exportgraphics(fig, [outBase '_particles.png'], 'Resolution', 200, ...
               'BackgroundColor', 'white');

fprintf(['Saved:\n  %s_particles.xlsx  (input for fitSpectraLorentzian.m)\n' ...
         '  %s_processed.mat\n  %s_particles.png\n'], outBase, outBase, outBase);

if isgraphics(fig), close(fig); end                          % image is saved as PNG

%% ---------------- Local functions ----------------
function [rows, strips] = selectParticles(img, imgSmooth, lo, hi, stripAspect, roiHalfSize, snap)
    % Click particles on the image with a crosshair.
    %   Left click          - add a particle
    %   Backspace / Delete
    %   or right click      - remove the last particle
    %   Enter (or Esc / q)  - finish
    [nRows, nStrip] = size(img);
    rows = zeros(0, 1);  strips = zeros(0, 1);  marks = {};

    fig = figure('Name', 'Select particles', 'NumberTitle', 'off', 'Color', 'w', ...
                 'Position', [100 40 650 1000]);
    if exist('theme', 'file'), theme(fig, 'light'); end
    ax = axes(fig);
    showParticleImage(ax, img, lo, hi, stripAspect);
    hold(ax, 'on');
    updateTitle();

    while true
        try
            [x, y, button] = ginput(1);                      % crosshair cursor
        catch
            break                                            % window was closed
        end
        if isempty(button) || any(button == [13 27 113])    % Enter / Esc / q
            break
        end
        if any(button == [3 8 127])                          % undo last
            if ~isempty(rows)
                delete(marks{end});  marks(end) = [];
                rows(end) = [];  strips(end) = [];
                updateTitle();
            end
            continue
        end
        if button ~= 1, continue; end

        r = round(y);  s = round(x);                         % x = strip, y = CCD row
        if r < 1 || r > nRows || s < 1 || s > nStrip, continue; end
        if snap                                              % nearest local maximum
            rr = max(1, r-2) : min(nRows, r+2);
            ss = max(1, s-3) : min(nStrip, s+3);
            [~, i] = max(reshape(imgSmooth(rr, ss), [], 1));
            [ir, is] = ind2sub([numel(rr) numel(ss)], i);
            r = rr(ir);  s = ss(is);
        end
        % Keep the whole box inside the scan
        r = min(max(r, roiHalfSize(1) + 1), nRows  - roiHalfSize(1));
        s = min(max(s, roiHalfSize(2) + 1), nStrip - roiHalfSize(2));
        if any(rows == r & strips == s)
            continue                                         % already selected
        end

        rows(end+1, 1) = r;  strips(end+1, 1) = s;          %#ok<AGROW>
        h1 = rectangle(ax, 'Position', [s-roiHalfSize(2)-0.5, r-roiHalfSize(1)-0.5, ...
                       2*roiHalfSize(2)+1, 2*roiHalfSize(1)+1], 'EdgeColor', 'c', 'LineWidth', 1);
        h2 = text(ax, s, r - roiHalfSize(1) - 0.8, num2str(numel(rows)), 'Color', 'c', ...
                  'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                  'VerticalAlignment', 'bottom');
        marks{end+1} = [h1 h2];                              %#ok<AGROW>
        updateTitle();
    end

    if isvalid(fig), close(fig); end

    function updateTitle()
        title(ax, {sprintf('Click particles in order (%d selected)', numel(rows)), ...
                   'Backspace / right-click = undo,   Enter = finish'}, 'Color', 'k');
    end
end

function showParticleImage(ax, img, lo, hi, stripAspect)
    % Particle image, vertical: x = strip (strip 1 on the left),
    % y = CCD row (row 1 at the top).
    imagesc(ax, img, [lo hi]);
    colormap(ax, blackOrangeWhite(256));
    set(ax, 'XDir', 'normal', 'YDir', 'reverse');
    daspect(ax, [1/stripAspect 1 1]);                        % physical aspect ratio
    xlabel(ax, 'Strip');  ylabel(ax, 'CCD row');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1);
end

function cmap = blackOrangeWhite(n)
    % Heat map: black -> orange -> white
    stops = [0 0 0; 1 0.5 0; 1 1 1];
    cmap  = interp1([0 0.5 1], stops, linspace(0, 1, n));
end

function F = measureParticles(cube, rows, strips, P)
    % Spectrum, brightness, spot size and quick Lorentzian fit per particle.
    nWl = size(cube, 1);  [~, nRows, nStrip] = size(cube);
    n = numel(rows);
    F.spectra    = zeros(nWl, n);
    F.brightness = nan(n, 1);
    F.area       = nan(n, 1);
    F.E0 = nan(n, 1);  F.fwhm = nan(n, 1);  F.R2 = nan(n, 1);  F.nPeaks = zeros(n, 1);
    F.atEdge = false(n, 1);

    hc = 1239.84193;
    inFit = P.wl >= P.fitWindow(1) & P.wl <= P.fitWindow(2);
    E = hc ./ P.wl(inFit);
    [E, eOrder] = sort(E);

    for k = 1:n
        % ---- Spectrum (box minus ring) ----
        r = rows(k) + (-P.roiHalfSize(1):P.roiHalfSize(1));
        c = strips(k) + (-P.roiHalfSize(2):P.roiHalfSize(2));
        box  = reshape(cube(:, r, c), nWl, []);
        spec = mean(box, 2, 'omitnan');
        if P.localBackground                     % ring, clipped to the scan
            rr = max(1, rows(k) - P.ringPad(1))   : min(nRows,  rows(k) + P.ringPad(1));
            cc = max(1, strips(k) - P.ringPad(2)) : min(nStrip, strips(k) + P.ringPad(2));
            inner = abs(rr' - rows(k)) <= P.roiHalfSize(1) + 1 & ...
                    abs(cc  - strips(k)) <= P.roiHalfSize(2) + 1;
            ring = reshape(cube(:, rr, cc), nWl, []);
            spec = spec - median(ring(:, ~inner(:)), 2, 'omitnan');
        end
        if strcmp(P.spectrumMode, 'sum')
            spec = spec * numel(r) * numel(c);
        end
        F.spectra(:, k) = spec;
        F.brightness(k) = mean(spec(P.inImage), 'omitnan');

        % ---- Spot area at half maximum (connected to the particle) ----
        F.area(k) = halfMaxArea(P.imgSmooth - P.bgLevel, rows(k), strips(k)) * P.stripAspect;

        % ---- Number of peaks and single-Lorentzian fit ----
        I = spec(inFit);  I = I(eOrder);
        good = ~isnan(I);
        if nnz(good) < 10, continue; end
        Is = movmean(I(good), 9);
        F.nPeaks(k) = nnz(islocalmax(Is, 'MinProminence', P.peakProminence * max(Is)));
        Eg = E(good);
        [~, iMax] = max(Is);
        lambdaMax = hc / Eg(iMax);
        F.atEdge(k) = lambdaMax < P.fitWindow(1) + P.edgeMargin || ...
                      lambdaMax > P.fitWindow(2) - P.edgeMargin;
        [F.E0(k), F.fwhm(k), F.R2(k)] = quickLorentzFit(E(good), I(good));
    end
end

function a = halfMaxArea(im, r0, s0)
    % Pixels above half the particle's peak that are connected to it.
    [nR, nS] = size(im);
    rr = max(1, r0-4) : min(nR, r0+4);
    ss = max(1, s0-8) : min(nS, s0+8);
    w  = im(rr, ss);
    above  = w >= im(r0, s0) / 2;
    region = false(size(w));
    region(r0 - rr(1) + 1, s0 - ss(1) + 1) = true;
    while true                                       % grow from the peak
        grown = conv2(double(region), ones(3), 'same') > 0 & above;
        if isequal(grown, region), break; end
        region = grown;
    end
    a = nnz(region);
end

function [E0, G, R2] = quickLorentzFit(E, I)
    % Same model as fitSpectraLorentzian.m: y0 + A*(G/2)^2/((E-E0)^2+(G/2)^2)
    E0 = NaN;  G = NaN;  R2 = NaN;
    scale = max(abs(I));
    if scale == 0, return; end
    In = I / scale;
    lor = @(p, x) p(1) + p(2) .* (abs(p(4))/2).^2 ./ ((x - p(3)).^2 + (abs(p(4))/2).^2);
    [pk, ip] = max(In);
    y0 = min(In);
    above = E(In >= y0 + (pk - y0)/2);
    g0 = max(above) - min(above);
    if g0 <= 0, g0 = 0.1; end
    p = fminsearch(@(p) sum((In - lor(p, E)).^2), [y0, pk - y0, E(ip), g0], ...
                   optimset('Display', 'off', 'MaxFunEvals', 5e3, 'MaxIter', 5e3));
    resid = In - lor(p, E);
    R2 = 1 - sum(resid.^2) / sum((In - mean(In)).^2);
    E0 = p(3);  G = abs(p(4));
    if p(2) <= 0 || E0 < E(1) || E0 > E(end)
        R2 = min(R2, 0);                             % no peak in the window = failed fit
    end
end

function s = loadScan(prompt, startFolder)
    if nargin < 2, startFolder = pwd; end
    [f, p] = uigetfile({'*.mat', 'MAT files (*.mat)'}, prompt, startFolder);
    if isequal(f, 0)
        error('No file selected.');
    end
    s = load(fullfile(p, f), 'wavelength', 'intensities', 'scanSize', 'metadata');
    if ~all(isfield(s, {'wavelength', 'intensities', 'scanSize', 'metadata'})) || isempty(s.scanSize)
        error('%s is not a converted scan file (run convertTDMStoMAT.m first).', f);
    end
    s.wavelength  = double(s.wavelength(:));
    s.intensities = double(s.intensities);
    s.file   = f;
    s.folder = p;
    fprintf('%-30s %s\n', prompt(12:end), f);
end

function a = estimateAspect(img)
    % Strip-step / row-size ratio from the half-width of the image
    % autocorrelation along each axis (particles are round in reality).
    x  = img - mean(img(:));
    ac = real(ifft2(abs(fft2(x)).^2));
    ac = ac / ac(1, 1);
    hwRow   = halfWidth(ac(:, 1));
    hwStrip = halfWidth(ac(1, :)');
    a = min(max(hwRow / hwStrip, 0.2), 5);
end

function hw = halfWidth(v)
    i = find(v < 0.5, 1);
    if isempty(i) || i < 2
        hw = 1;
    else
        hw = (i - 2) + (v(i-1) - 0.5) / (v(i-1) - v(i));
    end
end

function q = prctile_(x, p)
    % Percentile without the Statistics Toolbox
    x = sort(x(~isnan(x)));
    q = x(max(1, min(numel(x), round(p / 100 * numel(x)))));
end
