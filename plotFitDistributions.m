%% plotFitDistributions.m
% Plots the distributions of resonance energy and line width from the
% fit_results.csv file written by fitSpectraLorentzian.m. Only the CSV is
% needed. (fit_results.xlsx and the 'Fit Results' sheet of
% <sample>_particles.xlsx can also be selected.)
%
% Optional: with excludeAggregates = true, particles classified as
% aggregates by processSpectra.m are left out. The classification is read
% from the 'Particle Info' sheet of <sample>_particles.xlsx, looked for in
% the same folder as the CSV or the folder above.
%
% Figures:
%   1. Histogram of resonance energy (eV)
%   2. Cumulative distribution (CDF) of resonance energy
%   3. Histogram of line width (FWHM, meV)
%   4. Cumulative distribution (CDF) of line width
%
% Each figure is saved as a PNG next to the input file.
% Base MATLAB only; no toolboxes required.

clear; close all; clc;

%% ---------------- User settings ----------------
minR2       = 0.9;        % exclude fits with R^2 below this (set to -Inf to keep all)
fitWindow   = [500 950];  % nm, must match lambdaLimits in fitSpectraLorentzian.m
edgeTol     = 2;          % nm, exclude peaks this close to the window edges
excludeAggregates = false;% true = leave out aggregates (needs <sample>_particles.xlsx)
excludeCrowded    = false;% true = leave out crowded particles (same file)
widthUnits  = 'meV';      % line width units: 'meV' or 'eV'
energyBin   = 0.02;       % histogram bin width for resonance energy (eV); [] = automatic
widthBin    = 5;          % histogram bin width for line width (in widthUnits); [] = automatic
saveFigures = true;       % save each figure as a PNG

%% ---------------- Load fit results ----------------
[fileName, filePath] = uigetfile({'*.csv', 'Fit results CSV (*.csv)'; ...
                                  '*.xlsx;*.xls', 'Excel files (*.xlsx, *.xls)'}, ...
                                 'Select the fit_results.csv file');
if isequal(fileName, 0)
    error('No file selected.');
end
dataFile = fullfile(filePath, fileName);
[~, baseName, ext] = fileparts(fileName);

readOpts = {'VariableNamingRule', 'modify'};
if strcmpi(ext, '.csv')
    T = readtable(dataFile, readOpts{:});
else
    sheets = sheetnames(dataFile);
    if any(sheets == "Fit Results")
        sheet = "Fit Results";                       % <sample>_particles.xlsx
    else
        sheet = sheets(1);                           % fit_results.xlsx
    end
    T = readtable(dataFile, 'Sheet', sheet, readOpts{:});
    fprintf('Reading sheet "%s" of %s\n', sheet, fileName);
end
T = standardizeColumns(T, fileName);

%% ---------------- Particle classification (single / aggregate) ----------------
% Only used when excludeAggregates or excludeCrowded is on.
info = [];
if excludeAggregates || excludeCrowded
    info = findParticleInfo(dataFile, ext);
end
isAgg   = false(height(T), 1);
isCrowd = false(height(T), 1);
if ~isempty(info)
    % Spectrum k in the fit results = particle labelled "k". Check that both
    % files come from the same run (same particles, matching peak energies).
    [found, idx] = ismember(string(T.Spectrum), info.Particle);
    nFitted = nnz(~startsWith(info.Particle, "A"));        % A# = aggregates not fitted
    dE = abs(T.PeakEnergy_eV(found) - info.PeakEnergy_eV(idx(found)));
    if ~all(found) || nFitted ~= height(T) || median(dE, 'omitnan') > 0.01
        warning(['The Particle Info sheet does not match these fit results (re-run ' ...
                 'fitSpectraLorentzian.m on the current particles file?). ' ...
                 'Aggregates are not excluded.']);
    else
        isAgg   = info.Class(idx) == "aggregate";
        isCrowd = logical(info.Crowded(idx));
    end
end

% Keep successful, good-quality fits whose peak lies inside the window
hc = 1239.84193;                                   % eV*nm
peakNm  = hc ./ T.PeakEnergy_eV;
failed  = isnan(T.PeakEnergy_eV) | isnan(T.FWHM_eV);
noPeak  = ~failed & T.Amplitude <= 0;              % flat/degenerate fit
atEdge  = ~failed & ~noPeak & (peakNm < fitWindow(1) + edgeTol | peakNm > fitWindow(2) - edgeTol);
poorFit = ~failed & ~noPeak & ~atEdge & T.R2 < minR2;
aggr    = excludeAggregates & isAgg & ~(failed | noPeak | atEdge | poorFit);
crowd   = excludeCrowded & isCrowd & ~(failed | noPeak | atEdge | poorFit | aggr);
keep    = ~(failed | noPeak | atEdge | poorFit | aggr | crowd);

fprintf('Read %d spectra from %s (fitting window %g-%g nm).\n', ...
        height(T), fileName, fitWindow);
fprintf('Using %d of %d spectra.\n', nnz(keep), height(T));
reportExcluded(T.Spectrum, failed,  'failed fit');
reportExcluded(T.Spectrum, noPeak,  'no peak (amplitude <= 0)');
reportExcluded(T.Spectrum, atEdge,  sprintf('peak within %g nm of the fitting window edge', edgeTol));
reportExcluded(T.Spectrum, poorFit, sprintf('R^2 < %g', minR2));
reportExcluded(T.Spectrum, aggr,    'aggregate');
reportExcluded(T.Spectrum, crowd,   'crowded');
if nnz(keep) < 2
    error('Not enough spectra left to plot.');
end

energy = T.PeakEnergy_eV(keep);
width  = T.FWHM_eV(keep);
if strcmp(widthUnits, 'meV')
    width = 1000 * width;
end

%% ---------------- Figures ----------------
outBase = fullfile(filePath, baseName);

figs = gobjects(1, 4);
figs(1) = plotHistogram(energy, 'Energy of Resonance (eV)', 'Energy of Resonance', '%.3f eV', ...
                        energyBin, saveFigures, [outBase '_resonance_histogram.png']);
figs(2) = plotCDF(energy, 'Energy of Resonance (eV)', 'Energy of Resonance CDF', ...
                  saveFigures, [outBase '_resonance_cdf.png']);
figs(3) = plotHistogram(width, sprintf('Line Width (%s)', widthUnits), 'Line Width', ...
                        ['%.' num2str(3 * strcmp(widthUnits, 'eV')) 'f ' widthUnits], ...
                        widthBin, saveFigures, [outBase '_linewidth_histogram.png']);
figs(4) = plotCDF(width, sprintf('Line Width (%s)', widthUnits), 'Line Width CDF', ...
                  saveFigures, [outBase '_linewidth_cdf.png']);

if saveFigures
    fprintf('Figures saved to:\n  %s\n', filePath);
    close(figs(isgraphics(figs)));   % windows closed once the PNGs are written
end

%% ---------------- Local functions ----------------
function T = standardizeColumns(T, fileName)
    % Map the column names used by the different output files onto one set:
    % Spectrum, PeakEnergy_eV, FWHM_eV, R2, Amplitude.
    v = T.Properties.VariableNames;
    if ~ismember('Spectrum', v)
        T.Spectrum = (1:height(T))';
    end
    if ~ismember('R2', v) && ismember('R_squared', v)
        T.R2 = T.R_squared;
    end
    if ~ismember('FWHM_eV', v) && ismember('FWHM_meV', v)
        T.FWHM_eV = T.FWHM_meV / 1000;
    end
    if ~ismember('Amplitude', v)
        if ismember('PeakIntensity', v)
            T.Amplitude = T.PeakIntensity;           % 'Fit Results' sheet: offset + amplitude
        else
            T.Amplitude = ones(height(T), 1);        % amplitude check not possible
        end
    end
    required = {'PeakEnergy_eV', 'FWHM_eV', 'R2'};
    missing = setdiff(required, T.Properties.VariableNames);
    if ~isempty(missing)
        error('%s is missing columns: %s', fileName, strjoin(missing, ', '));
    end
end

function info = findParticleInfo(dataFile, ext)
    % 'Particle Info' sheet from processSpectra.m: in the selected Excel file,
    % or in a single *_particles.xlsx next to it or one folder up (fit
    % results are saved in a spectra_fits subfolder).
    info = [];
    candidates = string.empty;
    if ~strcmpi(ext, '.csv') && any(sheetnames(dataFile) == "Particle Info")
        candidates = string(dataFile);
    else
        here = fileparts(dataFile);
        for folder = [string(here), string(fileparts(here))]
            f = dir(fullfile(folder, '*_particles.xlsx'));
            f = f(~startsWith({f.name}, '~$'));              % skip Excel lock files
            if numel(f) == 1
                candidates = string(fullfile(f.folder, f.name));
                break
            elseif numel(f) > 1
                warning(['Several *_particles.xlsx files in %s; not using particle ' ...
                         'classification.'], folder);
                return
            end
        end
    end
    if isempty(candidates) || ~any(sheetnames(candidates) == "Particle Info")
        fprintf('No particle classification found; aggregates cannot be excluded.\n');
        return
    end
    info = readtable(candidates, 'Sheet', 'Particle Info', 'TextType', 'string');
    info.Particle = string(info.Particle);
    info.Class    = string(info.Class);
    fprintf('Particle classification from %s\n', candidates);
end

function reportExcluded(ids, mask, reason)
    if any(mask)
        fprintf('  Excluded %d (%s): spectra %s\n', nnz(mask), reason, ...
                strjoin(string(ids(mask)'), ', '));
    end
end

function ax = newSquareAxes()
    fig = figure('Color', 'w', 'Position', [100 100 600 600]);
    if exist('theme', 'file'), theme(fig, 'light'); end   % override dark theme
    ax = axes(fig, 'Color', 'w');
    hold(ax, 'on');
end

function styleAxes(ax, xLabel, yLabel, titleText)
    pbaspect(ax, [1 1 1]);                          % square plot area
    box(ax, 'on');                                  % full black border
    set(ax, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5, ...
            'FontSize', 12, 'TickDir', 'in');
    xlabel(ax, xLabel, 'Color', 'k');
    ylabel(ax, yLabel, 'Color', 'k');
    title(ax, titleText, 'Color', 'k');
    hold(ax, 'off');
end

function fig = plotHistogram(x, xLabel, titleText, fmt, binWidth, saveFig, outFile)
    ax = newSquareAxes();
    if isempty(binWidth)
        histogram(ax, x, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.8);
    else
        edges = (floor(min(x) / binWidth) : ceil(max(x) / binWidth)) * binWidth;
        if numel(edges) < 2, edges = [edges, edges(end) + binWidth]; end
        histogram(ax, x, edges, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.8);
    end
    ylim(ax, [0, max(ylim(ax)) * 1.2]);             % room for the statistics box
    stats = sprintf(['N = %d\nMean = ' fmt '\nMedian = ' fmt '\nSD = ' fmt], ...
                    numel(x), mean(x), median(x), std(x));
    text(ax, 0.97, 0.97, stats, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top', 'FontSize', 11, 'Color', 'k', ...
         'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 4);
    styleAxes(ax, xLabel, 'Count', titleText);
    fig = ax.Parent;
    saveFigure(fig, saveFig, outFile);
end

function fig = plotCDF(x, xLabel, titleText, saveFig, outFile)
    ax = newSquareAxes();
    xs = sort(x(:));
    n  = numel(xs);
    % Empirical CDF as a step function, starting from 0
    stairs(ax, [xs(1); xs], [0; (1:n)' / n], 'b-', 'LineWidth', 2);
    ylim(ax, [0 1.05]);
    pad = 0.05 * (xs(end) - xs(1));
    xlim(ax, [xs(1) - pad, xs(end) + pad]);
    grid(ax, 'on');
    styleAxes(ax, xLabel, 'Cumulative Probability', titleText);
    fig = ax.Parent;
    saveFigure(fig, saveFig, outFile);
end

function saveFigure(fig, saveFig, outFile)
    if saveFig
        exportgraphics(fig, outFile, 'Resolution', 300, 'BackgroundColor', 'white');
    end
end
