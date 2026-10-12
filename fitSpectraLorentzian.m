%% fitSpectraLorentzian.m
% Plots single-particle spectra from an Excel file and fits each spectrum
% with a Lorentzian line shape.
%
% Excel layout:
%   Column A        -> wavelength (nm)
%   Columns B, C... -> intensity for each particle
%
% For each particle this script:
%   - converts wavelength to photon energy (eV)
%   - plots intensity as blue points vs. energy
%   - fits a Lorentzian (red) over the 650-950 nm window
%   - sets each plot's axes automatically around the fitted peak
%   - saves the figure as a PNG
% A summary table of fit parameters is saved to fit_results.csv, and the
% peak intensity, FWHM, peak energy and R^2 of each spectrum are written to
% a "Fit Results" sheet in the original Excel file.
%
% Uses only base MATLAB (fminsearch); no toolboxes required.

clear; close all; clc;

%% ---------------- User settings ----------------
lambdaLimits = [500 950];     % wavelength window used for fitting (nm)
xRangeFWHM   = 3;             % x-axis spans the peak +/- this many FWHMs
yHeadroom    = 0.15;          % y-axis extends this fraction above the peak
saveFigures  = true;          % save each plot as a PNG
resultsSheet = 'Fit Results'; % name of the sheet added to the Excel file
hc           = 1239.84193;    % h*c in eV*nm  (E[eV] = hc / lambda[nm])

%% ---------------- Load data ----------------
[fileName, filePath] = uigetfile({'*.xlsx;*.xls', 'Excel files'}, ...
                                 'Select the spectral data file');
if isequal(fileName, 0)
    error('No file selected.');
end
dataFile = fullfile(filePath, fileName);

raw = readmatrix(dataFile, 'Sheet', 1); % spectra are on the first sheet; header rows skipped
raw = raw(~isnan(raw(:,1)), :);        % drop rows with no wavelength

wavelength  = raw(:,1);                % nm
intensities = raw(:,2:end);            % one column per particle
nSpectra    = size(intensities, 2);

%% ---------------- Convert to eV ----------------
energy  = hc ./ wavelength;                    % eV
fitLimits = sort(hc ./ lambdaLimits);          % [950 nm, 650 nm] -> [~1.305, ~1.907] eV

% Points inside the fitting window
inWindow = energy >= fitLimits(1) & energy <= fitLimits(2);

%% ---------------- Output folder ----------------
outDir = fullfile(filePath, 'spectra_fits');
if saveFigures && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ---------------- Lorentzian model ----------------
% p = [y0, A, E0, Gamma]  (Gamma = FWHM in eV)
lorentz = @(p, E) p(1) + p(2) .* (abs(p(4))/2).^2 ./ ((E - p(3)).^2 + (abs(p(4))/2).^2);

fitOpts = optimset('Display', 'off', 'MaxFunEvals', 2e4, 'MaxIter', 2e4, ...
                   'TolX', 1e-10, 'TolFun', 1e-12);

results = nan(nSpectra, 6);   % [y0, A, E0, Gamma, lambda0, R^2]
prevFig = [];                 % handle of the previously drawn figure

%% ---------------- Loop over particles ----------------
for k = 1:nSpectra
    y = intensities(:, k);
    valid = inWindow & ~isnan(y);
    E = energy(valid);
    I = y(valid);

    % Sort by energy so lines draw cleanly
    [E, order] = sort(E);
    I = I(order);

    % ---- Fit (on normalized data for better conditioning) ----
    scale = max(abs(I));
    if scale == 0 || numel(I) < 4
        warning('Spectrum %d: not enough usable data, skipping fit.', k);
        pFit = nan(1,4);
    else
        In = I / scale;

        % Initial guesses
        [peakVal, iPeak] = max(In);
        y0Guess = min(In);
        AGuess  = peakVal - y0Guess;
        E0Guess = E(iPeak);
        halfMax = y0Guess + AGuess/2;
        above   = E(In >= halfMax);
        GGuess  = max(above) - min(above);
        if GGuess <= 0, GGuess = 0.1; end

        p0   = [y0Guess, AGuess, E0Guess, GGuess];
        sse  = @(p) sum((In - lorentz(p, E)).^2);
        pFit = fminsearch(sse, p0, fitOpts);

        % Undo normalization
        pFit(1:2) = pFit(1:2) * scale;
        pFit(4)   = abs(pFit(4));
    end

    % Goodness of fit
    if all(~isnan(pFit))
        resid = I - lorentz(pFit, E);
        R2 = 1 - sum(resid.^2) / sum((I - mean(I)).^2);
    else
        R2 = NaN;
    end
    results(k,:) = [pFit, hc / pFit(3), R2];

    % ---- Plot ----
    fig = figure('Color', 'w', 'Position', [100 100 600 600]);
    if exist('theme', 'file')              % R2025a+: override dark theme
        theme(fig, 'light');
    end
    ax  = axes(fig, 'Color', 'w');
    hold(ax, 'on');

    % ---- Axis limits from the fitted peak ----
    % x: peak energy +/- xRangeFWHM * FWHM, kept inside the fitted data
    % y: from 0 (or the lowest visible point if negative) to above the peak
    if all(~isnan(pFit)) && pFit(3) > E(1) && pFit(3) < E(end)
        halfSpan = xRangeFWHM * pFit(4);
        xLimits  = [max(E(1), pFit(3) - halfSpan), min(E(end), pFit(3) + halfSpan)];
    else
        xLimits  = [E(1), E(end)];         % fit failed or peak outside window
    end
    shown = E >= xLimits(1) & E <= xLimits(2);
    EFine = linspace(xLimits(1), xLimits(2), 1000);
    yFit  = lorentz(pFit, EFine);
    yTop  = max([I(shown); yFit(:)]);
    yBot  = min([0; I(shown)]);
    yLimits = [yBot - yHeadroom * (yTop - yBot) * (yBot < 0), ...
               yTop + yHeadroom * (yTop - yBot)];

    hData = plot(ax, E(shown), I(shown), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    hFit  = plot(ax, EFine, yFit, 'r-', 'LineWidth', 2);

    xlim(ax, xLimits);
    ylim(ax, yLimits);
    pbaspect(ax, [1 1 1]);                 % square plot area
    box(ax, 'on');                         % full black border
    set(ax, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5, ...
            'FontSize', 12, 'TickDir', 'in');

    xlabel(ax, 'Energy (eV)', 'Color', 'k');
    ylabel(ax, 'Intensity', 'Color', 'k');
    title(ax, sprintf('Spectrum %d', k), 'Color', 'k');
    legend(ax, [hData, hFit], {'Data', 'Lorentzian Fit'}, ...
           'Location', 'northeast', 'Box', 'on', ...
           'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    hold(ax, 'off');

    if saveFigures
        exportgraphics(fig, fullfile(outDir, sprintf('Spectrum_%d.png', k)), ...
                       'Resolution', 300, 'BackgroundColor', 'white');
    end

    % Close the previous spectrum's window now that this one is drawn
    drawnow;
    if ~isempty(prevFig) && isvalid(prevFig)
        close(prevFig);
    end
    prevFig = fig;
end

%% ---------------- Summary table ----------------
fitTable = array2table(results, 'VariableNames', ...
    {'Offset', 'Amplitude', 'PeakEnergy_eV', 'FWHM_eV', 'PeakWavelength_nm', 'R2'});
fitTable.FWHM_meV = 1000 * fitTable.FWHM_eV;
fitTable = addvars(fitTable, (1:nSpectra)', 'Before', 1, 'NewVariableNames', 'Spectrum');
disp(fitTable);

if saveFigures
    writetable(fitTable, fullfile(outDir, 'fit_results.csv'));
    fprintf('Figures and fit_results.csv saved to:\n  %s\n', outDir);
end

%% ---------------- Write results to a new sheet in the Excel file ----------------
% Peak intensity = fitted Lorentzian height at E0 (offset + amplitude)
summaryTable = table((1:nSpectra)', ...
                     results(:,1) + results(:,2), ...
                     results(:,4), ...
                     results(:,3), ...
                     results(:,6), ...
    'VariableNames', {'Spectrum', 'PeakIntensity', 'FWHM_eV', 'PeakEnergy_eV', 'R_squared'});

try
    writetable(summaryTable, dataFile, 'Sheet', resultsSheet, ...
               'WriteMode', 'overwritesheet');
    fprintf('Fit summary written to sheet "%s" in:\n  %s\n', resultsSheet, dataFile);
catch err
    warning(['Could not write to "%s" (is it open in Excel?).\n%s\n' ...
             'Fit results are still available in fitTable and fit_results.csv.'], ...
            dataFile, err.message);
end

%% ---------------- Close the last figure window ----------------
% Kept open if figures are not being saved, so the plot is not lost.
if saveFigures && ~isempty(prevFig) && isvalid(prevFig)
    close(prevFig);
end
