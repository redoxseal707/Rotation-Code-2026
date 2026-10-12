%% convertTDMStoMAT.m
% Converts spectral data saved by LabVIEW as .tdms files into .mat files
% for use with the spectral fitting script (fitSpectraLorentzian.m).
%
% Expected TDMS layout:
%   Group "Spectra" -> one channel per spectrum (raw CCD counts)
%   Group "Info"    -> channel "wvlths" with the wavelength axis (nm)
%   File properties -> acquisition settings (exposure, gain, strips, ...)
%
% Each selected file.tdms is saved as file.mat containing:
%   wavelength   - [nWavelengths x 1] wavelength axis (nm)
%   intensities  - [nWavelengths x nSpectra] counts, one column per spectrum
%   scanSize     - [nRows nStrips] scan grid, so that
%                  cube = reshape(intensities, [], nRows, nStrips)
%                  gives intensities at (wavelength, CCD row, strip)
%   metadata     - file-level acquisition properties
%   info         - all other channels in the file (e.g. Info group)
%   sourceFile   - full path of the original .tdms file
%
% Requires readTDMS.m in the same folder (no toolboxes needed).

clear; clc;

%% ---------------- User settings ----------------
spectraGroup      = 'Spectra';   % group holding one channel per spectrum
wavelengthChannel = 'wvlths';    % channel holding the wavelength axis
outputFolder      = '';          % '' = save next to each .tdms file

%% ---------------- Select files ----------------
[fileNames, filePath] = uigetfile({'*.tdms', 'TDMS files (*.tdms)'}, ...
                                  'Select TDMS file(s) to convert', ...
                                  'MultiSelect', 'on');
if isequal(fileNames, 0)
    error('No file selected.');
end
fileNames = cellstr(fileNames);

%% ---------------- Convert each file ----------------
for f = 1:numel(fileNames)
    sourceFile = fullfile(filePath, fileNames{f});
    fprintf('Reading %s ... ', fileNames{f});
    tdms = readTDMS(sourceFile);

    % ---- Spectra ----
    groupNames = {tdms.groups.name};
    gSpec = find(strcmp(groupNames, spectraGroup), 1);
    if isempty(gSpec)
        % Fall back to the group with the most channels
        [~, gSpec] = max(arrayfun(@(g) numel(g.channels), tdms.groups));
        warning('No "%s" group in %s; using group "%s".', ...
                spectraGroup, fileNames{f}, groupNames{gSpec});
    end
    chans = tdms.groups(gSpec).channels;

    % Order channels by their LabVIEW array column, if recorded
    cols = arrayfun(@(c) getfield_or(c.properties, 'NI_ArrayColumn', NaN), chans);
    if all(~isnan(cols))
        [~, order] = sort(cols);
        chans = chans(order);
    end

    nPoints = arrayfun(@(c) numel(c.data), chans);
    if any(nPoints ~= nPoints(1))
        error('Spectra in %s have different lengths.', fileNames{f});
    end
    intensities = double([chans.data]);          % nWavelengths x nSpectra

    % ---- Wavelength axis ----
    wavelength = [];
    info = struct();
    for g = 1:numel(tdms.groups)
        if g == gSpec, continue; end
        gField = matlab.lang.makeValidName(tdms.groups(g).name);
        chNames = matlab.lang.makeUniqueStrings( ...
                    matlab.lang.makeValidName({tdms.groups(g).channels.name}));
        for c = 1:numel(tdms.groups(g).channels)
            ch = tdms.groups(g).channels(c);
            info.(gField).(chNames{c}) = ch.data;
            if isempty(wavelength) && strcmp(ch.name, wavelengthChannel)
                wavelength = double(ch.data(:));
            end
        end
    end
    if isempty(wavelength)
        error('No "%s" wavelength channel found in %s.', wavelengthChannel, fileNames{f});
    end
    if numel(wavelength) ~= size(intensities, 1)
        error('Wavelength axis (%d points) does not match spectra (%d points) in %s.', ...
              numel(wavelength), size(intensities, 1), fileNames{f});
    end

    % ---- Scan geometry (spectra are stored strip by strip) ----
    metadata = tdms.properties;
    scanSize = [];
    if all(isfield(metadata, {'strips', 'topPixel', 'bottomPixel'}))
        nRows   = metadata.bottomPixel - metadata.topPixel + 1;
        nStrips = metadata.strips;
        if nRows * nStrips == size(intensities, 2)
            scanSize = [nRows nStrips];
        end
    end

    % ---- Save ----
    outDir = outputFolder;
    if isempty(outDir), outDir = filePath; end
    [~, baseName] = fileparts(fileNames{f});
    outFile = fullfile(outDir, [baseName '.mat']);

    vars = {'wavelength', 'intensities', 'scanSize', 'metadata', 'info', 'sourceFile'};
    if numel(intensities) * 8 > 2^31
        save(outFile, vars{:}, '-v7.3');         % > 2 GB needs v7.3
    else
        save(outFile, vars{:});
    end

    fprintf('%d spectra x %d wavelengths (%.1f-%.1f nm) -> %s\n', ...
            size(intensities, 2), size(intensities, 1), ...
            min(wavelength), max(wavelength), outFile);
end

fprintf('Done. Converted %d file(s).\n', numel(fileNames));

%% ---------------- Local functions ----------------
function v = getfield_or(s, name, default)
    if isfield(s, name)
        v = double(s.(name));
    else
        v = default;
    end
end
