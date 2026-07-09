function [distProfileGFP, distProfileRFP, pairs] = mitoIntensityFS(pathGFP, pathRFP, outputDir, stageName, opts)
% mitoIntensityFS
%
% Single-experiment processing pipeline for GFP and/or mCherry spatial
% fluorescence intensity profiles from Fiji-generated single-cell CSV files.
%
% -------------------------------------------------------------------------
% OVERVIEW
% -------------------------------------------------------------------------
% This function is the single-experiment profile-processing layer of the
% pooled mitochondrial spatial-analysis pipeline.
%
% For one selected meiotic stage, the function:
%
%       - discovers GFP and/or mCherry profile CSV files
%       - parses treatment and timepoint information from filenames
%       - filters files by the requested experimental metadata
%       - pairs GFP and mCherry profiles when both channels are analysed
%       - normalises position to cell length
%       - optionally normalises fluorescence intensity per cell
%       - divides each cell into equally spaced position bins
%       - calculates the mean fluorescence intensity within each bin
%
% The resulting matrices contain one spatial intensity profile per cell and
% can be analysed directly or concatenated across experiments by
% mitoIntensityFS_pooled.
%
% -------------------------------------------------------------------------
% PIPELINE POSITION
% -------------------------------------------------------------------------
%
% runPooledPol
%       |
%       +-- mitoIntensityFS_pooled
%               |
%               +-- mitoIntensityFS
%
% This function performs the core single-cell profile processing.
%
% It does NOT pool experiments, calculate front/back polarisation, or
% calculate the Regional Segregation Score.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% pathGFP : character vector or string-compatible path
%     Directory containing GFP single-cell profile CSV files.
%
%     The path is used when opts.channelMode is:
%
%         "both"
%         "gfp"
%
% pathRFP : character vector or string-compatible path
%     Directory containing mCherry/RFP single-cell profile CSV files.
%
%     The path is used when opts.channelMode is:
%
%         "both"
%         "rfp"
%
% outputDir : character vector or string-compatible path
%     Directory used for pairing reports and optional profile figures.
%
% stageName : character vector or string-compatible stage name
%     Descriptive stage identifier used in console reporting, pairing
%     reports, and figure filenames.
%
% opts : structure containing file-selection, profile-processing, and
%     plotting options.
%
% -------------------------------------------------------------------------
% EXPECTED CSV DATA
% -------------------------------------------------------------------------
% Each profile CSV must contain at least two numeric data columns.
%
% The first two numeric columns are interpreted as:
%
%       Column 1 = spatial position
%       Column 2 = fluorescence intensity
%
% Additional columns are ignored by the profile-binning calculation.
%
% Files are loaded using importdata.
%
% -------------------------------------------------------------------------
% REQUIRED ON MATLAB PATH
% -------------------------------------------------------------------------
%
%       getAllFiles.m
%
% OPTIONAL:
%
%       shadedErrorBar.m
%
% If shadedErrorBar is available, single-experiment mean profile plots
% include SEM shading. Otherwise, the mean profile is plotted without the
% shaded SEM region.
%
% -------------------------------------------------------------------------
% CORE FEATURES
% -------------------------------------------------------------------------
%
% 1) CHANNEL-AWARE PROFILE PROCESSING
%    - opts.channelMode controls which channel directories are analysed:
%
%          "both"
%              Paired GFP and mCherry analysis.
%
%          "gfp"
%              GFP-only analysis.
%
%          "rfp"
%              mCherry/RFP-only analysis.
%
%    - The unused output profile matrix remains empty in single-channel mode.
%
%    - Channel mode is validated before file processing.
%
% 2) CSV FILE DISCOVERY
%    - CSV files are collected from the selected channel path using:
%
%          getAllFiles(<channelPath>, '*.csv')
%
%    - File discovery is performed independently for GFP and mCherry/RFP.
%
%    - No files are collected for channels disabled by opts.channelMode.
%
% 3) FILENAME METADATA PARSING
%    - Treatment, timepoint, channel, series number, ROI number, and a
%      normalised pairing key are parsed from each CSV filename.
%
%    - Treatment parsing:
%
%          filenames containing "H2O2"
%              -> treatment = "H2O2"
%
%          tokens beginning with "Veh" or "Vehicle"
%              -> treatment = "Veh"
%
%          otherwise
%              -> treatment = "UNKNOWN"
%
%    - Timepoint parsing currently recognises:
%           Nmin/NminT and exact TN tokens are recognised:
%
%          "30min"
%              -> timepoint = "30min"
%
%          "60min"
%              -> timepoint = "60min"
%
%          "5min", "0min", T0-style tokens, or T5-style tokens
%              -> timepoint = "5min"
%
%    - IMPORTANT:
%      Any filename without one of the explicitly recognised timepoint
%      tokens currently defaults to:
%
%          timepoint = "UNKNOWN"
%
%    - Series numbers are parsed from tokens of the form:
%
%          seriesN
%
%    - ROI numbers are parsed from tokens of the form:
%
%          roiN
%
%    - A leading token of the form ChannelN is used to record the channel
%      number when present.
%
% 4) CONDITION AND TIMEPOINT FILTERING
%    - opts.condition can restrict files by parsed treatment.
%
%    - opts.timepoint can restrict files by parsed timepoint.
%
%    - Common values:
%
%          opts.condition = "Veh"
%          opts.condition = "H2O2"
%
%          opts.timepoint = "5min"
%          opts.timepoint = "30min"
%          opts.timepoint = "60min"
%
%    - Empty strings disable filtering for the corresponding factor.
%
%    - Filtering is performed independently in the GFP and mCherry/RFP file
%      lists before paired-channel matching.
%
% 5) GFP / mCHERRY FILE PAIRING
%    - In channelMode="both", GFP and mCherry profile files can be matched
%      using normalised filename keys.
%
%    - Enabled with:
%
%          opts.matchByName = true
%
%    - Pairing keys are generated by:
%
%          removing the leading ChannelN token
%          removing recognised treatment tokens
%          removing recognised 5min, 30min, and 60min tokens
%          removing common parent-image extensions
%          converting the remaining filename core to lowercase
%          normalising spaces and hyphens to underscores
%
%    - Duplicate normalised keys within either channel trigger an error
%      because pairing would be ambiguous.
%
%    - Only keys present in both GFP and mCherry lists are retained.
%
%    - Unmatched files are dropped with a warning.
%
%    - If no common pairing keys are found, the function stops the selected
%      analysis with an informative error.
%
%    - A pairing report is generated containing:
%
%          PairIndex
%          MatchKey
%          GFP_File
%          RFP_File
%          Treatment
%          Timepoint
%          Series
%          ROI
%
%    - If:
%
%          opts.matchByName = false
%
%      GFP and mCherry file counts must be identical.
%
%    - In that mode, the current file-list order is used rather than
%      filename-key matching.
%
% 6) SPATIAL POSITION NORMALISATION
%    - Spatial position is normalised independently for every cell.
%
%    - For each profile:
%
%          x = x / max(x)
%
%    - The spatial axis is therefore expressed relative to the maximum cell
%      profile position.
%
%    - Spatial normalisation is performed in BOTH fluorescence intensity
%      modes.
%
%    - Profiles with max(x)=0 are returned as all-NaN binned profiles.
%
% 7) FLUORESCENCE INTENSITY MODES
%    - opts.intensityMode controls y-intensity handling.
%
%    - Supported modes:
%
%          "normalized"
%
%              Per-cell min-max scaling:
%
%                  y = (y - min(y)) / (max(y) - min(y))
%
%              Each cell profile is independently scaled before spatial
%              binning.
%
%              If max(y)=min(y), the profile intensity values are replaced
%              by zeros.
%
%          "raw"
%
%              Original fluorescence intensity values are retained.
%
%              No y-intensity normalisation is performed.
%
%    - The aliases:
%
%          "normalised"
%          "norm"
%
%      are converted internally to:
%
%          "normalized"
%
% 8) SPATIAL BINNING
%    - The normalised spatial interval 0 to 1 is divided into:
%
%          opts.no_points
%
%      equally spaced bins.
%
%    - Default:
%
%          no_points = 100
%
%    - Bin edges are generated using:
%
%          linspace(0, 1, no_points + 1)
%
%    - For every cell, the mean fluorescence intensity of all measurements
%      assigned to each spatial bin is calculated.
%
%    - Each output column therefore represents one cell:
%
%          rows    = spatial bins
%          columns = cells
%
%    - Bins without valid measurements remain NaN.
%
% 9) INVALID OR EMPTY PROFILE HANDLING
%    - Non-finite x and y measurements are removed before processing.
%
%    - A CSV returns an all-NaN profile when:
%
%          importdata does not return numeric data
%          no data field is available
%          the numeric table is empty
%          fewer than two numeric columns are available
%          no finite x/y pairs remain
%          max(x) equals zero
%
%    - Missing channel files after filtering generate warnings and empty
%      output matrices rather than fabricated profiles.
%
% 10) SINGLE-CHANNEL ANALYSIS
%    - GFP-only and mCherry-only modes bypass paired-channel matching.
%
%    - Every retained file is independently converted into a binned profile.
%
%    - The pairs output preserves the source file list and records the
%      selected mode.
%
% 11) OPTIONAL PROFILE PLOTTING
%    - Enabled with:
%
%          opts.makePlot = true
%
%    - The function plots the mean spatial profile for each available
%      channel.
%
%    - When shadedErrorBar is available, SEM shading is included.
%
%    - GFP mean profiles are displayed in green.
%    - mCherry/RFP mean profiles are displayed in magenta.
%
%    - The plot contains:
%
%          normalised cell-length x-axis
%          intensity-mode-specific y-axis label
%          stage, channel mode, and intensity mode in the figure title
%          channel legend
%
%    - The current single-experiment plotting function uses a square axis,
%      background grid, and boxed axes.
%
% 12) PROFILE Y-AXIS HANDLING
%    - opts.profileYLim can manually define the profile y-axis range.
%
%    - If profileYLim is empty:
%
%          normalized mode -> y-limits [0 1]
%
%          raw mode        -> automatic limits from all finite profile
%                             values with approximately 5% padding
%
%    - For non-negative raw intensity data, the automatic lower limit is 0.
%
% 13) CONSOLE REPORTING
%    - The function reports:
%
%          number of GFP files retained in GFP-only mode
%          number of mCherry files retained in RFP-only mode
%          number of matched GFP/mCherry pairs in paired mode
%          pairing-report output paths
%          saved profile-plot paths
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% distProfileGFP : numeric matrix
%     GFP spatial intensity profiles.
%
%     Size:
%
%         [opts.no_points × nGFPCells]
%
%     Each column represents one single-cell GFP profile.
%
%     Empty when GFP is not analysed or no GFP files are retained.
%
% distProfileRFP : numeric matrix
%     mCherry/RFP spatial intensity profiles.
%
%     Size:
%
%         [opts.no_points × nRFPCells]
%
%     Each column represents one single-cell mCherry profile.
%
%     Empty when mCherry/RFP is not analysed or no files are retained.
%
% pairs : structure or structure array
%     Source-file information associated with the returned profile columns.
%
%     Fields include:
%
%         gfp
%         rfp
%         mode
%
%     In paired mode, GFP and mCherry source files correspond to matched
%     filename keys.
%
% The function can additionally generate:
%
%   - PairingReport_<stage>*<condition>.csv
%       Paired-channel filename matching report.
%
%   - Profile*<stage>*<channelMode>*<intensityMode>.png
%
%   - Profile_<stage>*<channelMode>*<intensityMode>.fig
%
% -------------------------------------------------------------------------
% DESIGN PRINCIPLES
% -------------------------------------------------------------------------
% - Performs spatial normalisation independently for every cell so profiles
%   of different physical lengths can be compared using a shared relative
%   cell-length axis.
%
% - Keeps spatial normalisation separate from fluorescence intensity
%   normalisation.
%
% - Supports raw fluorescence analysis without changing the spatial-binning
%   framework.
%
% - Uses filename-derived metadata to select experimental conditions before
%   profile processing.
%
% - Uses filename-key matching in paired mode to avoid relying solely on file
%   list order.
%
% - Returns matrices in a format that can be directly concatenated across
%   experiments by mitoIntensityFS_pooled.
%
% - Does not perform biological summary calculations or statistical tests;
%   these are handled by higher-level functions.
%
% -------------------------------------------------------------------------
% EXAMPLE USAGE
% -------------------------------------------------------------------------
%
% opts = struct();
%
% opts.channelMode = "both";
% opts.condition = "H2O2";
% opts.timepoint = "60min";
%
% opts.matchByName = true;
% opts.no_points = 100;
% opts.intensityMode = "normalized";
%
% opts.makePlot = true;
% opts.profileYLim = [];
%
% [distG, distR, pairs] = mitoIntensityFS( ...
%     'V:\CSML\Raful\Experiment1\G\Horsetail', ...
%     'V:\CSML\Raful\Experiment1\R\Horsetail', ...
%     'V:\CSML\Raful\Experiment1\ProfileOutputs', ...
%     'Horsetail', ...
%     opts);
%
% -------------------------------------------------------------------------
% Author: Raful Navarro Espindola
% Project: Mitochondrial distribution and segregation during fission yeast
%          meiosis
% -------------------------------------------------------------------------

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'no_points'),   opts.no_points   = 100; end
    if ~isfield(opts,'makePlot'),    opts.makePlot    = true; end
    if ~isfield(opts,'matchByName'), opts.matchByName = true; end
    if ~isfield(opts,'condition'),   opts.condition   = ""; end
    if ~isfield(opts,'timepoint'),   opts.timepoint   = ""; end
    if ~isfield(opts,'channelMode'), opts.channelMode = "both"; end
    if ~isfield(opts,'intensityMode') || isempty(opts.intensityMode)
        opts.intensityMode = "normalized";
    end
    if ~isfield(opts,'profileYLim'), opts.profileYLim = []; end

    stageName = char(stageName);
    no_points = opts.no_points;

    channelMode = lower(string(opts.channelMode));
    if ~ismember(channelMode, ["both","gfp","rfp"])
        error("opts.channelMode must be 'both', 'gfp', or 'rfp'.");
    end

    intensityMode = lower(string(opts.intensityMode));
    if intensityMode == "normalised" || intensityMode == "norm"
        intensityMode = "normalized";
    end
    if ~ismember(intensityMode, ["normalized","raw"])
        error("opts.intensityMode must be 'normalized' or 'raw'.");
    end
    opts.intensityMode = intensityMode;

    runGFP = any(channelMode == ["both","gfp"]);
    runRFP = any(channelMode == ["both","rfp"]);

    edges = linspace(0, 1, no_points+1);

    if runGFP
        GFPfiles = getAllFiles(pathGFP, '*.csv');
        if strlength(string(opts.condition)) > 0 || strlength(string(opts.timepoint)) > 0
            GFPfiles = filterByMeta(GFPfiles, opts.condition, opts.timepoint);
        end
    else
        GFPfiles = {};
    end

    if runRFP
        RFPfiles = getAllFiles(pathRFP, '*.csv');
        if strlength(string(opts.condition)) > 0 || strlength(string(opts.timepoint)) > 0
            RFPfiles = filterByMeta(RFPfiles, opts.condition, opts.timepoint);
        end
    else
        RFPfiles = {};
    end

    if runGFP && ~runRFP
        if isempty(GFPfiles)
            warning('No GFP CSV files found after filtering.\n  %s', pathGFP);
            distProfileGFP = [];
            distProfileRFP = [];
            pairs = struct('gfp', {{}}, 'rfp', {{}});
            return;
        end

        nCells = numel(GFPfiles);
        fprintf('Stage %s: found %d GFP cell files.\n', stageName, nCells);

        distProfileGFP = nan(no_points, nCells);
        distProfileRFP = [];

        for i = 1:nCells
            distProfileGFP(:, i) = readAndBinProfile(GFPfiles{i}, edges, intensityMode);
        end

        pairs = struct('gfp', GFPfiles, 'rfp', {{}}, 'mode', 'gfp');

        if opts.makePlot
            makeProfilePlot(distProfileGFP, [], outputDir, stageName, channelMode, no_points, opts);
        end
        return;
    end

    if runRFP && ~runGFP
        if isempty(RFPfiles)
            warning('No RFP CSV files found after filtering.\n  %s', pathRFP);
            distProfileGFP = [];
            distProfileRFP = [];
            pairs = struct('gfp', {{}}, 'rfp', {{}});
            return;
        end

        nCells = numel(RFPfiles);
        fprintf('Stage %s: found %d RFP cell files.\n', stageName, nCells);

        distProfileGFP = [];
        distProfileRFP = nan(no_points, nCells);

        for i = 1:nCells
            distProfileRFP(:, i) = readAndBinProfile(RFPfiles{i}, edges, intensityMode);
        end

        pairs = struct('gfp', {{}}, 'rfp', RFPfiles, 'mode', 'rfp');

        if opts.makePlot
            makeProfilePlot([], distProfileRFP, outputDir, stageName, channelMode, no_points, opts);
        end
        return;
    end

    if isempty(GFPfiles) || isempty(RFPfiles)
        warning(['No CSV files found after filtering in paired mode. ' ...
                 'GFP=%d files, RFP=%d files.\n  %s\n  %s'], ...
                 numel(GFPfiles), numel(RFPfiles), pathGFP, pathRFP);

        distProfileGFP = [];
        distProfileRFP = [];
        pairs = struct('gfp', {{}}, 'rfp', {{}});
        return;
    end

    if opts.matchByName
        [GFPfiles, RFPfiles, report] = matchFiles_ChannelTifCsv(GFPfiles, RFPfiles, stageName);

        if ~exist(outputDir,'dir'), mkdir(outputDir); end
        condTag = char(string(opts.condition));
        if isempty(condTag), condTag = 'ALL'; end
        repFile = fullfile(outputDir, sprintf('PairingReport_%s_%s.csv', stageName, condTag));
        writetable(report, repFile);
        fprintf('Pairing report saved: %s\n', repFile);
    else
        if numel(GFPfiles) ~= numel(RFPfiles)
            error(['GFP and RFP file counts differ (GFP=%d, RFP=%d). ' ...
                   'Set opts.matchByName=true.'], numel(GFPfiles), numel(RFPfiles));
        end
    end

    nCells = numel(GFPfiles);
    if nCells == 0
        warning('No matched GFP/RFP file pairs found for stage "%s".', stageName);
        distProfileGFP = [];
        distProfileRFP = [];
        pairs = struct('gfp', {{}}, 'rfp', {{}});
        return;
    end

    fprintf('Stage %s: matched %d GFP/RFP cell files.\n', stageName, nCells);

    distProfileGFP = nan(no_points, nCells);
    distProfileRFP = nan(no_points, nCells);

    for i = 1:nCells
        distProfileGFP(:, i) = readAndBinProfile(GFPfiles{i}, edges, intensityMode);
        distProfileRFP(:, i) = readAndBinProfile(RFPfiles{i}, edges, intensityMode);
    end

    pairs = struct('gfp', GFPfiles, 'rfp', RFPfiles, 'mode', 'both');

    if opts.makePlot
        makeProfilePlot(distProfileGFP, distProfileRFP, outputDir, stageName, channelMode, no_points, opts);
    end
end

% ---------------- Helpers ----------------

function makeProfilePlot(distProfileGFP, distProfileRFP, outputDir, stageName, channelMode, no_points, opts)

    runGFP = ~isempty(distProfileGFP);
    runRFP = ~isempty(distProfileRFP);

    f = figure('Color','w');
    hold on;

    legendLabels = {};
    legendHandles = gobjects(0);

    if runGFP
        meanG = mean(distProfileGFP, 2, 'omitnan');
        nG = sum(~all(isnan(distProfileGFP),1));
        semG = std(distProfileGFP, 0, 2, 'omitnan') ./ sqrt(max(nG,1));
        if nG == 0, semG = nan(size(meanG)); end

        if exist('shadedErrorBar','file') == 2
            shadedErrorBar(1:no_points, meanG, semG, struct('color',[0.6 1 0.6],'LineWidth',1), 1);
            h = plot(1:no_points, meanG, 'Color',[0 0.6 0], 'LineWidth',2);
        else
            h = plot(1:no_points, meanG, 'Color',[0 0.6 0], 'LineWidth',2);
        end
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'GFP'; %#ok<AGROW>
    end

    if runRFP
        meanR = mean(distProfileRFP, 2, 'omitnan');
        nR = sum(~all(isnan(distProfileRFP),1));
        semR = std(distProfileRFP, 0, 2, 'omitnan') ./ sqrt(max(nR,1));
        if nR == 0, semR = nan(size(meanR)); end

        if exist('shadedErrorBar','file') == 2
            shadedErrorBar(1:no_points, meanR, semR, struct('color',[1 0.6 0.8],'LineWidth',1), 1);
            h = plot(1:no_points, meanR, 'Color',[0.8 0 0.8], 'LineWidth',2);
        else
            h = plot(1:no_points, meanR, 'Color',[0.8 0 0.8], 'LineWidth',2);
        end
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = 'RFP'; %#ok<AGROW>
    end

    [yLabel, modeTag] = intensityModeLabel(opts.intensityMode);
    yL = profileYLimits({distProfileGFP, distProfileRFP}, opts);

    axis square;
    xlim([1 no_points]);
    ylim(yL);
    title(sprintf('Meiocytes - %s (%s, %s)', ...
        stageName, upper(char(channelMode)), upper(char(modeTag))), ...
        'Interpreter','none');
    xlabel('Normalized Cell Length (binned)');
    ylabel(yLabel);
    grid on; box on;

    if ~isempty(legendHandles)
        legend(legendHandles, legendLabels, 'Location','best');
    end

    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end

    baseName = strrep(stageName, ' ', '_');

    pngName = fullfile(outputDir, sprintf('Profile_%s_%s_%s.png', ...
        baseName, char(channelMode), char(modeTag)));
    figName = fullfile(outputDir, sprintf('Profile_%s_%s_%s.fig', ...
        baseName, char(channelMode), char(modeTag)));
    saveas(f, pngName);
    saveas(f, figName);
    close(f);

    fprintf('Profile plot saved: %s\n', pngName);
end

function prof = readAndBinProfile(csvFile, edges, intensityMode)
    D = importdata(csvFile);
    if ~isstruct(D) || ~isfield(D,'data') || isempty(D.data) || size(D.data,2) < 2
        prof = nan(numel(edges)-1, 1);
        return;
    end
    data = D.data;

    x = data(:,1);
    y = data(:,2);

    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);

    if isempty(x) || isempty(y)
        prof = nan(numel(edges)-1, 1);
        return;
    end

    xMax = max(x);
    if xMax == 0
        prof = nan(numel(edges)-1, 1);
        return;
    end
    x = x / xMax;

    if string(intensityMode) == "normalized"
        yMin = min(y);
        yMax = max(y);
        if yMax == yMin
            y = zeros(size(y));
        else
            y = (y - yMin) / (yMax - yMin);
        end
    elseif string(intensityMode) == "raw"
        % Keep the original fluorescence intensity values.
        % x is still normalized to cell length so profiles can be binned.
    else
        error("Unknown intensityMode: %s", char(string(intensityMode)));
    end

    bin = discretize(x, edges);

    nBins = numel(edges)-1;
    prof = nan(nBins,1);
    valid = ~isnan(bin);
    if any(valid)
        sums = accumarray(bin(valid), y(valid), [nBins 1], @sum, NaN);
        cnts = accumarray(bin(valid), 1,       [nBins 1], @sum, NaN);
        prof = sums ./ cnts;
    end
end


function [yLabel, modeTag] = intensityModeLabel(intensityMode)
    intensityMode = lower(string(intensityMode));
    if intensityMode == "raw"
        yLabel = 'Raw intensity (a.u.)';
        modeTag = "raw";
    else
        yLabel = 'Normalized Intensity ((x-min)/range)';
        modeTag = "normalized";
    end
end

function yL = profileYLimits(profileCells, opts)
    intensityMode = lower(string(opts.intensityMode));

    if isfield(opts,'profileYLim') && ~isempty(opts.profileYLim)
        yL = opts.profileYLim;
        return;
    end

    if intensityMode ~= "raw"
        yL = [0 1];
        return;
    end

    vals = [];
    for k = 1:numel(profileCells)
        X = profileCells{k};
        if ~isempty(X)
            vals = [vals; X(:)]; %#ok<AGROW>
        end
    end

    vals = vals(isfinite(vals));

    if isempty(vals)
        yL = [0 1];
        return;
    end

    yMin = min(vals);
    yMax = max(vals);

    if yMax == yMin
        pad = max(abs(yMax) * 0.05, 1);
        yL = [yMin - pad, yMax + pad];
    else
        pad = 0.05 * (yMax - yMin);
        yL = [yMin - pad, yMax + pad];
    end

    if yMin >= 0
        yL(1) = 0;
    end
end

function [GFPout, RFPout, report] = matchFiles_ChannelTifCsv(GFPfiles, RFPfiles, stageName)
    gKey = cellfun(@pairKey, GFPfiles, 'UniformOutput', false);
    rKey = cellfun(@pairKey, RFPfiles, 'UniformOutput', false);

    if numel(unique(gKey)) ~= numel(gKey)
        error('Stage %s: Duplicate GFP keys after normalization. Pairing would be ambiguous.', stageName);
    end
    if numel(unique(rKey)) ~= numel(rKey)
        error('Stage %s: Duplicate RFP keys after normalization. Pairing would be ambiguous.', stageName);
    end

    common = intersect(gKey, rKey);
    if isempty(common)
        egG = gKey{1}; egR = rKey{1};
        error(['Stage %s: No matching pairs after normalization.\n' ...
               'Example normalized GFP key: "%s"\nExample normalized RFP key: "%s"\n' ...
               'Check that filenames share the same <core> after ChannelX_ and parent extension removal.\n'], ...
              stageName, egG, egR);
    end

    common = sort(common);

    gMap = containers.Map(gKey, GFPfiles);
    rMap = containers.Map(rKey, RFPfiles);

    GFPout = cell(size(common));
    RFPout = cell(size(common));
    for i = 1:numel(common)
        GFPout{i} = gMap(common{i});
        RFPout{i} = rMap(common{i});
    end

    gMeta = cellfun(@parseMetaFromCsvPath, GFPout, 'UniformOutput', false);

    gSeries = cellfun(@(m)m.series, gMeta);
    gRoi    = cellfun(@(m)m.roi, gMeta);
    gTp     = string(cellfun(@(m)char(m.timepoint), gMeta, 'UniformOutput', false));
    gTr     = string(cellfun(@(m)char(m.treatment), gMeta, 'UniformOutput', false));

    report = table((1:numel(common))', string(common(:)), string(GFPout(:)), string(RFPout(:)), ...
        gTr(:), gTp(:), gSeries(:), gRoi(:), ...
        'VariableNames', {'PairIndex','MatchKey','GFP_File','RFP_File','Treatment','Timepoint','Series','ROI'});

    if numel(common) < numel(gKey) || numel(common) < numel(rKey)
        warning('Stage %s: Matched %d pairs, but GFP had %d files and RFP had %d files. Unmatched files were dropped.', ...
            stageName, numel(common), numel(gKey), numel(rKey));
    end
end

function key = pairKey(fullPath)
    meta = parseMetaFromCsvPath(fullPath);
    key = char(meta.keyCore);
end

function filesOut = filterByMeta(filesIn, wantedTreatment, wantedTimepoint)

    wantedTreatment = string(wantedTreatment);
    wantedTimepoint = string(wantedTimepoint);

    keep = false(size(filesIn));

    for i = 1:numel(filesIn)
        meta = parseMetaFromCsvPath(filesIn{i});

        ok = true;
        if wantedTreatment ~= ""
            ok = ok && (meta.treatment == wantedTreatment);
        end
        if wantedTimepoint ~= ""
            ok = ok && (meta.timepoint == wantedTimepoint);
        end

        keep(i) = ok;
    end

    filesOut = filesIn(keep);
end

function meta = parseMetaFromCsvPath(fullPath)
% parseMetaFromCsvPath
%
% Parses experimental metadata from a single-cell profile CSV filename.
%
% Extracts:
%   channel     : ChannelN prefix, when present
%   treatment   : "Veh", "H2O2", or "UNKNOWN"
%   timepoint   : arbitrary minute-based timepoint (e.g. "5min", "60min",
%                 "180min") or "UNKNOWN"
%   series      : seriesN identifier
%   roi         : roiN identifier
%   keyCore     : normalized filename core used for GFP/mCherry pairing
%
% Timepoint handling:
%   - Recognises arbitrary Nmin tokens:
%         5min
%         30min
%         60min
%         180min
%
%   - Also recognises an optional trailing T:
%         180minT
%         60minT
%
%     These are assigned to the same numeric timepoint:
%         180minT -> "180min"
%
%   - Recognises exact T-number tokens:
%         T0
%         T5
%         T30
%         T60
%
%   - T0 and 0min are treated as the 5 min baseline for compatibility with
%     the current experimental nomenclature.
%
%   - If no recognised timepoint token is present, the timepoint remains
%     "UNKNOWN". Unknown filenames are NOT silently assigned to 5 min.

```
meta = struct( ...
    'channel',NaN, ...
    'treatment',"UNKNOWN", ...
    'timepoint',"UNKNOWN", ...
    'series',NaN, ...
    'roi',NaN, ...
    'keyCore',"");

[~, nameNoCsv, ~] = fileparts(fullPath);
s = string(nameNoCsv);

% Remove a parent-image extension if retained in the CSV filename
s = regexprep(s, ...
    '\.(tif|tiff|nd2|czi|lif|ome\.tif)$', ...
    '', ...
    'ignorecase');

% Split filename into underscore-delimited tokens
tok = split(s, "_");
tok = tok(tok ~= "");
lowTok = lower(tok);

% -------------------- Channel --------------------
if numel(tok) >= 1
    m = regexp(tok(1), ...
        '^Channel(\d+)$', ...
        'tokens', ...
        'once', ...
        'ignorecase');

    if ~isempty(m)
        meta.channel = str2double(m{1});
    end
end

% -------------------- Series and ROI --------------------
for k = 1:numel(tok)

    ms = regexp(tok(k), ...
        '^series(\d+)$', ...
        'tokens', ...
        'once', ...
        'ignorecase');

    if ~isempty(ms)
        meta.series = str2double(ms{1});
    end

    mr = regexp(tok(k), ...
        '^roi(\d+)$', ...
        'tokens', ...
        'once', ...
        'ignorecase');

    if ~isempty(mr)
        meta.roi = str2double(mr{1});
    end
end

% -------------------- Treatment --------------------
if any(strcmp(lowTok, "h2o2")) || any(contains(lowTok, "h2o2"))
    meta.treatment = "H2O2";

elseif any(startsWith(lowTok, "veh")) || ...
        any(startsWith(lowTok, "vehicle"))
    meta.treatment = "Veh";

else
    meta.treatment = "UNKNOWN";
end

% -------------------- Timepoint --------------------
% Prefer explicit Nmin / NminT tokens.
%
% Examples:
%   5min     -> 5min
%   60min    -> 60min
%   180min   -> 180min
%   180minT  -> 180min

for k = 1:numel(tok)

    timeTok = regexp(tok(k), ...
        '^(\d+)min(?:t)?$', ...
        'tokens', ...
        'once', ...
        'ignorecase');

    if ~isempty(timeTok)

        timeValue = str2double(timeTok{1});

        % Current pipeline convention:
        % explicit 0 min represents the 5 min baseline.
        if timeValue == 0
            timeValue = 5;
        end

        meta.timepoint = string(timeValue) + "min";
        break;
    end
end

% Fallback: exact T-number token.
%
% Examples:
%   T0  -> 5min
%   T5  -> 5min
%   T30 -> 30min
%   T60 -> 60min

if meta.timepoint == "UNKNOWN"

    for k = 1:numel(tok)

        timeTok = regexp(tok(k), ...
            '^T(\d+)$', ...
            'tokens', ...
            'once', ...
            'ignorecase');

        if ~isempty(timeTok)

            timeValue = str2double(timeTok{1});

            % Preserve baseline convention
            if timeValue == 0
                timeValue = 5;
            end

            meta.timepoint = string(timeValue) + "min";
            break;
        end
    end
end

% -------------------- Pairing key --------------------
% Remove channel, treatment, and timepoint tokens so matched GFP and
% mCherry files reduce to the same cell-specific filename core.

keep = true(size(tok));

% Remove leading ChannelN token
if numel(tok) >= 1 && startsWith(lower(tok(1)), "channel")
    keep(1) = false;
end

for k = 1:numel(tok)

    t = lower(tok(k));

    % Remove treatment token
    if startsWith(t, "veh") || ...
            startsWith(t, "vehicle") || ...
            t == "h2o2"

        keep(k) = false;
        continue;
    end

    % Remove arbitrary Nmin / NminT timepoint token
    isMinuteToken = ~isempty(regexp(t, ...
        '^\d+min(?:t)?$', ...
        'once'));

    % Remove exact T-number timepoint token
    isTToken = ~isempty(regexp(t, ...
        '^t\d+$', ...
        'once'));

    if isMinuteToken || isTToken
        keep(k) = false;
    end
end

core = lower(strjoin(tok(keep), "_"));
core = regexprep(core, '[\s\-]+', '_');

meta.keyCore = string(core);
```

end