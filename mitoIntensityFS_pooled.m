function [distProfileGFP, distProfileRFP, pooledMeta] = mitoIntensityFS_pooled(baseDirs, outputDir, stageName, opts)
% mitoIntensityFS_pooled
%
% Multi-experiment pooling layer for binned single-cell GFP and/or mCherry
% spatial intensity profiles.
%
% -------------------------------------------------------------------------
% OVERVIEW
% -------------------------------------------------------------------------
% This function processes one selected meiotic stage, treatment condition,
% and timepoint across multiple experiment directories.
%
% For each experiment, the function identifies the corresponding GFP and/or
% mCherry stage folders and calls mitoIntensityFS to:
%
%       - discover single-cell profile CSV files
%       - filter profiles by condition and timepoint
%       - pair GFP and mCherry files when required
%       - normalise cell-length position
%       - optionally normalise fluorescence intensity
%       - bin each single-cell profile
%
% The resulting profile matrices are concatenated across experiments while
% preserving one row of metadata for each pooled cell/profile.
%
% The function can optionally generate per-experiment profile plots through
% mitoIntensityFS and a pooled mean ± SEM profile plot across all experiments.
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
% This function is the multi-experiment pooling layer.
%
% It does NOT calculate front/back polarisation or the Regional Segregation
% Score. Those analyses are performed by runPooledPol after profiles have
% been pooled.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% baseDirs : cell array or string array
%     Experiment root directories.
%
%     Each experiment is expected to contain:
%
%         <baseDir>\G<stageName>*.csv
%         <baseDir>\R<stageName>*.csv
%
%     depending on opts.channelMode.
%
% outputDir : character vector or string-compatible path
%     Directory where pooled metadata, pooled profile plots, and
%     per-experiment report subdirectories are saved.
%
% stageName : character vector or string-compatible stage name
%     Stage folder to analyse.
%
%     Common values:
%
%         "Horsetail"
%         "M1"
%         "M2"
%         "M2l"
%         "Spores"
%
% opts : structure containing profile-processing and output options.
%
% -------------------------------------------------------------------------
% REQUIRED ON MATLAB PATH
% -------------------------------------------------------------------------
%
%       mitoIntensityFS.m
%       getAllFiles.m
%
% shadedErrorBar.m is optional for per-experiment profile plots generated
% by mitoIntensityFS.
%
% -------------------------------------------------------------------------
% CORE FEATURES
% -------------------------------------------------------------------------
%
% 1) MULTI-EXPERIMENT PROFILE POOLING
%    - Iterates through every directory supplied in baseDirs.
%
%    - Calls mitoIntensityFS once per experiment for the selected:
%
%          Stage × Condition × Timepoint
%
%    - Profile matrices returned from each experiment are concatenated
%      column-wise.
%
%    - Final matrix structure:
%
%          rows    = spatial profile bins
%          columns = pooled single-cell profiles
%
% 2) CHANNEL-AWARE FOLDER HANDLING
%    - opts.channelMode controls which folders are required:
%
%          "both"
%              Requires:
%                  <baseDir>\G<stageName>
%                  <baseDir>\R<stageName>
%
%          "gfp"
%              Requires only:
%                  <baseDir>\G<stageName>
%
%          "rfp"
%              Requires only:
%                  <baseDir>\R<stageName>
%
%    - Missing required channel folders cause the affected experiment to be
%      skipped with a warning.
%
%    - Unused channel paths are replaced with empty placeholders before
%      mitoIntensityFS is called.
%
% 3) CONDITION AND TIMEPOINT FILTERING
%    - opts.condition and opts.timepoint are passed directly to
%      mitoIntensityFS.
%
%    - Common condition values:
%
%          "Veh"
%          "H2O2"
%          ""
%
%    - Common timepoint values:
%
%          "5min"
%          "30min"
%          "60min"
%          ""
%
%    - An empty condition or timepoint disables filtering for that factor.
%
%    - Filename metadata parsing and filtering are performed by
%      mitoIntensityFS rather than by this pooling function.
%
% 4) NORMALISED OR RAW INTENSITY MODE
%    - opts.intensityMode is passed to mitoIntensityFS.
%
%    - Supported values:
%
%          "normalized"
%              Per-cell min-max fluorescence scaling before binning.
%
%          "raw"
%              Original fluorescence intensity values are retained.
%
%    - The aliases "normalised" and "norm" are converted internally to
%      "normalized".
%
%    - Cell-length position remains normalised in both intensity modes.
%
% 5) CONTROL OF PROFILE RESOLUTION
%    - opts.no_points defines the number of spatial bins in each profile.
%
%    - Default:
%
%          no_points = 100
%
%    - Every returned profile matrix therefore has:
%
%          opts.no_points rows
%
% 6) PER-EXPERIMENT PROCESSING
%    - Each experiment is assigned an experiment name using the final folder
%      name in its baseDir path.
%
%    - If the final folder name cannot be determined, the fallback label is:
%
%          ExpN
%
%    - A dedicated report directory is created for every experiment:
%
%          <outputDir>\PerExperimentReports<ExperimentName>
%
%    - mitoIntensityFS writes paired-channel matching reports to this
%      directory when channelMode="both" and matchByName=true.
%
%    - opts.perExperimentPlots controls whether mitoIntensityFS also creates
%      an individual profile figure for each experiment.
%
% 7) PROFILE CONCATENATION
%    - GFP profiles are concatenated into distProfileGFP.
%
%    - mCherry/RFP profiles are concatenated into distProfileRFP.
%
%    - In single-channel mode, the unused output profile matrix remains
%      empty.
%
%    - The number of new cells contributed by each experiment is determined
%      from the returned profile matrices.
%
%    - Experiments returning no cells after filtering are skipped.
%
% 8) POOLED CELL PROVENANCE
%    - One metadata row is generated for every pooled cell/profile.
%
%    - pooledMeta contains:
%
%          ExperimentPath
%          ExperimentName
%          Stage
%          Condition
%          Timepoint
%          CellIndexWithinExperiment
%          GFP_File
%          RFP_File
%          ChannelMode
%          IntensityMode
%
%    - ExperimentPath stores the full experiment root directory.
%
%    - ExperimentName stores the final experiment folder name.
%
%    - CellIndexWithinExperiment numbers pooled cells independently within
%      each experiment.
%
%    - GFP_File and RFP_File preserve the source CSV file path when
%      available.
%
%    - These metadata allow downstream analyses such as runPooledPol to
%      recover experiment identity after profiles have been concatenated.
%
% 9) OPTIONAL POOLED METADATA EXPORT
%    - Enabled with:
%
%          opts.saveMeta = true
%
%    - The pooled metadata table is saved as:
%
%          PooledMeta_<stage_condition_timepoint_channel_mode_intensity_mode>.csv
%
%    - The exact filename tag is assembled from stageName and the current
%      profile-processing options.
%
% 10) OPTIONAL POOLED PROFILE PLOT
%    - Enabled with:
%
%          opts.makePlot = true
%
%    - Generates a pooled profile figure using all concatenated cells.
%
%    - For each available channel, the plot displays:
%
%          pooled mean spatial profile
%          SEM shading
%
%    - GFP is displayed in green.
%    - mCherry is displayed in magenta.
%
%    - The legend reports the number of valid profile columns contributing
%      to each channel.
%
%    - The figure uses an A4 landscape-friendly layout with:
%
%          no background grid
%          no top or right box axes
%          outward ticks
%          enlarged axis labels and numeric tick labels
%
%    - The mCherry legend label is displayed as:
%
%          mCherry
%
%      rather than RFP.
%
% 11) PROFILE Y-AXIS HANDLING
%    - opts.profileYLim can manually define a common profile y-axis range.
%
%    - If profileYLim is empty:
%
%          normalized mode -> y-limits [0 1]
%
%          raw mode        -> automatic limits from all finite pooled
%                             profile values with approximately 5% padding
%
%    - For non-negative raw data, the automatic lower y-limit is set to 0.
%
% 12) VERBOSE STATUS REPORTING
%    - opts.verbose controls detailed console reporting.
%
%    - The function can report:
%
%          experiment currently being pooled
%          number of cells added from each experiment
%          combinations returning no cells
%          pooled metadata output path
%          final GFP and mCherry pooled cell counts
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% distProfileGFP : numeric matrix
%     Pooled GFP spatial profiles.
%
%     Size:
%
%         [opts.no_points × nGFPCells]
%
%     Empty when GFP is not analysed or no GFP profiles are returned.
%
% distProfileRFP : numeric matrix
%     Pooled mCherry/RFP spatial profiles.
%
%     Size:
%
%         [opts.no_points × nRFPCells]
%
%     Empty when mCherry/RFP is not analysed or no profiles are returned.
%
% pooledMeta : table
%     One metadata row per pooled cell/profile containing experiment,
%     stage, condition, timepoint, source-file, channel-mode, and
%     intensity-mode information.
%
% The function can additionally generate:
%
%   - PooledMeta_<tag>.csv
%       Pooled per-cell metadata.
%
%   - Pooled profile PNG and FIG files.
%
%   - PerExperimentReports<ExperimentName>\PairingReport_*.csv
%       GFP/mCherry matching reports generated by mitoIntensityFS.
%
%   - Optional per-experiment profile PNG and FIG files when
%       opts.perExperimentPlots=true.
%
% -------------------------------------------------------------------------
% DESIGN PRINCIPLES
% -------------------------------------------------------------------------
% - Uses mitoIntensityFS as the single source of truth for profile discovery,
%   metadata filtering, channel pairing, and profile binning.
%
% - The pooling layer concatenates profiles without reprocessing their
%   spatial or fluorescence intensity values.
%
% - Preserves full experiment paths so downstream analyses can distinguish
%   experiments even when final folder names are repeated.
%
% - Keeps pooled cell provenance aligned with the concatenated profile
%   matrices.
%
% - Supports paired and single-channel experiments through the same pooling
%   function.
%
% - Separates per-experiment plotting from pooled plotting.
%
% -------------------------------------------------------------------------
% EXAMPLE USAGE
% -------------------------------------------------------------------------
%
% baseDirs = {
%     'V:\CSML\Raful\Exp1'
%     'V:\CSML\Raful\Exp2'
%     'V:\CSML\Raful\Exp3'
%     'V:\CSML\Raful\Exp4'
% };
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
% opts.perExperimentPlots = false;
% opts.saveMeta = true;
% opts.verbose = true;
%
% opts.profileYLim = [];
%
% [poolG, poolR, poolMeta] = mitoIntensityFS_pooled( ...
%     baseDirs, ...
%     'V:\CSML\Raful\PooledAnalysis', ...
%     'Horsetail', ...
%     opts);
%
% -------------------------------------------------------------------------
% Author: Raful Navarro Espindola
% Project: Mitochondrial distribution and segregation during fission yeast
%          meiosis
% -------------------------------------------------------------------------

    % -------------------- defaults --------------------
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'channelMode'),         opts.channelMode = "both"; end
    if ~isfield(opts,'condition'),           opts.condition = ""; end
    if ~isfield(opts,'timepoint'),           opts.timepoint = ""; end
    if ~isfield(opts,'matchByName'),         opts.matchByName = true; end
    if ~isfield(opts,'makePlot'),            opts.makePlot = true; end
    if ~isfield(opts,'no_points'),           opts.no_points = 100; end
    if ~isfield(opts,'saveMeta'),            opts.saveMeta = true; end
    if ~isfield(opts,'verbose'),             opts.verbose = true; end
    if ~isfield(opts,'perExperimentPlots'),  opts.perExperimentPlots = false; end
    if ~isfield(opts,'intensityMode') || isempty(opts.intensityMode)
        opts.intensityMode = "normalized";
    end
    if ~isfield(opts,'profileYLim'), opts.profileYLim = []; end

    baseDirs = cellstr(string(baseDirs(:)));
    stageName = char(string(stageName));
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

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % -------------------- storage --------------------
    distProfileGFP = [];
    distProfileRFP = [];

    metaExperiment = strings(0,1);
    metaExperimentName = strings(0,1);
    metaStage = strings(0,1);
    metaCondition = strings(0,1);
    metaTimepoint = strings(0,1);
    metaCellIndexWithinExperiment = zeros(0,1);
    metaGFPFile = strings(0,1);
    metaRFPFile = strings(0,1);
    metaMode = strings(0,1);
    metaIntensityMode = strings(0,1);

    % Use a modified opts for mitoIntensityFS so we can suppress
    % per-experiment profile plots if desired.
    optsSingle = opts;
    optsSingle.makePlot = logical(opts.perExperimentPlots);

    % -------------------- loop over experiments --------------------
    for e = 1:numel(baseDirs)
        thisBase = baseDirs{e};

        % experiment label = folder name
        [~, expName, ~] = fileparts(thisBase);
        if expName == ""
            expName = sprintf("Exp%d", e);
        end

        if opts.verbose
            fprintf('\n=== Pooling experiment %d/%d: %s ===\n', e, numel(baseDirs), thisBase);
        end

        pathG = fullfile(thisBase, 'G', stageName);
        pathR = fullfile(thisBase, 'R', stageName);

        if runGFP && ~isfolder(pathG)
            warning('Skipping experiment "%s": missing GFP folder:\n  %s', expName, pathG);
            continue;
        end
        if runRFP && ~isfolder(pathR)
            warning('Skipping experiment "%s": missing RFP folder:\n  %s', expName, pathR);
            continue;
        end

        % Safe placeholders for unused paths in single-channel mode
        if ~runGFP, pathG = ''; end
        if ~runRFP, pathR = ''; end

        [~, expName, ~] = fileparts(thisBase);
        if expName == ""
            expName = sprintf("Exp%d", e);
        end

        expOutDir = fullfile(outputDir, 'PerExperimentReports', char(expName));
        if ~exist(expOutDir, 'dir')
            mkdir(expOutDir);
        end

        [distG, distR, pairs] = mitoIntensityFS(pathG, pathR, expOutDir, stageName, optsSingle);

        % Determine how many new cells were returned
        nNewG = size(distG, 2);
        nNewR = size(distR, 2);
        nNew = max(nNewG, nNewR);

        if nNew == 0
            if opts.verbose
                fprintf('  No cells returned for %s | %s | %s in %s\n', ...
                    string(stageName), string(opts.condition), string(opts.timepoint), expName);
            end
            continue;
        end

        % Concatenate pooled profiles
        if runGFP
            if isempty(distProfileGFP)
                distProfileGFP = distG;
            else
                distProfileGFP = [distProfileGFP, distG]; %#ok<AGROW>
            end
        end

        if runRFP
            if isempty(distProfileRFP)
                distProfileRFP = distR;
            else
                distProfileRFP = [distProfileRFP, distR]; %#ok<AGROW>
            end
        end

        % Collect per-cell file info
        gfpFiles = strings(nNew,1);
        rfpFiles = strings(nNew,1);

        if isstruct(pairs) && isfield(pairs,'gfp')
            if numel(pairs) == 1
                if ~isempty(pairs.gfp)
                    tmp = string(reshape(pairs.gfp(:), [], 1));
                    gfpFiles(1:min(numel(tmp),nNew)) = tmp(1:min(numel(tmp),nNew));
                end
            else
                tmp = string({pairs.gfp}');
                gfpFiles(1:min(numel(tmp),nNew)) = tmp(1:min(numel(tmp),nNew));
            end
        end

        if isstruct(pairs) && isfield(pairs,'rfp')
            if numel(pairs) == 1
                if ~isempty(pairs.rfp)
                    tmp = string(reshape(pairs.rfp(:), [], 1));
                    rfpFiles(1:min(numel(tmp),nNew)) = tmp(1:min(numel(tmp),nNew));
                end
            else
                tmp = string({pairs.rfp}');
                rfpFiles(1:min(numel(tmp),nNew)) = tmp(1:min(numel(tmp),nNew));
            end
        end

        metaExperiment = [metaExperiment; repmat(string(thisBase), nNew, 1)]; %#ok<AGROW>
        metaExperimentName = [metaExperimentName; repmat(string(expName), nNew, 1)]; %#ok<AGROW>
        metaStage = [metaStage; repmat(string(stageName), nNew, 1)]; %#ok<AGROW>
        metaCondition = [metaCondition; repmat(string(opts.condition), nNew, 1)]; %#ok<AGROW>
        metaTimepoint = [metaTimepoint; repmat(string(opts.timepoint), nNew, 1)]; %#ok<AGROW>
        metaCellIndexWithinExperiment = [metaCellIndexWithinExperiment; (1:nNew)']; %#ok<AGROW>
        metaGFPFile = [metaGFPFile; gfpFiles]; %#ok<AGROW>
        metaRFPFile = [metaRFPFile; rfpFiles]; %#ok<AGROW>
        metaMode = [metaMode; repmat(channelMode, nNew, 1)]; %#ok<AGROW>
        metaIntensityMode = [metaIntensityMode; repmat(intensityMode, nNew, 1)]; %#ok<AGROW>

        if opts.verbose
            fprintf('  Added %d pooled cells from %s\n', nNew, expName);
        end
    end

    % -------------------- metadata table --------------------
    pooledMeta = table( ...
        metaExperiment, ...
        metaExperimentName, ...
        metaStage, ...
        metaCondition, ...
        metaTimepoint, ...
        metaCellIndexWithinExperiment, ...
        metaGFPFile, ...
        metaRFPFile, ...
        metaMode, ...
        metaIntensityMode, ...
        'VariableNames', { ...
            'ExperimentPath', ...
            'ExperimentName', ...
            'Stage', ...
            'Condition', ...
            'Timepoint', ...
            'CellIndexWithinExperiment', ...
            'GFP_File', ...
            'RFP_File', ...
            'ChannelMode', ...
            'IntensityMode'});

    % -------------------- save metadata --------------------
    if opts.saveMeta && ~isempty(pooledMeta)
        metaFile = fullfile(outputDir, sprintf('PooledMeta_%s.csv', makeTag(stageName, opts)));
        writetable(pooledMeta, metaFile);

        if opts.verbose
            fprintf('Pooled metadata saved: %s\n', metaFile);
        end
    end

    % -------------------- pooled plot --------------------
    if opts.makePlot
        makePooledProfilePlot(distProfileGFP, distProfileRFP, outputDir, stageName, opts);
    end

    % -------------------- status --------------------
    if opts.verbose
        nG = size(distProfileGFP, 2);
        nR = size(distProfileRFP, 2);
        fprintf('\nPooled summary for %s | %s | %s:\n', ...
            string(stageName), string(opts.condition), string(opts.timepoint));
        fprintf('  GFP cells: %d\n', nG);
        fprintf('  RFP cells: %d\n', nR);
    end
end

% ========================================================================
function makePooledProfilePlot(distProfileGFP, distProfileRFP, outputDir, stageName, opts)

    no_points = opts.no_points;
    channelMode = lower(string(opts.channelMode));

    runGFP = ~isempty(distProfileGFP);
    runRFP = ~isempty(distProfileRFP);

    if ~runGFP && ~runRFP %#ok<*BDSCA>
        warning('No pooled profiles available to plot for stage "%s".', stageName);
        return;
    end

    % A4 landscape figure
    f = figure('Color','w', ...
        'Units','centimeters', ...
        'Position',[2 2 29.7 21]);

    ax = axes(f);
    hold(ax, 'on');

    % A4 export settings
    set(f, 'PaperUnits', 'centimeters');
    set(f, 'PaperSize', [29.7 21]);
    set(f, 'PaperPosition', [0 0 29.7 21]);
    set(f, 'PaperPositionMode', 'manual');

    legendLabels = {};
    legendHandles = gobjects(0);

    x = 1:no_points;

    if runGFP
        meanG = mean(distProfileGFP, 2, 'omitnan');
        nG = sum(~all(isnan(distProfileGFP),1));
        semG = std(distProfileGFP, 0, 2, 'omitnan') ./ sqrt(max(nG,1));
        if nG == 0, semG = nan(size(meanG)); end

        xv = [x, fliplr(x)];
        yv = [transpose(meanG - semG), fliplr(transpose(meanG + semG))];
        fill(xv, yv, [0.85 1.00 0.85], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.45);

        h = plot(x, meanG, 'Color', [0.00 0.60 0.00], 'LineWidth', 2);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = sprintf('GFP (n=%d)', nG); %#ok<AGROW>
    end

    if runRFP
        meanR = mean(distProfileRFP, 2, 'omitnan');
        nR = sum(~all(isnan(distProfileRFP),1));
        semR = std(distProfileRFP, 0, 2, 'omitnan') ./ sqrt(max(nR,1));
        if nR == 0, semR = nan(size(meanR)); end

        xv = [x, fliplr(x)];
        yv = [transpose(meanR - semR), fliplr(transpose(meanR + semR))];
        fill(xv, yv, [1.00 0.88 0.95], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.45);

        h = plot(x, meanR, 'Color', [0.80 0.00 0.80], 'LineWidth', 2);
        legendHandles(end+1) = h; %#ok<AGROW>
        legendLabels{end+1} = sprintf('mCherry (n=%d)', nR); %#ok<AGROW>
    end

    [yLabel, modeTag] = intensityModeLabel(opts.intensityMode);
    yL = profileYLimits({distProfileGFP, distProfileRFP}, opts);

    % Do not use axis square here; A4 landscape gives the profile more room
    xlim([1 no_points]);
    ylim(yL);

    % Remove background grid, top x-axis, and right y-axis
    grid off;
    box off;

    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;

    % Numeric axis labels
    ax.FontSize = 20;
    ax.FontWeight = 'normal';

    % A4-friendly plotting area
    ax.Position = [0.12 0.16 0.76 0.76];

    % Axis labels: 26 pt, bold
    xl = xlabel('Normalized Cell Length (binned)');
    yl = ylabel(yLabel);

    xl.FontSize = 26;
    xl.FontWeight = 'bold';

    yl.FontSize = 26;
    yl.FontWeight = 'bold';

    % No figure title
    title('');

    % Keep legend, but make it cleaner
    if ~isempty(legendHandles)
        legend(legendHandles, legendLabels, ...
            'Location', 'northeast', ...
            'Box', 'off', ...
            'FontSize', 18);
    end

    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end

    stem = makeTag(stageName, opts);
    pngName = fullfile(outputDir, sprintf('PooledProfile_%s.png', stem));
    figName = fullfile(outputDir, sprintf('PooledProfile_%s.fig', stem));

    exportgraphics(f, pngName, 'Resolution', 300);
    savefig(f, figName);
    close(f);

    fprintf('Pooled profile plot saved: %s\n', pngName);
end

% ========================================================================

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

function stem = makeTag(stageName, opts)
    cond = string(opts.condition);
    tp   = string(opts.timepoint);
    mode = string(opts.channelMode);

    if isfield(opts,'intensityMode')
        intMode = string(opts.intensityMode);
    else
        intMode = "normalized";
    end

    if cond == "", cond = "ALL"; end
    if tp == "",   tp = "ALL";   end
    if mode == "", mode = "both"; end
    if intMode == "", intMode = "normalized"; end

    stem = sprintf('%s_%s_%s_%s_%s', ...
        strrep(char(stageName), ' ', '_'), ...
        char(cond), char(tp), char(mode), char(intMode));

    stem = regexprep(stem, '[^\w\-]', '_');
end