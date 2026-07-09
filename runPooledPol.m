function pooledResults = runPooledPol(baseDirs, opts)
% runPooledPol
%
% Multi-experiment analysis pipeline for mitochondrial spatial intensity
% profiles, front/back polarisation, and regional GFP-mCherry segregation
% across meiotic stages, treatment conditions, and timepoints.
%
% -------------------------------------------------------------------------
% OVERVIEW
% -------------------------------------------------------------------------
% This function is the top-level coordinator for the pooled mitochondrial
% spatial-profile analysis pipeline.
%
% For each selected:
%
%       Timepoint × Condition × Meiotic Stage
%
% the function calls mitoIntensityFS_pooled to process and concatenate
% single-cell GFP and/or mCherry intensity profiles across experiments.
%
% The pooled profiles are then used to:
%
%       - calculate front/back polarisation ratios
%       - generate stage-wise polarisation summaries
%       - generate polarisation timecourse plots
%       - generate pooled GFP profile overlays
%       - compare GFP profiles across timepoints and conditions
%       - calculate an experiment-level Regional GFP-mCherry
%         Segregation Score
%
% The function supports GFP-only, mCherry-only, or paired GFP/mCherry
% analyses and can use either per-cell normalised or raw intensity profiles.
%
% -------------------------------------------------------------------------
% PIPELINE HIERARCHY
% -------------------------------------------------------------------------
%
% runPooledPol
%       |
%       +-- mitoIntensityFS_pooled
%               |
%               +-- mitoIntensityFS
%
% Responsibilities:
%
%   runPooledPol
%       Coordinates timepoints, conditions, stages, pooled polarisation
%       analysis, Regional Segregation Score calculation, and global plots.
%
%   mitoIntensityFS_pooled
%       Repeats the single-experiment analysis across multiple experiment
%       directories, concatenates profiles, and preserves cell provenance.
%
%   mitoIntensityFS
%       Discovers CSV files, filters by treatment/timepoint, pairs channels,
%       normalises spatial position, optionally normalises intensity, and
%       bins each single-cell profile.
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% baseDirs : cell array, string array, or compatible list of paths
%     Experiment root directories.
%
%     Each experiment directory is expected to contain stage-specific
%     channel folders using the structure:
%
%         <baseDir>\G<stageName>*.csv
%         <baseDir>\R<stageName>*.csv
%
%     For the standard meiotic stages:
%
%         <baseDir>\G\Horsetail*.csv
%         <baseDir>\G\M1*.csv
%         <baseDir>\G\M2*.csv
%         <baseDir>\G\M2l*.csv
%         <baseDir>\G\Spores*.csv
%
%         <baseDir>\R\Horsetail*.csv
%         <baseDir>\R\M1*.csv
%         <baseDir>\R\M2*.csv
%         <baseDir>\R\M2l*.csv
%         <baseDir>\R\Spores*.csv
%
%     Only the folders required by opts.channelMode are used.
%
% opts : structure containing analysis and plotting options.
%
% -------------------------------------------------------------------------
% REQUIRED ON MATLAB PATH
% -------------------------------------------------------------------------
%
%       mitoIntensityFS_pooled.m
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
% 1) MULTI-DIMENSIONAL ANALYSIS COORDINATION
%    - Iterates through selected:
%
%          timepoints
%          treatment conditions
%          meiotic stages
%
%    - For every Timepoint × Condition × Stage combination, pooled profiles
%      are generated across all supplied experiment directories.
%
%    - Default timepoint keys:
%
%          ["5min","30min","60min"]
%
%    - Default conditions:
%
%          ["Veh","H2O2"]
%
%    - Default stages:
%
%          ["Horsetail","M1","M2","M2l","Spores"]
%
%    - Custom stage subsets are supported through opts.stages.
%
%      For example:
%
%          opts.stages = "Horsetail";
%
% 2) CHANNEL-AWARE ANALYSIS
%    - opts.channelMode controls which profiles are analysed:
%
%          "both"   paired GFP and mCherry
%          "gfp"    GFP only
%          "rfp"    mCherry/RFP only
%
%    - Channel-dependent outputs are skipped automatically when the
%      required channel is not available.
%
%    - The Regional Segregation Score requires both GFP and mCherry.
%
%    - GFP profile overlays are generated only when GFP is enabled.
%
% 3) NORMALISED OR RAW INTENSITY PROFILES
%    - opts.intensityMode controls y-intensity processing:
%
%          "normalized"
%              Per-cell min-max intensity scaling is performed by
%              mitoIntensityFS before profiles are pooled.
%
%          "raw"
%              Original fluorescence intensity values are retained.
%
%    - Spatial position is normalised to cell length in both modes.
%
%    - The aliases "normalised" and "norm" are converted internally to
%      "normalized".
%
% 4) POOLED PROFILE GENERATION
%    - Calls mitoIntensityFS_pooled separately for every selected
%      Timepoint × Condition × Stage combination.
%
%    - Single-cell binned profiles from all experiments are concatenated
%      into pooled GFP and/or mCherry profile matrices.
%
%    - Per-cell metadata from the pooling layer are retained for downstream
%      experiment-level analyses.
%
%    - Missing folders or failed combinations are reported and skipped
%      without terminating the complete pooled analysis.
%
% 5) FRONT/BACK POLARISATION ANALYSIS
%    - Calculates front and back mean intensities from configurable profile
%      bin ranges.
%
%    - Default regions:
%
%          Front = bins 1:50
%          Back  = bins 51:100
%
%    - For each valid cell profile:
%
%          Polarisation = Front mean / Back mean
%
%    - Channel-specific outputs are generated as:
%
%          Front_GFP
%          Back_GFP
%          Pol_GFP
%
%          Front_RFP
%          Back_RFP
%          Pol_RFP
%
%      depending on opts.channelMode.
%
%    - Non-finite values and profiles with Back mean = 0 are excluded from
%      the polarisation result.
%
% 6) POLARISATION METADATA AND TABLE EXPORT
%    - Polarisation values are combined with pooled profile metadata when
%      available.
%
%    - Additional analysis metadata are added:
%
%          StageFolder
%          StagePretty
%          Condition
%          ConditionPretty
%          Timepoint
%          ChannelModeUsed
%          IntensityModeUsed
%
%    - A polarisation CSV is written for every valid
%      Timepoint × Condition × Stage combination.
%
%    - All valid polarisation tables are concatenated into one master table.
%
% 7) POOLED PROFILE PLOTS
%    - opts.makeProfilePlots controls pooled GFP/mCherry profile plots for
%      individual Timepoint × Condition × Stage combinations.
%
%    - Profiles show:
%
%          pooled mean profile
%          SEM shading
%
%    - GFP is displayed in green.
%    - mCherry is displayed in magenta.
%    - The pooled profile figures use an A4 landscape-friendly layout.
%
%    - opts.profileYLim can enforce a common y-axis range.
%
% 8) STAGE-WISE POLARISATION SUMMARY PLOTS
%    - Enabled with:
%
%          opts.makeSummaryPlots = true
%
%    - Generates one summary plot for each timepoint and condition.
%
%    - Individual cell polarisation values are displayed as jittered points.
%    - Stage medians are displayed as horizontal lines.
%    - Cell counts are annotated near the medians.
%
%    - In paired mode, GFP and mCherry polarisation values are shown
%      side-by-side within each stage.
%
%    - opts.ylimSummary controls the y-axis limits for these figures.
%
%    - opts.ylim is retained as a backward-compatible general y-limit.
%
% 9) POOLED POLARISATION TIMECOURSE
%    - Enabled with:
%
%          opts.makeTimecoursePlot = true
%
%    - Generates timepoint panels containing Vehicle and H2O2 values across
%      the selected meiotic stages.
%
%    - Individual cell values and pooled medians are displayed.
%
%    - When GFP is enabled, the timecourse plot uses Pol_GFP.
%    - Pol_RFP is used when the analysis is RFP-only.
%
%    - Therefore, channelMode="both" currently generates a GFP
%      polarisation timecourse rather than separate GFP and RFP timecourses.
%
%    - opts.ylimTimecourse controls the y-axis range.
%
% 10) GFP STAGE OVERLAY PLOTS
%    - Enabled with:
%
%          opts.makeGFPStageOverlay = true
%
%    - For each Timepoint × Condition combination, pooled GFP profiles from
%      all selected stages are superimposed.
%
%    - Stage profiles use a light-to-dark green gradient.
%    - Mean profiles are displayed with a black underlay for visibility.
%    - SEM regions are shaded without outlines.
%
%    - These plots are skipped when GFP is not enabled.
%
% 11) GFP TIMEPOINT / CONDITION SUMMARY PLOTS
%    - Enabled with:
%
%          opts.makeGFPStageTimeCondSummary = true
%
%    - Generates one figure per meiotic stage.
%
%    - Each figure overlays pooled GFP profiles from all selected
%      timepoint and condition combinations.
%
%    - Vehicle profiles use a green gradient.
%    - H2O2 profiles use a grey gradient.
%    - Colour shade changes across timepoints.
%
%    - Mean profiles are displayed with a black underlay and coloured
%      foreground line.
%
%    - SEM shading is adjusted to retain visibility for light profiles.
%
% 12) REGIONAL GFP-mCHERRY SEGREGATION SCORE
%    - Enabled with:
%
%          opts.makeRegionalSegScore = true
%
%    - Requires:
%
%          channelMode="both"
%          no_points=100
%
%    - The score is calculated PER EXPERIMENT for every:
%
%          Timepoint × Condition × Stage
%
%    - Cells are assigned to experiments using pooled metadata.
%
%    - Experiment identity is preferentially extracted from the first
%      YYYYMMDD-style date found in ExperimentPath.
%
%    - If no date is found, the final experiment folder name is used.
%      "ExpN" is used as a final fallback.
%
%    - For each experiment:
%
%          1. Calculate the mean GFP profile.
%          2. Calculate the mean mCherry profile.
%          3. Partition both profiles into stage-dependent regions.
%
%    - Region definitions:
%
%          Horsetail, M1, M2:
%
%              Region 1 = bins 1:50
%              Region 2 = bins 51:100
%
%          M2l, Spores:
%
%              Q1 = bins 1:25
%              Q2 = bins 26:50
%              Q3 = bins 51:75
%              Q4 = bins 76:100
%
%    - Within each region:
%
%          GFP_Area = sum(mean GFP profile bins)
%          RFP_Area = sum(mean mCherry profile bins)
%
%          Signed difference = RFP_Area - GFP_Area
%
%          Absolute regional difference =
%              abs(RFP_Area - GFP_Area)
%
%    - Final score:
%
%          Regional Segregation Score =
%              sum(Absolute regional differences)
%
%    - Higher scores indicate stronger regional separation between GFP and
%      mCherry profiles.
%
%    - Lower scores indicate greater regional overlap or mixing.
%
%    - Each dot in the Regional Segregation Score plot represents one
%      experiment.
%
%    - The horizontal summary line in this plot represents the mean of the
%      experiment-level scores.
%
%    - opts.ylimRegionalSeg can enforce the y-axis range.
%
% 13) OUTPUT DIRECTORY AND FIGURE CONTROL
%    - opts.outputDir defines the analysis root directory.
%
%    - If no output directory is supplied, the default is:
%
%          <first baseDir>\PooledPolarisation
%
%    - Stage-specific outputs are organised as:
%
%          <outputDir><timepoint><condition><stage>
%
%    - Global summary outputs are saved to:
%
%          <outputDir>\SummaryPlots
%
%    - Major plot families can be independently enabled or disabled using:
%
%          makeProfilePlots
%          makeSummaryPlots
%          makeTimecoursePlot
%          makeRegionalSegScore
%          makeGFPStageOverlay
%          makeGFPStageTimeCondSummary
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% pooledResults : table
%     Combined polarisation table containing one row per retained valid
%     cell/profile result.
%
%     Depending on channelMode, the table contains:
%
%         Front_GFP
%         Back_GFP
%         Pol_GFP
%
%     and/or:
%
%         Front_RFP
%         Back_RFP
%         Pol_RFP
%
%     together with pooled cell provenance and analysis metadata when
%     available.
%
% The pipeline can additionally generate:
%
%   - PooledPol_<stage>*<condition>*<timepoint>*<mode>.csv
%       Per-combination front/back polarisation tables.
%
%   - PooledPolarisation_AllResults*<mode>.csv
%       Combined polarisation table.
%
%   - RegionalSegregationScore_PerExperiment.csv
%       One Regional Segregation Score per experiment, stage, condition,
%       and timepoint.
%
%   - RegionalSegregationScore_PerExperiment_Segments.csv
%       Region-level GFP and mCherry profile sums and differences.
%
%   - RegionalSegregationScore_PerExperiment_DotPlot_fromCSV.png/.fig
%
%   - Pooled stage polarisation summary plots.
%
%   - Pooled polarisation timecourse plots.
%
%   - Pooled GFP stage-overlay plots.
%
%   - Pooled GFP timepoint/condition summary plots.
%
%   - Pooled GFP/mCherry profile plots generated through
%       mitoIntensityFS_pooled.
%
%   - Per-experiment pairing reports and optional profile plots generated
%       through the daughter functions.
%
% -------------------------------------------------------------------------
% DESIGN PRINCIPLES
% -------------------------------------------------------------------------
% - Separates single-experiment profile processing from multi-experiment
%   pooling and high-level biological summaries.
%
% - Preserves cell and experiment provenance wherever pooled metadata are
%   available.
%
% - Supports paired and single-channel analyses through one shared pipeline.
%
% - Keeps spatial normalisation independent from intensity normalisation:
%   cell length is always normalised, while fluorescence intensity can remain
%   raw.
%
% - Uses cell-level profiles for front/back polarisation visualisation and
%   descriptive pooled summaries.
%
% - Calculates the Regional Segregation Score at the experiment level to
%   preserve experiment as the unit represented by each score.
%
% - The current function performs descriptive analysis and visualisation;
%   it does not perform inferential statistical hypothesis testing.
%
% -------------------------------------------------------------------------
% EXAMPLE USAGE
% -------------------------------------------------------------------------
%
% baseDirs = {
%     'V:\CSML\Raful\Experiment1'
%     'V:\CSML\Raful\Experiment2'
%     'V:\CSML\Raful\Experiment3'
%     'V:\CSML\Raful\Experiment4'
% };
%
% opts = struct();
%
% opts.outputDir = 'V:\CSML\Raful\PooledPolarisation';
%
% opts.stages = ["Horsetail","M1","M2","M2l","Spores"];
% opts.timeKeys = ["5min","60min"];
% opts.conditions = ["Veh","H2O2"];
% opts.condPretty = ["Vehicle","H2O2"];
%
% opts.channelMode = "both";
% opts.intensityMode = "normalized";
%
% opts.no_points = 100;
% opts.binsFront = 1:50;
% opts.binsBack = 51:100;
%
% opts.makeProfilePlots = true;
% opts.perExperimentPlots = false;
% opts.makeSummaryPlots = true;
% opts.makeTimecoursePlot = true;
% opts.makeRegionalSegScore = true;
% opts.makeGFPStageOverlay = true;
% opts.makeGFPStageTimeCondSummary = true;
%
% opts.profileYLim = [];
% opts.ylimSummary = [0 10];
% opts.ylimTimecourse = [0 3];
% opts.ylimRegionalSeg = [];
%
% pooledResults = runPooledPol(baseDirs, opts);
%
% -------------------------------------------------------------------------
% Author: Raful Navarro Espindola
% Project: Mitochondrial distribution and segregation during fission yeast
%          meiosis
% -------------------------------------------------------------------------

    if nargin < 2 || isempty(opts), opts = struct(); end

    % -------------------- defaults --------------------
    if ~isfield(opts,'outputDir') || isempty(opts.outputDir)
        opts.outputDir = fullfile(char(string(baseDirs(1))), 'PooledPolarisation');
    end
    if ~isfield(opts,'binsFront'),          opts.binsFront = 1:50; end
    if ~isfield(opts,'binsBack'),           opts.binsBack  = 51:100; end
    if ~isfield(opts,'matchByName'),        opts.matchByName = true; end
    if ~isfield(opts,'timeKeys') || isempty(opts.timeKeys)
        opts.timeKeys = ["5min","30min","60min"];
    end
    if ~isfield(opts,'conditions') || isempty(opts.conditions)
        opts.conditions = ["Veh","H2O2"];
    end
    if ~isfield(opts,'condPretty') || isempty(opts.condPretty)
        opts.condPretty = ["Vehicle","H2O2"];
    end
    if ~isfield(opts,'channelMode') || isempty(opts.channelMode)
        opts.channelMode = "both";
    end
    if ~isfield(opts,'makeRegionalSegScore'), opts.makeRegionalSegScore = true; end
    if ~isfield(opts,'ylimRegionalSeg'), opts.ylimRegionalSeg = []; end
    if ~isfield(opts,'no_points'),          opts.no_points = 100; end
    if ~isfield(opts,'makeProfilePlots'),   opts.makeProfilePlots = true; end
    if ~isfield(opts,'perExperimentPlots'), opts.perExperimentPlots = false; end
    if ~isfield(opts,'saveMeta'),           opts.saveMeta = true; end
    if ~isfield(opts,'verbose'),            opts.verbose = true; end
    if ~isfield(opts,'makeSummaryPlots'),   opts.makeSummaryPlots = true; end
    if ~isfield(opts,'makeTimecoursePlot'), opts.makeTimecoursePlot = true; end
    if ~isfield(opts,'intensityMode') || isempty(opts.intensityMode)
        opts.intensityMode = "normalized";
    end
    if ~isfield(opts,'profileYLim'), opts.profileYLim = []; end

    % Backward-compatible general y-limit option
    if ~isfield(opts,'ylim'), opts.ylim = []; end

    % Independent y-limits for different plot types
    if ~isfield(opts,'ylimSummary'), opts.ylimSummary = opts.ylim; end
    if ~isfield(opts,'ylimTimecourse'), opts.ylimTimecourse = opts.ylim; end

    if ~isfield(opts,'makeGFPStageOverlay'), opts.makeGFPStageOverlay = true; end
    if ~isfield(opts,'makeGFPStageTimeCondSummary'), opts.makeGFPStageTimeCondSummary = true; end

    baseDirs = string(baseDirs(:));
    opts.timeKeys = string(opts.timeKeys);
    opts.conditions = string(opts.conditions);
    opts.condPretty = string(opts.condPretty);
    opts.channelMode = lower(string(opts.channelMode));
    opts.intensityMode = lower(string(opts.intensityMode));
    if opts.intensityMode == "normalised" || opts.intensityMode == "norm"
        opts.intensityMode = "normalized";
    end

    if ~ismember(opts.channelMode, ["both","gfp","rfp"])
        error("opts.channelMode must be 'both', 'gfp', or 'rfp'.");
    end

    if ~ismember(opts.intensityMode, ["normalized","raw"])
        error("opts.intensityMode must be 'normalized' or 'raw'.");
    end

    [~, runGFP, runRFP] = parseChannelMode_local(opts.channelMode);

    % -------------------- stage selection --------------------
    % Default: run all standard meiotic stages.
    % To run only Horsetail, use:
    %   opts.stages = "Horsetail";
    if ~isfield(opts,'stages') || isempty(opts.stages)
        opts.stages = ["Horsetail","M1","M2","M2l","Spores"];
    end

    stages = string(opts.stages);

    pretty = strings(size(stages));
    shortLb = strings(size(stages));

    for ii = 1:numel(stages)
        switch stages(ii)
            case "Horsetail"
                pretty(ii) = "Horsetail";
                shortLb(ii) = "Ht";
            case "M1"
                pretty(ii) = "Meiosis I";
                shortLb(ii) = "M1";
            case "M2"
                pretty(ii) = "Meiosis II";
                shortLb(ii) = "M2";
            case "M2l"
                pretty(ii) = "Late Meiosis II";
                shortLb(ii) = "LM2";
            case "Spores"
                pretty(ii) = "Spores";
                shortLb(ii) = "S";
            otherwise
                pretty(ii) = stages(ii);
                shortLb(ii) = stages(ii);
        end
    end

    outRoot = char(opts.outputDir);
    if ~exist(outRoot,'dir'), mkdir(outRoot); end

    nT = numel(opts.timeKeys);
    nC = numel(opts.conditions);
    nS = numel(stages);

    allTables = cell(nT, nC, nS);
    gPolCell = cell(nT, nC, nS);
    rPolCell = cell(nT, nC, nS);
    gProfCell = cell(nT, nC, nS);
    rProfCell = cell(nT, nC, nS);
    profMetaCell = cell(nT, nC, nS);

    fprintf('\n============================================\n');
    fprintf('   POOLED POLARISATION ANALYSIS\n');
    fprintf('============================================\n');
    fprintf('Experiments: %d\n', numel(baseDirs));
    fprintf('Output: %s\n', outRoot);
    fprintf('Channel mode: %s\n', char(opts.channelMode));
    fprintf('Intensity mode: %s\n', char(opts.intensityMode));

    % -------------------- main loop --------------------
    for t = 1:nT
        tpKey = opts.timeKeys(t);

        fprintf('\n############################################\n');
        fprintf('   TIMEPOINT: %s\n', tpKey);
        fprintf('############################################\n');

        for c = 1:nC
            cond = opts.conditions(c);
            condName = opts.condPretty(min(c, numel(opts.condPretty)));

            fprintf('\n---- CONDITION: %s ----\n', condName);

            for s = 1:nS
                stageName = stages(s);
                prettyStage = pretty(s);

                fprintf('\nStage: %s\n', prettyStage);

                mitoOpts = struct();
                mitoOpts.no_points = opts.no_points;
                mitoOpts.makePlot = opts.makeProfilePlots;
                mitoOpts.matchByName = opts.matchByName;
                mitoOpts.condition = cond;
                mitoOpts.timepoint = tpKey;
                mitoOpts.channelMode = char(opts.channelMode);
                mitoOpts.saveMeta = opts.saveMeta;
                mitoOpts.verbose = opts.verbose;
                mitoOpts.perExperimentPlots = opts.perExperimentPlots;
                mitoOpts.intensityMode = opts.intensityMode;
                mitoOpts.profileYLim = opts.profileYLim;

                stageOutDir = fullfile(outRoot, char(tpKey), char(cond), char(stageName));
                if ~exist(stageOutDir,'dir'), mkdir(stageOutDir); end

                try
                    [distG, distR, pooledMeta] = mitoIntensityFS_pooled( ...
                        baseDirs, stageOutDir, char(stageName), mitoOpts);

                    if runGFP && ~isempty(distG)
                        gProfCell{t,c,s} = distG;
                    end

                    if runRFP && ~isempty(distR)
                        rProfCell{t,c,s} = distR;
                    end

                    if ~isempty(pooledMeta)
                        profMetaCell{t,c,s} = pooledMeta;
                    end

                    Tpol = compareFrontBack_tableOnly_local(distG, distR, opts);

                    if isempty(Tpol)
                        fprintf('  No valid polarisation values for %s | %s | %s\n', ...
                            stageName, cond, tpKey);
                        continue;
                    end

                    Tout = addMetadataToPolTable(Tpol, pooledMeta, stageName, prettyStage, cond, condName, tpKey, opts);

                    allTables{t,c,s} = Tout;

                    if runGFP && ismember('Pol_GFP', Tout.Properties.VariableNames)
                        gPolCell{t,c,s} = Tout.Pol_GFP;
                    end
                    if runRFP && ismember('Pol_RFP', Tout.Properties.VariableNames)
                        rPolCell{t,c,s} = Tout.Pol_RFP;
                    end

                    outCsv = fullfile(stageOutDir, sprintf('PooledPol_%s_%s_%s_%s.csv', ...
                        cleanTag(stageName), cleanTag(cond), cleanTag(tpKey), ...
                        cleanTag(opts.intensityMode)));
                    writetable(Tout, outCsv);
                    fprintf('  Saved: %s\n', outCsv);

                catch ME
                    warning('Skipping %s | %s | %s because: %s', ...
                        stageName, cond, tpKey, ME.message);
                    continue;
                end
            end

            if opts.makeSummaryPlots
                makePooledStageSummaryPlot(outRoot, tpKey, cond, condName, shortLb, ...
                    squeeze(gPolCell(t,c,:)), squeeze(rPolCell(t,c,:)), opts);
            end

            if opts.makeGFPStageOverlay
                makePooledGFPStageOverlayPlot(outRoot, tpKey, cond, condName, ...
                    stages, pretty, reshape(gProfCell(t,c,:), 1, []), opts);
            end
        end
    end

    % -------------------- concatenate and save all results --------------------
    pooledResults = table();
    for t = 1:nT
        for c = 1:nC
            for s = 1:nS
                if ~isempty(allTables{t,c,s})
                    pooledResults = [pooledResults; allTables{t,c,s}]; %#ok<AGROW>
                end
            end
        end
    end

    allCsv = fullfile(outRoot, sprintf('PooledPolarisation_AllResults_%s.csv', cleanTag(opts.intensityMode)));
    if ~isempty(pooledResults)
        writetable(pooledResults, allCsv);
        fprintf('\nSaved combined pooled polarisation table:\n%s\n', allCsv);
    else
        warning('No pooled polarisation results were generated.');
    end

    if opts.makeRegionalSegScore
        makeRegionalSegScorePlot(outRoot, stages, shortLb, ...
            gProfCell, rProfCell, profMetaCell, opts);
    end

    % -------------------- global plots --------------------
    if opts.makeTimecoursePlot && ~isempty(pooledResults)
        makePooledTimecoursePlot(outRoot, pooledResults, shortLb, opts);
    end
    if opts.makeGFPStageTimeCondSummary
        makeGFPStageTimeCondSummaryPlots(outRoot, stages, pretty, gProfCell, opts);
    end
end

% ========================================================================
function makeRegionalSegScorePlot(outRoot, stages, shortLb, ...
    gProfCell, rProfCell, profMetaCell, opts)
% makeRegionalSegScorePlot
%
% PER-EXPERIMENT REGIONAL GFP-RFP SEGREGATION SCORE
% -------------------------------------------------
% This output quantifies regional segregation between the GFP and RFP
% profiles, but it does so PER EXPERIMENT rather than on the fully pooled
% dataset.
%
% For each timepoint x condition x stage x experiment:
%
%   1. Extract the cells belonging to one experiment using pooledMeta.
%
%   2. Compute the experiment-level mean GFP profile.
%
%   3. Compute the experiment-level mean RFP profile.
%
%   4. Partition the profiles into biologically meaningful regions:
%
%        Horsetail, M1, M2:
%            Region 1 = bins 1:50
%            Region 2 = bins 51:100
%
%        M2l, Spores:
%            Region 1 = bins 1:25
%            Region 2 = bins 26:50
%            Region 3 = bins 51:75
%            Region 4 = bins 76:100
%
%   5. For each region, calculate:
%
%            GFP_area = sum(meanGFP(region))
%            RFP_area = sum(meanRFP(region))
%
%   6. Calculate the absolute regional difference:
%
%            abs(RFP_area - GFP_area)
%
%   7. Sum the regional differences:
%
%            Regional segregation score = sum(abs regional differences)
%
% Interpretation:
%   Each dot in the output plot represents one experiment.
%
%   Higher values indicate stronger regional separation between the GFP
%   and RFP profiles.
%
%   Lower values indicate greater regional overlap or mixing.
%
% This figure should be interpreted as an experiment-level summary, not as
% a cell-level pooled score.

    [~, runGFP, runRFP] = parseChannelMode_local(opts.channelMode);

    if ~(runGFP && runRFP)
        warning("Skipping Regional segregation score because both GFP and RFP are required.");
        return;
    end

    if opts.no_points ~= 100
        error("Regional segregation score currently expects opts.no_points = 100.");
    end

    timepoints = string(opts.timeKeys);
    conditions = string(opts.conditions);
    condPretty = string(opts.condPretty);

    nT = numel(timepoints);
    nC = numel(conditions);
    nS = numel(stages);

    scoreTable = table();
    segmentOut = table();

    for t = 1:nT
        for c = 1:nC
            for s = 1:nS

                G = gProfCell{t,c,s};
                R = rProfCell{t,c,s};
                M = profMetaCell{t,c,s};

                if isempty(G) || isempty(R) || isempty(M)
                    continue;
                end

                % Keep only the matched number of cells across GFP, RFP, and metadata
                nCellsTotal = min([size(G,2), size(R,2), height(M)]);
                G = G(:,1:nCellsTotal);
                R = R(:,1:nCellsTotal);
                M = M(1:nCellsTotal,:);

                % Get experiment identity from the YYYYMMDD date in the full experiment path.
                % This avoids collapsing different experiments that share the same final
                % folder name but belong to different acquisition dates.
                if ismember('ExperimentPath', M.Properties.VariableNames)

                    expPaths = string(M.ExperimentPath);
                    expNames = strings(size(expPaths));

                    for ii = 1:numel(expPaths)
                        expNames(ii) = extractExperimentDateID(expPaths(ii), ii);
                    end

                elseif ismember('ExperimentName', M.Properties.VariableNames)

                    % Fallback only. This may collapse experiments if folder names repeat.
                    warning(['ExperimentPath was not found in pooledMeta. ', ...
                        'Falling back to ExperimentName, which may collapse experiments.']);

                    expNames = string(M.ExperimentName);

                else
                    warning('No ExperimentPath or ExperimentName found in pooledMeta for %s | %s | %s.', ...
                        stages(s), conditions(c), timepoints(t));
                    continue;
                end

                uniqueExps = unique(expNames, 'stable');

                for e = 1:numel(uniqueExps)
                    thisExp = uniqueExps(e);
                    idx = expNames == thisExp;

                    if ~any(idx)
                        continue;
                    end

                    meanG = mean(G(:,idx), 2, 'omitnan');
                    meanR = mean(R(:,idx), 2, 'omitnan');

                    [score, segmentTable] = computeRegionalSegScore(meanG, meanR, stages(s));

                    nCellsExp = sum(idx);

                    % One row per experiment-level score
                    tmpScore = table( ...
                        timepoints(t), ...
                        conditions(c), ...
                        condPretty(c), ...
                        string(stages(s)), ...
                        string(shortLb(s)), ...
                        thisExp, ...
                        nCellsExp, ...
                        score, ...
                        'VariableNames', { ...
                            'Timepoint', ...
                            'Condition', ...
                            'ConditionPretty', ...
                            'Stage', ...
                            'StageShort', ...
                            'ExperimentName', ...
                            'NCells', ...
                            'RegionalSegregationScore'});

                    scoreTable = [scoreTable; tmpScore]; %#ok<AGROW>

                    % Segment-level output for QC / optional future stacked plots
                    nSeg = height(segmentTable);

                    segmentTable.Timepoint = repmat(timepoints(t), nSeg, 1);
                    segmentTable.Condition = repmat(conditions(c), nSeg, 1);
                    segmentTable.ConditionPretty = repmat(condPretty(c), nSeg, 1);
                    segmentTable.Stage = repmat(string(stages(s)), nSeg, 1);
                    segmentTable.StageShort = repmat(string(shortLb(s)), nSeg, 1);
                    segmentTable.ExperimentName = repmat(thisExp, nSeg, 1);
                    segmentTable.NCells = repmat(nCellsExp, nSeg, 1);
                    segmentTable.RegionalSegregationScore = repmat(score, nSeg, 1);

                    segmentTable = segmentTable(:, { ...
                        'Timepoint', ...
                        'Condition', ...
                        'ConditionPretty', ...
                        'Stage', ...
                        'StageShort', ...
                        'ExperimentName', ...
                        'NCells', ...
                        'RegionLabel', ...
                        'BinStart', ...
                        'BinEnd', ...
                        'GFP_Area', ...
                        'RFP_Area', ...
                        'Signed_RFP_minus_GFP', ...
                        'Absolute_Regional_Difference', ...
                        'RegionalSegregationScore'});

                    segmentOut = [segmentOut; segmentTable]; %#ok<AGROW>
                end
            end
        end
    end

    outDir = fullfile(outRoot, 'SummaryPlots');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    if isempty(scoreTable)
        warning('No per-experiment Regional segregation score values were generated.');
        return;
    end

    % Save experiment-level scores
    outCsv = fullfile(outDir, 'RegionalSegregationScore_PerExperiment.csv');
    writetable(scoreTable, outCsv);
    fprintf('Saved per-experiment Regional segregation score table:\n%s\n', outCsv);

    % Save segment-level values
    if ~isempty(segmentOut)
        outSegCsv = fullfile(outDir, 'RegionalSegregationScore_PerExperiment_Segments.csv');
        writetable(segmentOut, outSegCsv);
        fprintf('Saved per-experiment Regional segregation segment table:\n%s\n', outSegCsv);
    end

    % -------------------- shared y-axis --------------------
    if isfield(opts,'ylimRegionalSeg') && ~isempty(opts.ylimRegionalSeg)
        yL = opts.ylimRegionalSeg;
    else
        finiteVals = scoreTable.RegionalSegregationScore;
        finiteVals = finiteVals(isfinite(finiteVals));

        if isempty(finiteVals)
            yL = [0 1];
        else
            yL = [0 1.10 * max(finiteVals)];
        end
    end

    % -------------------- dot plot: one continuous stage/timepoint plot --------------------
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

    % Dot colours for Regional SegScore plot
    vehCol = [0.60 0.30 0.75];   % Vehicle: purple
    h2Col  = [0.60 0.60 0.60];   % H2O2: grey

    condCols = zeros(nC,3);
    for c = 1:nC
        if conditions(c) == "Veh"
            condCols(c,:) = vehCol;
        elseif conditions(c) == "H2O2"
            condCols(c,:) = h2Col;
        else
            condCols(c,:) = [0.35 0.35 0.35];
        end
    end

    % Make timepoint labels prettier: "5min" -> "5 min"
    tpLabels = string(regexprep(cellstr(timepoints), '(\d+)\s*min', '$1 min'));

    % Layout:
    % For each stage:
    %   Vehicle 5 min | H2O2 5 min | Vehicle 60 min | H2O2 60 min
    timeGap = 1.0;     % increase this slightly if you want more separation between 5 and 60 min
    stageGap = 1.5;    % separation between stages

    groupWidth = nT*nC + (nT-1)*timeGap;
    stageStarts = (0:nS-1) * (groupWidth + stageGap);

    jitter = 0.07;
    meanHalfWidth = 0.18;

    tickPos = [];
    tickLabels = strings(0);

    for s = 1:nS

        for t = 1:nT
            for c = 1:nC

                idx = scoreTable.Timepoint == timepoints(t) & ...
                    scoreTable.Stage == string(stages(s)) & ...
                    scoreTable.Condition == conditions(c);

                vals = scoreTable.RegionalSegregationScore(idx);
                vals = vals(isfinite(vals));

                if isempty(vals)
                    continue;
                end

                xCenter = stageStarts(s) + (t-1)*(nC + timeGap) + c;
                xVals = xCenter + jitter * randn(size(vals));

                scatter(xVals, vals, 100, ...
                    'filled', ...
                    'MarkerFaceColor', condCols(c,:), ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 0.5, ...
                    'MarkerFaceAlpha', 0.75);

                meanVal = mean(vals, 'omitnan');

                line([xCenter-meanHalfWidth, xCenter+meanHalfWidth], ...
                    [meanVal, meanVal], ...
                    'Color', 'k', ...
                    'LineWidth', 2.2);
            end

            % One tick label per timepoint, centred over Vehicle/H2O2 pair
            tickPos(end+1) = stageStarts(s) + (t-1)*(nC + timeGap) + (nC+1)/2;
            tickLabels(end+1) = tpLabels(t); %#ok<AGROW>
        end
    end

    % Axis limits and labels
    xlim([0.5, stageStarts(end) + groupWidth + 0.5]);
    ylim(yL);

    xticks(tickPos);
    xticklabels(tickLabels);
    xtickangle(0);

    % Add stage labels underneath the timepoint labels
    yRange = diff(yL);
    stageLabelY = yL(1) - 0.09 * yRange;

    for s = 1:nS
        stageMid = stageStarts(s) + (1 + groupWidth)/2;

        text(stageMid, stageLabelY, shortLb(s), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontWeight', 'bold', ...
            'FontSize', 24, ...
            'Clipping', 'off');
    end
    % -------------------- axis styling --------------------
    % No x-axis label
    xlabel('');

    % Y-axis label handle
    yl = ylabel('Regional GFP-mCherry segregation score');

    grid off;
    box off;

    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;

    % Tick labels: numeric y-axis labels + "5 min"/"60 min" x-axis labels
    ax.FontSize = 18;
    ax.FontWeight = 'normal';

    % Now force only the y-axis label to be large and bold
    yl.FontSize = 24;
    yl.FontWeight = 'bold';

    % Make sure labels fit inside the A4 page
    ax.Position = [0.12 0.22 0.74 0.68];

 % -------------------- legend --------------------
    %hLeg = gobjects(nC,1);
    %for c = 1:nC
    %    hLeg(c) = plot(nan, nan, 'o', ...
    %        'MarkerFaceColor', condCols(c,:), ...
    %        'MarkerEdgeColor', 'k', ...
    %        'LineStyle', 'none', ...
    %        'MarkerSize', 7);
    %end

    %legend(hLeg, condPretty(1:nC), ...
    %    'Location', 'northeastoutside');

    % -------------------- save --------------------
    outPng = fullfile(outDir, 'RegionalSegregationScore_PerExperiment_DotPlot_fromCSV.png');
    outFig = fullfile(outDir, 'RegionalSegregationScore_PerExperiment_DotPlot_fromCSV.fig');

    exportgraphics(f, outPng, 'Resolution', 300);
    savefig(f, outFig);
    close(f);

    fprintf('Saved RegSeg plot from CSV:\n%s\n', outPng);
end

% ========================================================================
function expID = extractExperimentDateID(pathIn, fallbackIdx)
% extractExperimentDateID
%
% Extracts the first YYYYMMDD-style date from a full experiment path and
% uses it as the experiment identifier.
%
% Example:
%   V:\...\20260321 H2O2Ac WTmcp5GFP Fzo1 WT\VA134xVA192...
%
% returns:
%   "20260321"
%
% If no date is found, it falls back to the final folder name. If that also
% fails, it returns "ExpN".

    pathIn = string(pathIn);

    % Look for dates such as 20260321, 20260402, etc.
    tok = regexp(char(pathIn), '20\d{6}', 'match', 'once');

    if ~isempty(tok)
        expID = string(tok);
        return;
    end

    % Fallback: use final folder name
    [~, folderName] = fileparts(char(pathIn));

    if ~isempty(folderName)
        expID = string(folderName);
    else
        expID = "Exp" + string(fallbackIdx);
    end
end

% ========================================================================
function [score, segmentTable] = computeRegionalSegScore(meanG, meanR, stageName)
% computeRegionalSegScore
%
% Calculates the Regional GFP-RFP segregation score for one pair of mean
% profiles.
%
% In the current pipeline, meanG and meanR correspond to the mean GFP and
% RFP profiles from ONE experiment for one stage x condition x timepoint.
%
% The calculation is region-first:
%
%   1. Partition the 100-bin profile into biologically meaningful regions.
%
%   2. Calculate the area under the GFP and RFP mean profiles within each
%      region.
%
%   3. Subtract the regional areas:
%
%            signed difference = RFP_area - GFP_area
%
%   4. Take the absolute value of each regional difference.
%
%   5. Sum those absolute regional differences to obtain the final score.
%
% Interpretation:
%   Higher scores indicate stronger regional GFP-RFP segregation.
%   Lower scores indicate greater regional GFP-RFP overlap or mixing.

    stageName = string(stageName);

    switch stageName
        case {"Horsetail","M1","M2"}
            regionLabels = ["Front half"; "Back half"];
            binStarts = [1; 51];
            binEnds   = [50; 100];

        case {"M2l","Spores"}
            regionLabels = ["Q1"; "Q2"; "Q3"; "Q4"];
            binStarts = [1; 26; 51; 76];
            binEnds   = [25; 50; 75; 100];

        otherwise
            error('Unknown stage name for regional score: %s', stageName);
    end

    nRegions = numel(regionLabels);

    gfpArea = nan(nRegions,1);
    rfpArea = nan(nRegions,1);
    signedDiff = nan(nRegions,1);
    absDiff = nan(nRegions,1);

    for k = 1:nRegions
        bins = binStarts(k):binEnds(k);

        gfpArea(k) = sum(meanG(bins), 'omitnan');
        rfpArea(k) = sum(meanR(bins), 'omitnan');

        signedDiff(k) = rfpArea(k) - gfpArea(k);
        absDiff(k) = abs(signedDiff(k));
    end

    score = sum(absDiff, 'omitnan');

    segmentTable = table( ...
        regionLabels, ...
        binStarts, ...
        binEnds, ...
        gfpArea, ...
        rfpArea, ...
        signedDiff, ...
        absDiff, ...
        'VariableNames', { ...
            'RegionLabel', ...
            'BinStart', ...
            'BinEnd', ...
            'GFP_Area', ...
            'RFP_Area', ...
            'Signed_RFP_minus_GFP', ...
            'Absolute_Regional_Difference'});
end

% ========================================================================
function makeGFPStageTimeCondSummaryPlots(outRoot, stages, pretty, gProfCell, opts)

    [~, runGFP, ~] = parseChannelMode_local(opts.channelMode);

    if ~runGFP
        warning("Skipping GFP stage/time-condition summary plots because channelMode does not include GFP.");
        return;
    end

    timepoints = string(opts.timeKeys);
    conditions = string(opts.conditions);
    condPretty = string(opts.condPretty);

    nT = numel(timepoints);
    nC = numel(conditions);
    nS = numel(stages);

    x = 1:opts.no_points;

    % ---- Vehicle palette: green gradient ----
    greenStart = [0.85 0.97 0.85];
    greenEnd   = [0.00 0.45 0.00];
    vehCols = [ ...
        linspace(greenStart(1), greenEnd(1), nT)', ...
        linspace(greenStart(2), greenEnd(2), nT)', ...
        linspace(greenStart(3), greenEnd(3), nT)' ...
    ];

    % ---- H2O2 palette: grey gradient ----
    greyStart = [0.80 0.80 0.80];
    greyEnd   = [0.35 0.35 0.35];
    h2Cols = [ ...
        linspace(greyStart(1), greyEnd(1), nT)', ...
        linspace(greyStart(2), greyEnd(2), nT)', ...
        linspace(greyStart(3), greyEnd(3), nT)' ...
    ];

    outDir = fullfile(outRoot, 'SummaryPlots');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    for s = 1:nS
        f = figure('Color','w');
        hold on;

        h = gobjects(0);
        labels = {};
        profilesForY = {};

        for t = 1:nT
            for c = 1:nC
                G = gProfCell{t,c,s};

                if isempty(G)
                    continue;
                end

                profilesForY{end+1} = G; %#ok<AGROW>

                meanG = mean(G, 2, 'omitnan');
                nG = sum(~all(isnan(G),1));
                semG = std(G, 0, 2, 'omitnan') ./ sqrt(max(nG,1));

                xv = [x, fliplr(x)];
                yv = [transpose(meanG - semG), fliplr(transpose(meanG + semG))];

                % ---- choose colour by condition and timepoint ----
                if conditions(c) == "Veh"
                    curveCol = vehCols(t,:);
                else
                    curveCol = h2Cols(t,:);
                end

                % ---- make lighter profiles slightly more visible ----
                faceAlpha = 0.08;
                semCol = curveCol;

                if conditions(c) == "Veh" && t == 1
                    semCol = curveCol * 0.75;
                    faceAlpha = 0.20;
                elseif conditions(c) == "Veh" && t == 2
                    semCol = curveCol * 0.90;
                    faceAlpha = 0.12;
                elseif conditions(c) ~= "Veh" && t == 1
                    semCol = curveCol * 0.90;
                    faceAlpha = 0.12;
                end

                % ---- SEM shading ----
                fill(xv, yv, semCol, ...
                    'FaceAlpha', faceAlpha, ...
                    'EdgeColor', 'none');

                % ---- black outline ----
                plot(x, meanG, ...
                    'Color', 'k', ...
                    'LineWidth', 3.0, ...
                    'LineStyle', '-');

                % ---- coloured mean line (always solid) ----
                h(end+1) = plot(x, meanG, ...
                    'Color', curveCol, ...
                    'LineWidth', 1.8, ...
                    'LineStyle', '-'); %#ok<AGROW>

                labels{end+1} = sprintf('%s %s (n=%d)', ...
                    timepoints(t), condPretty(c), nG); %#ok<AGROW>
            end
        end

        if isempty(h)
            close(f);
            warning('No GFP pooled profiles available for stage %s', stages(s));
            continue;
        end

        axis square;
        box on;
        grid on;
        xlim([1 opts.no_points]);

        applyProfileYLim_local(profilesForY, opts);

        xlabel('Normalized Cell Length (binned)');
        ylabel(profileYLabel_local(opts, "GFP"));

        [modeTag, modeTitle] = intensityModeTagTitle_local(opts);
        title(sprintf('Pooled GFP profiles across timepoints and conditions — %s — %s', ...
            pretty(s), modeTitle), ...
            'FontWeight', 'bold', 'Interpreter', 'none');

        legend(h, labels, 'Location', 'northeastoutside');

        stem = sprintf('PooledGFP_TimeCond_%s_%s', cleanTag(stages(s)), cleanTag(modeTag));
        saveas(f, fullfile(outDir, [stem '.png']));
        saveas(f, fullfile(outDir, [stem '.fig']));
        close(f);

        fprintf('Saved pooled GFP time/condition summary plot:\n%s\n', ...
            fullfile(outDir, [stem '.png']));
    end
end

% ========================================================================
function makePooledGFPStageOverlayPlot(outRoot, tpKey, cond, condName, ...
    stages, pretty, gProfByStage, opts)

    [~, runGFP, ~] = parseChannelMode_local(opts.channelMode);

    if ~runGFP
        warning("Skipping GFP stage overlay plot because channelMode does not include GFP.");
        return;
    end

    nStages = numel(stages);
    x = 1:opts.no_points;

    % ---- light-to-dark green gradient ----
    greenStart = [0.85 0.97 0.85];   % light green
    greenEnd   = [0.00 0.45 0.00];   % dark green

    cols = [ ...
        linspace(greenStart(1), greenEnd(1), nStages)', ...
        linspace(greenStart(2), greenEnd(2), nStages)', ...
        linspace(greenStart(3), greenEnd(3), nStages)' ...
        ];

    % ---- use solid lines for all stages ----
    lineStyles = repmat({'-'}, 1, nStages);

    % ---- SEM appearance ----
    % Keep the SEM clean by removing outlines.
    % Make the lighter green stages more visible against the white background.
    semCols = cols;
    semAlpha = 0.08 * ones(1, nStages);

    if nStages >= 1
        semCols(1,:) = cols(1,:) * 0.75;   % Horsetail
        semAlpha(1) = 0.22;
    end

    if nStages >= 2
        semCols(2,:) = cols(2,:) * 0.82;   % M1
        semAlpha(2) = 0.16;
    end

    f = figure('Color','w');
    hold on;

    h = gobjects(0);
    labels = {};
    profilesForY = {};

    for s = 1:nStages
        G = gProfByStage{s};

        if isempty(G)
            continue;
        end

        profilesForY{end+1} = G; %#ok<AGROW>

        meanG = mean(G, 2, 'omitnan');
        nG = sum(~all(isnan(G), 1));
        semG = std(G, 0, 2, 'omitnan') ./ sqrt(max(nG,1));

        xv = [x, fliplr(x)];
        yv = [transpose(meanG - semG), fliplr(transpose(meanG + semG))];

        % ---- SEM shading without outline ----
        fill(xv, yv, semCols(s,:), ...
            'FaceAlpha', semAlpha(s), ...
            'EdgeColor', 'none');

        % ---- black underlay for the mean profile ----
        plot(x, meanG, ...
            'Color', 'k', ...
            'LineWidth', 3.0, ...
            'LineStyle', lineStyles{s});

        % ---- coloured mean profile on top ----
        h(end+1) = plot(x, meanG, ...
            'Color', cols(s,:), ...
            'LineWidth', 1.8, ...
            'LineStyle', lineStyles{s}); %#ok<AGROW>

        labels{end+1} = sprintf('%s (n=%d)', pretty(s), nG); %#ok<AGROW>
    end

    if isempty(h)
        close(f);
        warning('No GFP pooled profiles available for %s | %s', cond, tpKey);
        return;
    end

    axis square;
    box on;
    grid on;
    xlim([1 opts.no_points]);

    applyProfileYLim_local(profilesForY, opts);

    xlabel('Normalized Cell Length (binned)');
    ylabel(profileYLabel_local(opts, "GFP"));

    [modeTag, modeTitle] = intensityModeTagTitle_local(opts);
    title(sprintf('Pooled GFP profiles by stage — %s — %s — %s', ...
        condName, tpKey, modeTitle), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    legend(h, labels, 'Location', 'northeastoutside');

    outDir = fullfile(outRoot, 'SummaryPlots');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    stem = sprintf('PooledGFPStageOverlay_%s_%s_%s', cleanTag(cond), cleanTag(tpKey), cleanTag(modeTag));
    saveas(f, fullfile(outDir, [stem '.png']));
    saveas(f, fullfile(outDir, [stem '.fig']));
    close(f);

    fprintf('Saved pooled GFP stage overlay plot:\n%s\n', ...
        fullfile(outDir, [stem '.png']));
end

% ========================================================================
function Tout = addMetadataToPolTable(Tpol, pooledMeta, stageName, prettyStage, cond, condName, tpKey, opts)

    nPol = height(Tpol);

    if ~isempty(pooledMeta) && height(pooledMeta) >= nPol
        metaUse = pooledMeta(1:nPol,:);
        Tout = [metaUse, Tpol];
    else
        Tout = Tpol;
        Tout.ExperimentPath = strings(nPol,1);
        Tout.ExperimentName = strings(nPol,1);
        Tout.CellIndexWithinExperiment = nan(nPol,1);
        Tout.GFP_File = strings(nPol,1);
        Tout.RFP_File = strings(nPol,1);
    end

    Tout.StageFolder = repmat(string(stageName), nPol, 1);
    Tout.StagePretty = repmat(string(prettyStage), nPol, 1);
    Tout.Condition = repmat(string(cond), nPol, 1);
    Tout.ConditionPretty = repmat(string(condName), nPol, 1);
    Tout.Timepoint = repmat(string(tpKey), nPol, 1);
    Tout.ChannelModeUsed = repmat(string(opts.channelMode), nPol, 1);
    Tout.IntensityModeUsed = repmat(string(opts.intensityMode), nPol, 1);
end

% ========================================================================
function T = compareFrontBack_tableOnly_local(distProfileGFP, distProfileRFP, opts)
% Channel-aware front/back polarisation calculation without plotting.

    [~, runGFP, runRFP] = parseChannelMode_local(opts.channelMode);

    if runGFP
        if isempty(distProfileGFP) || size(distProfileGFP,1) ~= opts.no_points
            T = table(); return;
        end
    end

    if runRFP
        if isempty(distProfileRFP) || size(distProfileRFP,1) ~= opts.no_points
            T = table(); return;
        end
    end

    if runGFP && runRFP
        n = min(size(distProfileGFP,2), size(distProfileRFP,2));
        G = distProfileGFP(:,1:n);
        R = distProfileRFP(:,1:n);

        [fG, bG, polG] = computePolarization_local(G, opts.binsFront, opts.binsBack);
        [fR, bR, polR] = computePolarization_local(R, opts.binsFront, opts.binsBack);

        valid = isfinite(fG) & isfinite(bG) & isfinite(polG) & ...
                isfinite(fR) & isfinite(bR) & isfinite(polR);

        fG = fG(valid); bG = bG(valid); polG = polG(valid);
        fR = fR(valid); bR = bR(valid); polR = polR(valid);

        if isempty(polG)
            T = table(); return;
        end

        T = table(fG, bG, polG, fR, bR, polR, ...
            'VariableNames', {'Front_GFP','Back_GFP','Pol_GFP', ...
                              'Front_RFP','Back_RFP','Pol_RFP'});
        return;
    end

    if runGFP
        [fG, bG, polG] = computePolarization_local(distProfileGFP, opts.binsFront, opts.binsBack);
        valid = isfinite(fG) & isfinite(bG) & isfinite(polG);
        fG = fG(valid); bG = bG(valid); polG = polG(valid);

        if isempty(polG)
            T = table(); return;
        end

        T = table(fG, bG, polG, ...
            'VariableNames', {'Front_GFP','Back_GFP','Pol_GFP'});
        return;
    end

    if runRFP
        [fR, bR, polR] = computePolarization_local(distProfileRFP, opts.binsFront, opts.binsBack);
        valid = isfinite(fR) & isfinite(bR) & isfinite(polR);
        fR = fR(valid); bR = bR(valid); polR = polR(valid);

        if isempty(polR)
            T = table(); return;
        end

        T = table(fR, bR, polR, ...
            'VariableNames', {'Front_RFP','Back_RFP','Pol_RFP'});
        return;
    end

    T = table();
end

% ========================================================================
function [frontVals, backVals, polVals] = computePolarization_local(X, binsFront, binsBack)
    frontVals = mean(X(binsFront,:), 1, 'omitnan')';
    backVals  = mean(X(binsBack,:),  1, 'omitnan')';

    valid = isfinite(frontVals) & isfinite(backVals) & (backVals ~= 0);
    frontVals = frontVals(valid);
    backVals  = backVals(valid);

    polVals = frontVals ./ backVals;
end

% ========================================================================
function makePooledStageSummaryPlot(outRoot, tpKey, cond, condName, shortLb, gPol, rPol, opts)

    [~, runGFP, runRFP] = parseChannelMode_local(opts.channelMode);

    nStages = numel(shortLb);
    colG = [0 0.6 0];
    colR = [0.8 0 0.8];
    medCol = [0.10 0.10 0.10];

    f = figure('Color','w');
    hold on;

    dx = 0.18;
    jitter = 0.045;

    for s = 1:nStages
        xG = s - dx;
        xR = s + dx;

        if runGFP
            g = gPol{s};
            if ~isempty(g)
                scatter(xG + jitter*randn(size(g)), g, 14, 'filled', ...
                    'MarkerFaceColor', colG, 'MarkerEdgeColor', 'none', ...
                    'MarkerFaceAlpha', 0.25);
                m = median(g, 'omitnan');
                line([xG-0.08 xG+0.08], [m m], 'Color', medCol, 'LineWidth', 2.2);
                text(xG, m, sprintf('  n=%d', numel(g)), ...
                    'FontSize', 8, 'VerticalAlignment','bottom');
            end
        end

        if runRFP
            r = rPol{s};
            if ~isempty(r)
                scatter(xR + jitter*randn(size(r)), r, 14, 'filled', ...
                    'MarkerFaceColor', colR, 'MarkerEdgeColor', 'none', ...
                    'MarkerFaceAlpha', 0.25);
                m = median(r, 'omitnan');
                line([xR-0.08 xR+0.08], [m m], 'Color', medCol, 'LineWidth', 2.2);
                text(xR, m, sprintf('  n=%d', numel(r)), ...
                    'FontSize', 8, 'VerticalAlignment','bottom');
            end
        end
    end

    h = gobjects(0);
    labels = {};
    if runGFP
        h(end+1) = plot(nan,nan,'o','MarkerFaceColor',colG,'MarkerEdgeColor','none'); %#ok<AGROW>
        labels{end+1} = 'GFP'; %#ok<AGROW>
    end
    if runRFP
        h(end+1) = plot(nan,nan,'o','MarkerFaceColor',colR,'MarkerEdgeColor','none'); %#ok<AGROW>
        labels{end+1} = 'RFP'; %#ok<AGROW>
    end
    if ~isempty(h)
        legend(h, labels, 'Location','northeastoutside');
    end

    xlim([0.5 nStages+0.5]);
    xticks(1:nStages);
    xticklabels(shortLb);
    ylabel('Polarisation: front mean / back mean');
    [modeTag, modeTitle] = intensityModeTagTitle_local(opts);
    title(sprintf('Pooled polarisation by stage — %s — %s — %s', ...
        condName, tpKey, modeTitle), ...
        'FontWeight','bold', 'Interpreter','none');
    grid on; box on; axis square;
    set(gca,'FontSize',12,'LineWidth',1);

    if isfield(opts,'ylimSummary') && ~isempty(opts.ylimSummary)
        ylim(opts.ylimSummary);
    elseif isfield(opts,'ylim') && ~isempty(opts.ylim)
        ylim(opts.ylim);
    end

    outDir = fullfile(outRoot, 'SummaryPlots');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    stem = sprintf('PooledStageSummary_%s_%s_%s_%s', cleanTag(cond), cleanTag(tpKey), cleanTag(opts.channelMode), cleanTag(modeTag));
    saveas(f, fullfile(outDir, [stem '.png']));
    saveas(f, fullfile(outDir, [stem '.fig']));
    close(f);
end

% ========================================================================
function makePooledTimecoursePlot(outRoot, pooledResults, shortLb, opts)

    [~, runGFP, runRFP] = parseChannelMode_local(opts.channelMode);

    timepoints = string(opts.timeKeys);
    conditions = string(opts.conditions);
    nT = numel(timepoints);
    nS = numel(shortLb);

    f = figure('Color','w');
    tl = tiledlayout(f, 1, nT, 'Padding','compact', 'TileSpacing','compact');

    dx = 0.18;
    jitter = 0.045;
    colVeh = [0 0.6 0];
    colH2  = [0.1 0.1 0.1];
    medCol = [0.10 0.10 0.10];

    for t = 1:nT
        nexttile;
        hold on;

        tp = timepoints(t);

        for s = 1:nS
            stageLabel = shortLb(s);

            for c = 1:min(2, numel(conditions))
                cond = conditions(c);

                if c == 1
                    x = s - dx;
                    col = colVeh;
                else
                    x = s + dx;
                    col = colH2;
                end

                idx = pooledResults.Timepoint == tp & pooledResults.Condition == cond;

                % Match by short stage label using StageFolder.
                switch stageLabel
                    case "Ht"
                        idx = idx & pooledResults.StageFolder == "Horsetail";
                    case "M1"
                        idx = idx & pooledResults.StageFolder == "M1";
                    case "M2"
                        idx = idx & pooledResults.StageFolder == "M2";
                    case "LM2"
                        idx = idx & pooledResults.StageFolder == "M2l";
                    case "S"
                        idx = idx & pooledResults.StageFolder == "Spores";
                end

                if runGFP && ismember('Pol_GFP', pooledResults.Properties.VariableNames)
                    vals = pooledResults.Pol_GFP(idx);
                elseif runRFP && ismember('Pol_RFP', pooledResults.Properties.VariableNames)
                    vals = pooledResults.Pol_RFP(idx);
                else
                    vals = [];
                end

                vals = vals(isfinite(vals));
                if isempty(vals), continue; end

                scatter(x + jitter*randn(size(vals)), vals, 12, 'filled', ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor','none', ...
                    'MarkerFaceAlpha', 0.25);

                m = median(vals, 'omitnan');
                line([x-0.07 x+0.07], [m m], 'Color', medCol, 'LineWidth', 2.2);
            end
        end

        xlim([0.5 nS+0.5]);
        xticks(1:nS);
        xticklabels(shortLb);
        grid on; box on;
        set(gca,'FontSize',11,'LineWidth',1);
        title(tp, 'FontWeight','bold', 'Interpreter','none');

        if t == 1
            if runGFP
                ylabel('GFP polarisation: front/back');
            elseif runRFP
                ylabel('RFP polarisation: front/back');
            end
        end

        if isfield(opts,'ylimTimecourse') && ~isempty(opts.ylimTimecourse)
            ylim(opts.ylimTimecourse);
        elseif isfield(opts,'ylim') && ~isempty(opts.ylim)
            ylim(opts.ylim);
        end
    end

    [modeTag, modeTitle] = intensityModeTagTitle_local(opts);
    title(tl, sprintf('Pooled polarisation timecourse (%s)', modeTitle), ...
        'FontWeight','bold', 'FontSize',14);

    h1 = plot(nan,nan,'o','MarkerFaceColor',colVeh,'MarkerEdgeColor','none');
    h2 = plot(nan,nan,'o','MarkerFaceColor',colH2,'MarkerEdgeColor','none');
    legend([h1 h2], cellstr(conditions(1:min(2,numel(conditions)))), ...
        'Location','southoutside', 'Orientation','horizontal');

    outDir = fullfile(outRoot, 'SummaryPlots');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    stem = sprintf('PooledTimecourse_%s_%s', cleanTag(opts.channelMode), cleanTag(modeTag));
    saveas(f, fullfile(outDir, [stem '.png']));
    saveas(f, fullfile(outDir, [stem '.fig']));
    close(f);
end

% ========================================================================
function [modeTag, modeTitle] = intensityModeTagTitle_local(opts)
    if isfield(opts,'intensityMode')
        modeTag = lower(string(opts.intensityMode));
    else
        modeTag = "normalized";
    end

    if modeTag == "normalised" || modeTag == "norm"
        modeTag = "normalized";
    end

    if modeTag == "raw"
        modeTitle = 'raw intensity';
    else
        modeTag = "normalized";
        modeTitle = 'min-max normalized intensity';
    end
end

% ========================================================================
function yLabel = profileYLabel_local(opts, channelName)
    [modeTag, ~] = intensityModeTagTitle_local(opts);
    channelName = char(string(channelName));

    if modeTag == "raw"
        yLabel = sprintf('Raw %s intensity (a.u.)', channelName);
    else
        yLabel = sprintf('Normalized %s intensity ((x-min)/range)', channelName);
    end
end

% ========================================================================
function applyProfileYLim_local(profileCells, opts)
    [modeTag, ~] = intensityModeTagTitle_local(opts);

    if isfield(opts,'profileYLim') && ~isempty(opts.profileYLim)
        ylim(opts.profileYLim);
        return;
    end

    if modeTag ~= "raw"
        ylim([0 1]);
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
        ylim([0 1]);
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

    ylim(yL);
end

% ========================================================================
function [channelMode, runGFP, runRFP] = parseChannelMode_local(channelModeIn)
    channelMode = lower(string(channelModeIn));
    if ~ismember(channelMode, ["both","gfp","rfp"])
        error("channelMode must be 'both', 'gfp', or 'rfp'.");
    end
    runGFP = any(channelMode == ["both","gfp"]);
    runRFP = any(channelMode == ["both","rfp"]);
end

% ========================================================================
function tag = cleanTag(x)
    tag = char(string(x));
    tag = strrep(tag, ' ', '_');
    tag = regexprep(tag, '[^\w\-]', '_');
end
