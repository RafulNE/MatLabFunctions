function analyse_H2O2perturbation(masterCsvOrDir, opts)
% analyse_H2O2perturbation
%
% Comprehensive analysis and plotting pipeline for H2O2-induced
% mitochondrial perturbation experiments using a configurable mitochondrial
% dye, with optional GFP-based ratio analysis.
%
% The function supports multi-experiment meiotic datasets, stage-free
% datasets such as interphase experiments, and simplified legacy datasets
% containing dye Mean and Integrated Density measurements.
%
% -------------------------------------------------------------------------
% OVERVIEW
% -------------------------------------------------------------------------
% This function loads one or multiple Fiji-generated ROI measurement tables,
% harmonizes table schemas across experiments, parses experimental factors
% from filenames, applies quality-control filtering, and generates raw and
% Vehicle-normalised analyses of mitochondrial fluorescence measurements.
%
% In full analysis mode, the function also calculates the ratio:
%
%       GFP_RawIntDen / <dye>_RawIntDen
%
% Data can optionally be stratified by meiotic stage. Statistical inference
% is performed using per-experiment medians as the unit of replication,
% while plots retain individual ROI/cell values, experiment-level medians,
% and pooled medians.
%
% Two plotting layouts are available:
%
%       "single"  - all timepoints and treatment groups on one axis
%       "panels"  - separate timepoint panels
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% masterCsvOrDir : (string)
%     Path to either:
%       - a single CSV file, or
%       - a directory containing multiple CSV files.
%
% opts : Name-value options defined in the arguments block.
%
% -------------------------------------------------------------------------
% CORE FEATURES
% -------------------------------------------------------------------------
%
% 1) MULTI-EXPERIMENT INPUT AND TABLE HARMONISATION
%    - Accepts a single CSV or multiple CSV files from a folder.
%    - Sorts CSV paths for reproducible experiment ordering.
%    - Harmonises columns before concatenating tables.
%    - Automatically assigns Experiment labels.
%    - If a CSV filename begins with an 8-digit date (YYYYMMDD), the
%      experiment is labelled:
%
%          Exp_YYYYMMDD
%
%      Otherwise, labels fall back to Exp1, Exp2, Exp3, etc.
%
% 2) DYE-AGNOSTIC ANALYSIS
%    - The mitochondrial dye is selected using opts.dyeName.
%    - Examples include:
%
%          "MTxRos"
%          "TMRE"
%
%    - Full analysis mode expects:
%
%          <dye>_Mean
%          <dye>_IntDen
%          <dye>_RawIntDen
%
%    - Full mode additionally calculates:
%
%          Ratio_GFPraw_over_<dye>_RawIntDen
%
% 3) FULL AND LEGACY / MINIMAL ANALYSIS MODES
%    - Full mode analyses GFP and mitochondrial dye measurements.
%    - Legacy/minimal mode is enabled with:
%
%          minimalDyeOnly=true
%
%    - Minimal mode analyses only:
%
%          <dye>_Mean
%          <dye>_IntDen
%
%    - Existing normalised dye columns can optionally be reused with:
%
%          useProvidedNorm=true
%
%    - Legacy datasets can be forced into a single timepoint and manually
%      assigned Vehicle and treated groups using filename-based matching.
%
% 4) AUTOMATIC FILENAME PARSING
%    - Parses the File column to determine:
%
%          Condition   : Vehicle, H2O2, or Unparsed
%          Dose_mM     : H2O2 concentration when present
%          Time_min    : acquisition time
%          Series      : series number
%
%    - Supports filenames containing additional prefixes such as genotype
%      or experiment date.
%    - Vehicle filenames containing "Veh" or "Vehicle" are recognised.
%    - H2O2 doses are extracted from explicit mM labels.
%    - Times expressed in minutes or hours are recognised.
%    - Files without an explicit time token default internally to Time_min=0.
%
% 5) QUALITY CONTROL AND STAGE FILTERING
%    - Optionally removes non-finite measurements.
%    - Full mode supports minimum GFP and dye RawIntDen thresholds.
%    - Stage exclusions are applied separately after measurement QC.
%    - The number of retained and removed rows is reported to the console.
%    - A parsed pre-QC table is saved for filename-parsing inspection.
%
% 6) OPTIONAL TIMEPOINT SELECTION
%    - opts.timePoints can restrict the analysis to selected Time_min values.
%    - Timepoint filtering is applied before:
%
%          ROI counting
%          normalisation
%          statistics
%          plotting
%
%    - The requested timepoint order is also respected in plots.
%
% 7) OPTIONAL STAGE STRATIFICATION
%    - Enabled with:
%
%          byStage=true
%
%    - Supports ordered meiotic stages including:
%
%          Horsetail
%          Meiosis I
%          Meiosis II
%          Late Meiosis II
%          Spores
%
%    - Additional stage categories are retained and appended automatically.
%    - Selected stages can be excluded using opts.excludeStages.
%    - Blank, NaN, or zero stage entries trigger an error when byStage=true.
%
%    - When byStage=false, a dummy stage label (default "ALL") is generated
%      internally so the same statistical pipeline can analyse stage-free
%      datasets such as interphase experiments.
%
% 8) VEHICLE NORMALISATION
%    - By default, each metric is divided by the median Vehicle value from
%      the corresponding timepoint.
%
%    - With:
%
%          normByStage=true
%
%      and:
%
%          byStage=true
%
%      normalisation is instead performed independently for each
%      Time_min × Stage combination.
%
%    - Standard output variables are:
%
%          Norm_<dye>_Mean
%          Norm_<dye>_IntDen
%          Norm_Ratio_GFPraw_over_<dye>_RawIntDen
%
%    - In minimal mode, existing <dye>_Mean_Norm and <dye>_IntDen_Norm
%      columns can be reused instead of recalculating normalisation.
%
% 9) ROI COUNT SUMMARIES
%    - Saves retained ROI/cell counts after QC and stage filtering.
%    - With byStage=true, counts are summarised by:
%
%          Time × Stage × PlotGroup
%
%      and:
%
%          Time × Stage × PlotGroup × Experiment
%
%    - Without stage stratification, counts are summarised by:
%
%          Time × PlotGroup
%
% 10) VISUALISATION
%    - Generates plots for raw and normalised metrics.
%    - Individual ROI/cell values are shown as jittered points.
%    - Experiments are distinguished using fixed colour-lightening factors.
%    - Per-experiment medians are shown as large outlined points.
%    - Pooled ROI/cell medians are shown as horizontal black lines.
%
%    - opts.plotLayout controls the figure structure:
%
%          "single"  - grouped single-axis layout
%          "panels"  - separate timepoint panels
%
%    - The single-axis layout supports:
%
%          time-only or full x-axis labels
%          short or compact time labels
%          optional x-axis tick removal
%          optional grid, legend, and figure title
%          configurable figure dimensions
%          configurable tick and axis-label font sizes
%          configurable point and median sizes
%          configurable axis and median line widths
%
%    - Figures are automatically exported in opts.figFormat when
%      opts.saveFigs=true.
%
% 11) PLOT-ONLY OUTLIER REMOVAL
%    - Optional outlier filtering can be applied only to plotted values.
%    - Supported methods are:
%
%          "percentile"
%          "iqr"
%          "absolute"
%
%    - Plot-only outlier removal does NOT modify:
%
%          the parsed analysis table
%          ROI counts
%          normalisation
%          statistical analysis
%
%    - The number of plotted values retained and the applied limits are
%      reported to the console.
%
% 12) EXPERIMENT-LEVEL NONPARAMETRIC STATISTICS
%    - Enabled with:
%
%          doStats=true
%
%    - The statistical unit of replication is the Experiment median.
%    - Analyses are performed per:
%
%          Time × Stage
%
%      or per Time when byStage=false and stage-free statistics are allowed.
%
%    - For exactly two treatment groups:
%
%          Wilcoxon rank-sum test
%
%    - For more than two treatment groups:
%
%          Kruskal-Wallis overall test
%          Dunn-Sidak post hoc pairwise comparisons
%
%      when the MATLAB multiple-comparison procedure is available.
%
%    - Vehicle-versus-dose comparisons are extracted for:
%
%          Vehicle vs H2O2 10 mM
%          Vehicle vs H2O2 50 mM
%          Vehicle vs H2O2 100 mM
%
%    - Spearman correlation is used to test for a dose-dependent trend.
%    - Significant Vehicle-versus-dose comparisons can be annotated with
%      stars on normalised plots.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% The function can generate:
%
%   - <fileStem>_parsed_beforeQC.csv
%       Parsed experimental factors before QC filtering.
%
%   - ROI count summary CSV files.
%
%   - Raw metric figures.
%
%   - Vehicle-normalised metric figures.
%
%   - <fileStem>_STATS_summary_<metric>.csv
%       Overall tests, dose trends, experiment counts, and Vehicle-versus-
%       dose p-values.
%
%   - <fileStem>_STATS_pairwise_<metric>.csv
%       Pairwise treatment comparisons.
%
%   - <fileStem>_parsed_withRatios_andNorm.csv
%       Final filtered table containing parsed factors, derived ratios, and
%       normalised measurements.
%
%   - Console reports describing experiment mapping, QC filtering,
%       timepoint selection, plot-only outlier filtering, and Vehicle
%       normalisation values.
%
% -------------------------------------------------------------------------
% DESIGN PRINCIPLES
% -------------------------------------------------------------------------
% - Statistical inference uses experiment-level medians rather than treating
%   individual ROIs/cells as independent biological replicates.
% - Individual ROI/cell distributions remain visible in the figures.
% - Dye naming is configurable and shared across metric generation,
%   normalisation, figure labelling, and output filenames.
% - The same analysis framework supports stage-stratified meiotic datasets
%   and stage-free datasets.
% - Legacy dataset overrides are retained for compatibility with older,
%   simplified Fiji measurement tables.
% - Plot-only filtering is isolated from quantitative analysis and statistics.
%
% -------------------------------------------------------------------------
% EXAMPLE USAGE: FULL DATASET
% -------------------------------------------------------------------------
%
% analyse_H2O2perturbation( ...
%     "path/to/RAW Tables", ...
%     outDir="path/to/output", ...
%     dyeName="TMRE", ...
%     timePoints=[0 60], ...
%     byStage=true, ...
%     normByStage=false, ...
%     plotLayout="single", ...
%     timeLabelMode="short", ...
%     doStats=true, ...
%     saveFigs=true ...
% );
%
% -------------------------------------------------------------------------
% EXAMPLE USAGE: LEGACY / MINIMAL DATASET
% -------------------------------------------------------------------------
%
% analyse_H2O2perturbation( ...
%     "path/to/legacy_dataset.csv", ...
%     outDir="path/to/output", ...
%     dyeName="MTxRos", ...
%     minimalDyeOnly=true, ...
%     useProvidedNorm=true, ...
%     forceSingleTime=true, ...
%     singleTime_min=0, ...
%     forcePlotGroups=true, ...
%     treatedKey="H2O2", ...
%     treatedLabel="H2O2 10 mM", ...
%     byStage=true, ...
%     doStats=true, ...
%     saveFigs=true ...
% );
%
% -------------------------------------------------------------------------
% Author: Raful Navarro Espindola
% Project: Mitochondrial dynamics during meiotic oxidative perturbation
% -------------------------------------------------------------------------
arguments
    masterCsvOrDir (1,1) string

    opts.csvPattern (1,1) string = "*.csv"     % used if you pass a folder
    opts.expNameFromFile (1,1) logical = true  % Exp label from csv filename
    opts.outDir (1,1) string = ""
    opts.saveFigs (1,1) logical = true
    opts.figFormat (1,1) string = "png"
    opts.fileStem (1,1) string = ""
    opts.dyeName (1,1) string = "MTxRos"   % e.g. "MTxRos" or "TMRE"

    opts.timePoints (1,:) double = []
    opts.timeZeroLabel (1,1) string = "5 min"
    opts.plotLayout (1,1) string = "single"
    opts.showGrid (1,1) logical = false
    opts.showLegend (1,1) logical = false
    opts.showFigureTitle (1,1) logical = false
    opts.timeLabelMode (1,1) string = "short"  % "short" or "compact"

    % X-axis label style
    opts.xTickLabelMode (1,1) string = "timeOnly"   % "timeOnly" or "full"

    % Plot-only outlier removal
    opts.removePlotOutliers (1,1) logical = false
    opts.plotOutlierMethod (1,1) string = "percentile"  % "percentile", "iqr", or "absolute"
    opts.plotOutlierPercentiles (1,2) double = [1 99]
    opts.plotOutlierIQRMultiplier (1,1) double = 3
    opts.plotOutlierLimits (1,2) double = [-Inf Inf]

    % Figure formatting
    opts.figSizeCm (1,2) double = [42 21.0]   % A4 landscape-friendly
    opts.xTickAngle (1,1) double = 0
    opts.xTickFontSize (1,1) double = 16
    opts.yTickFontSize (1,1) double = 24
    opts.yLabelFontSize (1,1) double = 30
    opts.yLabelFontWeight (1,1) string = "bold"
    opts.axisLineWidth (1,1) double = 1.5
    opts.hideXAxisTicks (1,1) logical = false
    opts.singlePointSize (1,1) double = 60
    opts.experimentMedianSize (1,1) double = 210
    opts.pooledMedianLineWidth (1,1) double = 1.65
    opts.pooledMedianHalfWidth (1,1) double = 0.35

    % Stage stratification
    opts.byStage (1,1) logical = false
    opts.stageColumn (1,1) string = "Stage"
    opts.stageOrder (1,:) string = ["Horsetail","Meiosis I","Meiosis II","Late Meiosis II","Spores","UNSTAGED","DELETE"]
    opts.excludeStages (1,:) string = ["DELETE"]
    opts.stageShortLabels (1,1) logical = true
    opts.expShadeMin (1,1) double = 0.10
    opts.expShadeMax (1,1) double = 0.80

    % QC
    opts.excludeNonFinite (1,1) logical = true
    opts.minGFPraw (1,1) double = 0
    opts.minDyeRaw (1,1) double = 0

    % Stats
    opts.doStats (1,1) logical = true
    opts.statsAlpha (1,1) double = 0.05
    opts.annotateStats (1,1) logical = true
    opts.statsYOffsetFrac (1,1) double = 0.06   % vertical offset for stars (fraction of y-range)
    opts.allowStatsWithoutStage (1,1) logical = true
    opts.noStageLabel (1,1) string = "ALL"
    opts.normByStage (1,1) logical = false

    % -------------------- Legacy / Minimal mode --------------------
    opts.minimalDyeOnly (1,1) logical = false   % analyse only dye Mean/IntDen
    opts.useProvidedNorm (1,1) logical = true   % use *_Norm columns if present

    % Force grouping (for old datasets)
    opts.forceSingleTime (1,1) logical = false
    opts.singleTime_min (1,1) double  = 0

    opts.forcePlotGroups (1,1) logical = false
    opts.vehicleLabel (1,1) string = "Vehicle"
    opts.treatedLabel (1,1) string = "H2O2 10 mM"
    opts.treatedKey (1,1) string = "H2O2"   % matches filenames like H2O2_series2

end

dye = opts.dyeName;

% Validate plot layout option
opts.plotLayout = lower(opts.plotLayout);

if ~ismember(opts.plotLayout, ["single","panels"])
    error("opts.plotLayout must be either 'single' or 'panels'.");
end

opts.xTickLabelMode = lower(opts.xTickLabelMode);
opts.plotOutlierMethod = lower(opts.plotOutlierMethod);

if ~ismember(opts.xTickLabelMode, ["timeonly","full"])
    error("opts.xTickLabelMode must be either 'timeOnly' or 'full'.");
end

if ~ismember(opts.plotOutlierMethod, ["percentile","iqr","absolute"])
    error("opts.plotOutlierMethod must be 'percentile', 'iqr', or 'absolute'.");
end

opts.timeLabelMode = lower(opts.timeLabelMode);

if ~ismember(opts.timeLabelMode, ["short","compact"])
    error("opts.timeLabelMode must be either 'short' or 'compact'.");
end

% Default output stem follows the selected dye unless explicitly provided
if opts.fileStem == ""
    opts.fileStem = "H2O2_" + dye;
end

% -------------------- IO (file OR folder) --------------------
if ~(isfile(masterCsvOrDir) || isfolder(masterCsvOrDir))
    error("Input must be a CSV file or a folder of CSVs: %s", masterCsvOrDir);
end

% Default outDir
if opts.outDir == ""
    if isfolder(masterCsvOrDir)
        opts.outDir = masterCsvOrDir;
    else
        opts.outDir = fileparts(masterCsvOrDir);
    end
end
if opts.outDir == ""
    opts.outDir = pwd;
end
if ~isfolder(opts.outDir)
    mkdir(opts.outDir);
end

% -------------------- Load one or many CSVs --------------------
csvPaths = strings(0,1);

if isfolder(masterCsvOrDir)
    D = dir(fullfile(masterCsvOrDir, opts.csvPattern));
    if isempty(D)
        error("No CSVs found in folder: %s (pattern: %s)", masterCsvOrDir, opts.csvPattern);
    end
    csvPaths = string(fullfile({D.folder}, {D.name}))';
else
    if ~isfile(masterCsvOrDir)
        error("Input is neither a folder nor a file: %s", masterCsvOrDir);
    end
    csvPaths = masterCsvOrDir;
end

% Reproducible experiment ordering
csvPaths = sort(csvPaths(:));

Tall = table();
for k = 1:numel(csvPaths)
    Tk = readtable(csvPaths(k), "TextType","string");

    % Add Experiment label
    if opts.expNameFromFile
        [~, baseName] = fileparts(csvPaths(k));

        % Prefer leading date if present, e.g. 20260121_Manual_ROI...
        dateTok = regexp(baseName, "^\d{8}", "match", "once");

        if ~isempty(dateTok)
            expName = "Exp_" + string(dateTok);
        else
            expName = "Exp" + string(k);
        end
    else
        expName = "Exp" + string(k);
    end

    Tk.Experiment = repmat(expName, height(Tk), 1);
    % --- Harmonize columns across CSVs (handles Stage added later) ---
    if isempty(Tall)
        % First file sets the baseline schema
        Tall = Tk;
    else
        % Add any columns missing in Tk that exist in Tall
        missInTk = setdiff(Tall.Properties.VariableNames, Tk.Properties.VariableNames);
        for m = missInTk
            v = m{1};
            TallClass = class(Tall.(v));
            Tk.(v) = defaultColumn(height(Tk), TallClass);
        end

        % Add any columns missing in Tall that exist in Tk (e.g., new Stage column)
        missInTall = setdiff(Tk.Properties.VariableNames, Tall.Properties.VariableNames);
        for m = missInTall
            v = m{1};
            TkClass = class(Tk.(v));
            Tall.(v) = defaultColumn(height(Tall), TkClass);
        end

        % Reorder Tk to match Tall column order
        Tk = Tk(:, Tall.Properties.VariableNames);

        % Now safe to append
        Tall = [Tall; Tk];
    end
end
T = Tall;
T.Experiment = categorical(string(T.Experiment));

fprintf("Loaded %d CSV(s), total rows: %d\n", numel(csvPaths), height(T));

fprintf("\nExperiment mapping:\n");
for k = 1:numel(csvPaths)
    fprintf("  Exp%d → %s\n", k, csvPaths(k));
end

vars = string(T.Properties.VariableNames);

if opts.minimalDyeOnly
    % ---- LEGACY / MINIMAL REQUIREMENTS ----
    required = ["File","ROI_Index", opts.stageColumn, ...
                dye + "_Mean", dye + "_IntDen"];

    % If you want to use provided norm, require those too
    if opts.useProvidedNorm
        required = [required, dye + "_Mean_Norm", dye + "_IntDen_Norm"];
    end

else
    % ---- FULL PIPELINE REQUIREMENTS ----
    required = ["File","ROI_Index","Area_um2","GFP_Mean","GFP_IntDen","GFP_RawIntDen", ...
                dye + "_Mean", dye + "_IntDen", dye + "_RawIntDen"];
end

missing = setdiff(required, vars);
if ~isempty(missing)
    error("Missing required columns: %s", strjoin(missing, ", "));
end

% -------------------- Stage column (only when byStage=true) --------------------
if opts.byStage
    stCol = opts.stageColumn;

    if ~ismember(stCol, string(T.Properties.VariableNames))
        error("byStage=true but column '%s' not found in CSV(s).", stCol);
    end

    % Normalize to clean strings (NO UNSTAGED fallback)
    S = string(T.(stCol));
    S = strtrim(S);

    % Hard fail if any empty/NaN/0 entries (since you prefer toggling byStage)
    bad = (S=="" | S=="NaN" | S=="0" | ismissing(S));
    if any(bad)
        error("Stage column '%s' contains blank/NaN/0 entries. Fix in Fiji or run with byStage=false.", stCol);
    end

    % Categorical with your preferred order (extras appended so nothing breaks)
    present = unique(S);
    pref = opts.stageOrder(ismember(opts.stageOrder, present));
    extras = present(~ismember(present, opts.stageOrder));
    stageCats = [pref(:); sort(extras(:))];

    T.Stage = categorical(S, stageCats);
end

% -------------------- Dummy stage for no-stage / interphase datasets --------------------
% This allows the same stats function to run when byStage=false.
% The stats will then be computed per Time_min only, with Stage = "ALL".
if ~opts.byStage
    T.Stage = categorical(repmat(opts.noStageLabel, height(T), 1), opts.noStageLabel);
end

% -------------------- Parse filename into factors --------------------
% Robust parsing of:
%   Condition: Veh vs H2O2
%   Dose_mM: 0 for Vehicle, numeric for H2O2 when present
%   Time_min: defaults to 0 if absent
%   Series: numeric when present
%
% This parser is token/substring-based, so it tolerates filenames such as:
%   H2O2_100mM_60min_series1.tif
%   Fzo1D_H2O2_60min_series1.tif
%   20260121_H2O2_100mM_60min_series1.tif
%   Veh_series1.tif
%   Vehicle_30min_series2.tif

[Condition, Dose_mM, Time_min, Series] = parseH2O2Filename(string(T.File));

T.Condition = categorical(Condition, ["Veh","H2O2","Unparsed"]);
T.Dose_mM   = Dose_mM;
T.Time_min  = Time_min;
T.Series    = Series;

% -------------------- Derived groups + ratio --------------------
% Group label for plotting: Vehicle vs H2O2 dose (within time)
plotGroup = strings(height(T),1);

for i = 1:height(T)
    if T.Condition(i) == "Veh"
        plotGroup(i) = opts.vehicleLabel;

    elseif T.Condition(i) == "H2O2"
        if isfinite(T.Dose_mM(i))
            plotGroup(i) = "H2O2 " + string(T.Dose_mM(i)) + " mM";
        else
            plotGroup(i) = "H2O2 unknown dose";
        end

    else
        plotGroup(i) = "Unparsed";
    end
end

T.PlotGroup = categorical(plotGroup);

% Explicit group order (within each timepoint facet)
groupOrder = [opts.vehicleLabel, ...
              "H2O2 10 mM", ...
              "H2O2 50 mM", ...
              "H2O2 100 mM", ...
              "H2O2 unknown dose", ...
              "Unparsed"];

existingGroups = string(categories(T.PlotGroup));

present = groupOrder(ismember(groupOrder, existingGroups));
extras  = existingGroups(~ismember(existingGroups, groupOrder));

present = present(:);
extras  = extras(:);

T.PlotGroup = reordercats(T.PlotGroup, cellstr([present; extras]));

% Ratio requested: GFP_RawIntDen / <dye>_RawIntDen
meanVar = dye + "_Mean";
intDenVar = dye + "_IntDen";
rawDyeVar = dye + "_RawIntDen";

if ~opts.minimalDyeOnly
    ratioVarName = "Ratio_GFPraw_over_" + dye + "_RawIntDen";
    T.(ratioVarName) = T.GFP_RawIntDen ./ T.(rawDyeVar);
end

% Save parsed table before QC filtering.
% This is useful for checking whether filenames were parsed correctly.
outParsedBeforeQC = fullfile(opts.outDir, opts.fileStem + "_parsed_beforeQC.csv");
writetable(T, outParsedBeforeQC);
fprintf("Saved parsed-before-QC table:\n  %s\n", outParsedBeforeQC);

% -------------------- QC filtering (optional) --------------------
keep = true(height(T),1);

if opts.minimalDyeOnly
    % ---------- LEGACY / MINIMAL QC ----------
    
    if opts.excludeNonFinite
        keep = keep & ...
            isfinite(T.(meanVar)) & ...
            isfinite(T.(intDenVar));
        
        if opts.useProvidedNorm
            legacyNormMean   = meanVar + "_Norm";
            legacyNormIntDen = intDenVar + "_Norm";
            
            if ismember(legacyNormMean, string(T.Properties.VariableNames))
                keep = keep & isfinite(T.(legacyNormMean));
            end
            if ismember(legacyNormIntDen, string(T.Properties.VariableNames))
                keep = keep & isfinite(T.(legacyNormIntDen));
            end
        end
    end

else
    % ---------- FULL PIPELINE QC ----------
    
    keep = keep & (T.GFP_RawIntDen >= opts.minGFPraw);
    keep = keep & (T.(dye + "_RawIntDen") >= opts.minDyeRaw);

    if opts.excludeNonFinite
        keep = keep & ...
            isfinite(T.GFP_RawIntDen) & ...
            isfinite(T.(dye + "_RawIntDen")) & ...
            isfinite(T.(meanVar)) & ...
            isfinite(T.(intDenVar));
        
        if exist("ratioVarName","var")
            keep = keep & isfinite(T.(ratioVarName));
        end
    end
end

% Apply QC filtering
T_beforeFiltering = T;
T = T(keep, :);

% Apply stage filtering separately after QC
if opts.byStage && ~isempty(opts.excludeStages)
    stageKeep = ~ismember(string(T.Stage), opts.excludeStages);
    T = T(stageKeep, :);
end

fprintf("Loaded %d rows; kept %d after QC/stage filtering (removed %d).\n", ...
    height(T_beforeFiltering), height(T), height(T_beforeFiltering) - height(T));

% -------------------- Legacy overrides (single timepoint / single dose) --------------------
if opts.forceSingleTime
    T.Time_min = repmat(opts.singleTime_min, height(T), 1);
end

if opts.forcePlotGroups
    pg = repmat(opts.vehicleLabel, height(T), 1);

    if ismember("File", T.Properties.VariableNames)
        isTreat = startsWith(string(T.File), opts.treatedKey, 'IgnoreCase', true);
        pg(isTreat) = opts.treatedLabel;
    else
        error("Cannot force groups: 'File' column not found.");
    end

    T.PlotGroup = categorical(pg, [opts.vehicleLabel, opts.treatedLabel]);

    % Compatibility fields (used by normalisation and printing)
    T.Condition = categorical(pg == opts.treatedLabel, [false true], ["Veh","H2O2"]);
    T.Dose_mM = nan(height(T),1);
    T.Dose_mM(T.Condition=="Veh") = 0;
    T.Dose_mM(T.Condition=="H2O2") = 10; % descriptive only
end

% -------------------- Optional timepoint filtering --------------------
% This affects counts, normalization, statistics, and plotting.
if ~isempty(opts.timePoints)
    wantedTimes = opts.timePoints(:)';

    nBeforeTimeFilter = height(T);
    T = T(ismember(T.Time_min, wantedTimes), :);

    if iscategorical(T.PlotGroup)
        T.PlotGroup = removecats(T.PlotGroup);
    end

    if iscategorical(T.Condition)
        T.Condition = removecats(T.Condition);
    end

    if opts.byStage && iscategorical(T.Stage)
        T.Stage = removecats(T.Stage);
    end

    fprintf("Timepoint filter: kept %d of %d rows using Time_min = [%s].\n", ...
        height(T), nBeforeTimeFilter, num2str(wantedTimes));

    if isempty(T)
        error("No rows remain after timepoint filtering. Check opts.timePoints and parsed Time_min values.");
    end
end

% -------------------- Counts summary (N) --------------------
% N = number of ROIs/cells kept after QC (and stage filtering if enabled)

if opts.byStage
    % Per time × stage × group (pooled across experiments)
    N_stage_group_time = groupsummary(T, ["Time_min","Stage","PlotGroup"], "IncludeEmptyGroups", true);
    N_stage_group_time.Properties.VariableNames(strcmp(N_stage_group_time.Properties.VariableNames,"GroupCount")) = "N_ROIs";

    outCounts1 = fullfile(opts.outDir, opts.fileStem + "_N_byStage_byGroup_byTime.csv");
    writetable(N_stage_group_time, outCounts1);

    % Optional: per time × stage × group × experiment (useful QC)
    N_stage_group_time_exp = groupsummary(T, ["Time_min","Stage","PlotGroup","Experiment"], "IncludeEmptyGroups", true);
    N_stage_group_time_exp.Properties.VariableNames(strcmp(N_stage_group_time_exp.Properties.VariableNames,"GroupCount")) = "N_ROIs";

    outCounts2 = fullfile(opts.outDir, opts.fileStem + "_N_byStage_byGroup_byTime_byExperiment.csv");
    writetable(N_stage_group_time_exp, outCounts2);

else
    % If stages are off, still save the old-style count table per time × group
    N_group_time = groupsummary(T, ["Time_min","PlotGroup"], "IncludeEmptyGroups", true);
    N_group_time.Properties.VariableNames(strcmp(N_group_time.Properties.VariableNames,"GroupCount")) = "N_ROIs";

    outCounts = fullfile(opts.outDir, opts.fileStem + "_N_byGroup_byTime.csv");
    writetable(N_group_time, outCounts);
end

% -------------------- Normalize to Vehicle median (per timepoint) --------------------
% Standard normalized variable names used downstream
normMeanVar   = "Norm_" + meanVar;
normIntDenVar = "Norm_" + intDenVar;

% Pre-create normalized columns (table-safe)
T.(normMeanVar)   = nan(height(T),1);
T.(normIntDenVar) = nan(height(T),1);

% Ratio norm only exists in full mode
if ~opts.minimalDyeOnly
    normRatioVar = "Norm_" + ratioVarName;
    T.(normRatioVar) = nan(height(T),1);
end

% ---- If legacy dataset already provides normalized columns, just reuse them ----
legacyNormMean   = meanVar + "_Norm";
legacyNormIntDen = intDenVar + "_Norm";

if opts.minimalDyeOnly && opts.useProvidedNorm ...
        && ismember(legacyNormMean, string(T.Properties.VariableNames)) ...
        && ismember(legacyNormIntDen, string(T.Properties.VariableNames))

    T.(normMeanVar)   = T.(legacyNormMean);
    T.(normIntDenVar) = T.(legacyNormIntDen);

else
    % ---- Otherwise compute normalization to Vehicle median per timepoint ----
    timeLevels = unique(T.Time_min(~isnan(T.Time_min)));

    if opts.normByStage && opts.byStage
        stageLevels = categories(T.Stage);
    else
        stageLevels = {"ALL"};
    end

    for ti = 1:numel(timeLevels)
        tval = timeLevels(ti);

        for si = 1:numel(stageLevels)

            if opts.normByStage && opts.byStage
                stName = string(stageLevels{si});
                idxBase = (T.Time_min == tval) & (string(T.Stage) == stName);
                normLabel = sprintf("time %g min, stage %s", tval, stName);
            else
                idxBase = (T.Time_min == tval);
                normLabel = sprintf("time %g min", tval);
            end

            idxVeh = idxBase & (string(T.PlotGroup) == opts.vehicleLabel);

            if ~any(idxVeh)
                warning("No Vehicle data for %s — skipping normalization for this group.", normLabel);
                continue;
            end

            vehMed_Mean   = median(T.(meanVar)(idxVeh),   'omitnan');
            vehMed_IntDen = median(T.(intDenVar)(idxVeh), 'omitnan');

            if isfinite(vehMed_Mean)
                T.(normMeanVar)(idxBase) = T.(meanVar)(idxBase) ./ vehMed_Mean;
            end

            if isfinite(vehMed_IntDen)
                T.(normIntDenVar)(idxBase) = T.(intDenVar)(idxBase) ./ vehMed_IntDen;
            end

            if ~opts.minimalDyeOnly
                vehMed_Ratio = median(T.(ratioVarName)(idxVeh), 'omitnan');

                if isfinite(vehMed_Ratio)
                    T.(normRatioVar)(idxBase) = T.(ratioVarName)(idxBase) ./ vehMed_Ratio;
                end
            end
        end
    end

    fprintf("Norm_%s_Mean NaNs: time0=%d, time60=%d\n", dye, ...
        sum(isnan(T.(normMeanVar)) & T.Time_min==0), ...
        sum(isnan(T.(normMeanVar)) & T.Time_min==60));

    for tval = unique(T.Time_min(~isnan(T.Time_min)))'
        isVeh = (T.Condition=="Veh") & (T.Time_min==tval);

        fprintf("\nTime %g:\n", tval);
        fprintf("  Veh median %s_Mean = %g\n",  dye, median(T.(meanVar)(isVeh), 'omitnan'));
        fprintf("  Veh median %s_IntDen = %g\n", dye, median(T.(intDenVar)(isVeh), 'omitnan'));
        if ~opts.minimalDyeOnly
            fprintf("  Veh median Ratio = %g\n", ...
                median(T.(ratioVarName)(isVeh), 'omitnan'));
        end
    end
end

% -------------------- Color map (distinct treatment hues) --------------------
if ~isfield(opts, "colorMap") || isempty(opts.colorMap)
    cm = containers.Map;

    % Vehicle = purple
    cm("Vehicle") = [0.45 0.20 0.70];

    % H2O2 doses (distinct hues)
    cm("H2O2 10 mM")  = [0.85 0.10 0.10];  % red
    cm("H2O2 50 mM")  = [0.95 0.55 0.10];  % orange
    cm("H2O2 100 mM") = [0.98 0.85 0.15];  % yellow

    opts.colorMap = cm;
end

    function plotMetric(T, varName, yLabel, stem, opts)
        % Chooses between the original panel layout and the new single-axis layout.

        if opts.plotLayout == "single"
            plotSingleGrouped(T, varName, yLabel, stem, opts);
        else
            plotByTime(T, varName, yLabel, stem, opts);
        end

    end

function plotByTime(T, varName, yLabel, stem, opts)

    % ---- Time panels ordering ----
    timeLevels = orderTimeLevels(T, opts);

    fig = figure('Color','w','Position',[100 100 1700 450]);
    tl = tiledlayout(1, numel(timeLevels)+1, ...
    "TileSpacing","compact", ...
    "Padding","tight");

    % ---- Preferred stage order + labels ----
    if opts.byStage
        stagePref = ["Horsetail","Meiosis I","Meiosis II","Late Meiosis II","Spores"];

        stageCats = string(categories(T.Stage));   % should be string, but normalize anyway
        stageCats = stageCats(:);                  % force column

        prefPresent = stagePref(ismember(stagePref, stageCats));
        prefPresent = prefPresent(:);              % force column

        extras = setdiff(stageCats, stagePref);    % extras as string
        extras = string(extras);
        extras = extras(:);                        % force column

        stageLevels = [prefPresent; extras];

        % Optional: remove DELETE if it somehow exists
        stageLevels = stageLevels(stageLevels ~= "DELETE");
    else
        stageLevels = "ALL";
    end

    % ---- Treatment order within each stage ----
    groupCats = string(categories(T.PlotGroup));  % already reordered earlier
    treatLevels = groupCats;                      % e.g. Vehicle, H2O2 10 mM, ...

    % ---- Fixed experiment shade factors (global) ----
    expsAll = categories(T.Experiment);
    nEall = numel(expsAll);

    if nEall <= 1
        aAll = opts.expShadeMin;
    else
        aAll = linspace(opts.expShadeMin, opts.expShadeMax, nEall);
    end

    % helper: get shade factor for a given experiment category
    % (returns scalar a in [0..1], where larger = lighter)
    getA = @(expCat) aAll(find(strcmp(expsAll, expCat), 1, 'first'));

    % Legend placeholders
    legendAx = [];
    hTreat = gobjects(0);
    treatNames = strings(0);
    hExp = gobjects(0);
    expNames = strings(0);


    for ti = 1:numel(timeLevels)
        tval = timeLevels(ti);
        ax = nexttile; hold(ax,'on');

        subT = T(T.Time_min == tval, :);
        if isempty(subT)
            axis(ax,'off'); 
            continue;
        end

        if ~opts.byStage
            % ==================== ORIGINAL LAYOUT (no stage grouping) ====================
            groups = categories(subT.PlotGroup);

            for gi = 1:numel(groups)
                gname = string(groups{gi});
                maskG = subT.PlotGroup == groups{gi};

                % pooled values (overall median line)
                valsG_all = subT.(varName)(maskG);
                valsG_all = valsG_all(isfinite(valsG_all));

                % base color
                base = [0.5 0.5 0.5];
                if isfield(opts,"colorMap") && ~isempty(opts.colorMap) && isKey(opts.colorMap, gname)
                    base = opts.colorMap(gname);
                end

                % experiments contributing
                exps = categories(subT.Experiment(maskG));
                nE = numel(exps);

                % scatter by experiment (background)
                for ei = 1:nE
                    maskE = maskG & (subT.Experiment == exps{ei});
                    vals = subT.(varName)(maskE);
                    vals = vals(isfinite(vals));
                    if isempty(vals), continue; end

                    a = getA(exps{ei});
                    c_bg = base*(1-a) + [1 1 1]*a;

                    xj = gi + 0.18*(rand(size(vals)) - 0.5);
                    s = scatter(ax, xj, vals, 8, 'filled', ...
                        'MarkerFaceColor', c_bg, ...
                        'MarkerEdgeColor', c_bg*0.6, ...
                        'LineWidth', 0.3);
                    s.MarkerFaceAlpha = 0.65;
                end

                % per-experiment medians (big dots)
                for ei = 1:numel(exps)
                    expName = exps{ei};
                    maskE = maskG & (subT.Experiment == expName);
                    valsE = subT.(varName)(maskE);
                    valsE = valsE(isfinite(valsE));
                    if isempty(valsE), continue; end

                    medVal = median(valsE);

                    % get fixed shade factor for this experiment
                    a = getA(expName);
                    c_med = base*(1-a) + [1 1 1]*a;

                    % small horizontal separation between experiments
                    if numel(exps) == 1
                        xMed = gi;
                    else
                        offsets = linspace(-0.10, 0.10, numel(exps));
                        xMed = gi + offsets(ei);
                    end

                    scatter(ax, xMed, medVal, 30, 'filled', ...
                        'MarkerFaceColor', c_med, ...
                        'MarkerEdgeColor', 'k', ...
                        'LineWidth', 1);
                end
                
                % overall median line
                if ~isempty(valsG_all)
                    medv = median(valsG_all);
                    plot(ax, [gi-0.30 gi+0.30], [medv medv], 'k-', 'LineWidth', 1.5);
                end

                % capture legend handles once (first tile only)
                if ti == 1
                    if isfield(opts,"colorMap") && ~isempty(opts.colorMap) && isKey(opts.colorMap, gname)
                        if ~any(treatNames == gname)
                            h = scatter(ax, nan, nan, 80, 'filled', ...
                                'MarkerFaceColor', base, 'MarkerEdgeColor','k','LineWidth',0.75);
                            hTreat(end+1) = h; %#ok<AGROW>
                            treatNames(end+1) = gname; %#ok<AGROW>
                        end
                    end
                end
            end

            xlim(ax, [0.5 numel(groups)+0.5]);
            xticks(ax, 1:numel(groups));
            xticklabels(ax, groups);
            xtickangle(ax, 25);

        else
            % ==================== NEW LAYOUT: Time panel > Stage groups > Treatment subgroups ====================

            % choose experiments list from this timepoint (for legend later)
            % but medians are computed per (stage,treat,exp), so per combo.
            allExpsHere = categories(subT.Experiment);
            nEmax = numel(allExpsHere);

            % Build x positions: each Stage is a "block" containing K treatments
            K = numel(treatLevels);
            stageGap = 1.0;  % spacing between stage blocks (tune)
            treatGap = 1.50;  % spacing between treatments inside stage (tune)

            % Precompute x centers for each (stage, treat)
            % Index scheme: block start = (si-1)*(K*treatGap + stageGap) + 1
            xCenter = nan(numel(stageLevels), K);
            for si = 1:numel(stageLevels)
                block0 = (si-1)*(K*treatGap + stageGap) + 1;
                for ki = 1:K
                    xCenter(si,ki) = block0 + (ki-1)*treatGap;
                end
            end

            % Plot each stage block
            for si = 1:numel(stageLevels)
                stName = stageLevels(si);

                % Subset for this stage at this time
                subS = subT(subT.Stage == stName, :);
                if isempty(subS), continue; end

                for ki = 1:K
                    gname = treatLevels(ki);
                    maskG = (subS.PlotGroup == gname);

                    valsG_all = subS.(varName)(maskG);
                    valsG_all = valsG_all(isfinite(valsG_all));

                    % base color by treatment
                    base = [0.5 0.5 0.5];
                    if isfield(opts,"colorMap") && ~isempty(opts.colorMap) && isKey(opts.colorMap, gname)
                        base = opts.colorMap(gname);
                    end

                    % experiments present in THIS stage+treat
                    exps = categories(subS.Experiment(maskG));
                    nE = numel(exps);

                    xc = xCenter(si,ki);

                    % background scatter per experiment (shade)
                    for ei = 1:nE
                        maskE = maskG & (subS.Experiment == exps{ei});
                        vals = subS.(varName)(maskE);
                        vals = vals(isfinite(vals));
                        if isempty(vals), continue; end

                        a = getA(exps{ei});
                        c_bg = base*(1-a) + [1 1 1]*a;

                        xj = xc + 0.22*(rand(size(vals)) - 0.5);
                        s = scatter(ax, xj, vals, 8, 'filled', ...
                            'MarkerFaceColor', c_bg, ...
                            'MarkerEdgeColor', c_bg*0.6, ...
                            'LineWidth', 0.3);
                        s.MarkerFaceAlpha = 0.65;
                    end

                    % per-experiment medians (big dots)
                    for ei = 1:numel(exps)
                        expName = exps{ei};
                        maskE = maskG & (subS.Experiment == expName);
                        valsE = subS.(varName)(maskE);
                        valsE = valsE(isfinite(valsE));
                        if isempty(valsE), continue; end

                        medVal = median(valsE);

                        % get fixed shade factor for this experiment
                        a = getA(expName);
                        c_med = base*(1-a) + [1 1 1]*a;

                        % small horizontal separation between experiments
                        if numel(exps) == 1
                            xMed = xc;
                        else
                            offsets = linspace(-0.10, 0.10, numel(exps));
                            xMed = xc + offsets(ei);
                        end

                        scatter(ax, xMed, medVal, 30, 'filled', ...
                            'MarkerFaceColor', c_med, ...
                            'MarkerEdgeColor', 'k', ...
                            'LineWidth', 1);
                    end

                    % overall median line (pooled across experiments)
                    if ~isempty(valsG_all)
                        medv = median(valsG_all);
                        plot(ax, [xc-0.30 xc+0.30], [medv medv], 'k-', 'LineWidth', 1.5);
                    end

                    % ---- Stats annotation (stars vs Vehicle) ----
                    if opts.byStage && isfield(opts, "statsCache") && opts.annotateStats

                        % Only annotate doses (not Vehicle)
                        thisDose = plotGroupToDose(gname);
                        if thisDose ~= 0 && ~isnan(thisDose)

                            key = matlab.lang.makeValidName(string(varName));
                            if isfield(opts.statsCache, key) && isfield(opts.statsCache.(key), "annot")
                                Ann = opts.statsCache.(key).annot;

                                % Find this time+stage row
                                m = (Ann.Time_min == tval) & (Ann.Stage == string(stName));
                                if any(m)
                                    row = find(m,1,'first');
                                    p10  = Ann.p10(row);
                                    p50  = Ann.p50(row);
                                    p100 = Ann.p100(row);

                                    pHere = NaN;
                                    if thisDose == 10,  pHere = p10;  end
                                    if thisDose == 50,  pHere = p50;  end
                                    if thisDose == 100, pHere = p100; end

                                    stars = pToStars(pHere);
                                    if stars ~= ""

                                        yl = ylim(ax);
                                        yr = diff(yl);
                                        yOff = opts.statsYOffsetFrac * yr;

                                        yTop = max(valsG_all);
                                        if ~isfinite(yTop)
                                            yTop = yl(2) - 2*yOff;
                                        end

                                        % Prevent clipping at top
                                        yStar = min(yTop + yOff, yl(2) - 0.5*yOff);

                                        text(ax, xc, yStar, stars, ...
                                            'HorizontalAlignment','center', ...
                                            'VerticalAlignment','bottom', ...
                                            'FontSize', 12, ...
                                            'FontWeight','bold');
                                    end
                                end
                            end
                        end
                    end

                    % Capture treatment legend handles once (first tile only)
                    if ti == 1 && si == 1
                        if ~any(treatNames == gname)
                            h = scatter(ax, nan, nan, 70, 'filled', ...
                                'MarkerFaceColor', base, 'MarkerEdgeColor','k','LineWidth',0.75);
                            hTreat(end+1) = h; %#ok<AGROW>
                            treatNames(end+1) = gname; %#ok<AGROW>
                        end
                    end
                end

                % draw a faint separator between stage blocks
                if si < numel(stageLevels)
                    xSep = xCenter(si,end) + 0.5*(stageGap);
                    xline(ax, xSep, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
                end
            end

            % X ticks: label each stage at the center of its block
            stageCenters = mean(xCenter, 2);
            xticks(ax, stageCenters);
            if isfield(opts,"stageShortLabels") && opts.stageShortLabels
                stageLabels = stageLevels;
                stageLabels(stageLabels=="Horsetail") = "Ht";
                stageLabels(stageLabels=="Meiosis I") = "M1";
                stageLabels(stageLabels=="Meiosis II") = "M2";
                stageLabels(stageLabels=="Late Meiosis II") = "LM2";
                stageLabels(stageLabels=="Spores") = "S";
                xticklabels(ax, stageLabels);
            else
                xticklabels(ax, stageLevels);
            end
            xtickangle(ax, 0);

            % Add minor tick labels for treatments under each stage block (optional)
            % (kept off by default to avoid clutter)
        end

        % ---- Title per time panel ----
        if tval == 0
            title(ax, "Time 5 min");
        else
            title(ax, sprintf("Time %d min", tval));
        end

        if ti == 1
            ylabel(ax, yLabel);
        else
            ylabel(ax, "");
        end
        grid(ax,'on'); box(ax,'on');
        set(ax,'FontSize',11);

        % Expand y-lims slightly
        yl = ylim(ax); yr = diff(yl);
        if yr > 0
            ylim(ax, [yl(1) - 0.05*yr, yl(2) + 0.08*yr]);
        end

        % Store axis for legend placement
        if isempty(legendAx), legendAx = ax; end
    end

    title(tl, yLabel, 'FontWeight','bold');

    % ---- Combined legend (Treatment + Experiment, fixed shades) ----
    if ~isempty(hTreat)

        % Create a dedicated legend tile on the far right
        axLeg = nexttile(tl, numel(timeLevels)+1);
        axis(axLeg, 'off'); hold(axLeg, 'on');

        % Get global experiment list
        expsAll = categories(T.Experiment);
        nEall = numel(expsAll);

        % Fixed shade factors (must match getA)
        if nEall <= 1
            aAll = opts.expShadeMin;
        else
            aAll = linspace(opts.expShadeMin, opts.expShadeMax, nEall);
        end

        % Dummy handles for experiments (grey shades)
        hExp = gobjects(nEall,1);
        for ei = 1:nEall
            a = aAll(ei);
            c = [0 0 0]*(1-a) + [1 1 1]*a;
            hExp(ei) = scatter(axLeg, nan, nan, 26, 'filled', ...
                'MarkerFaceColor', c, ...
                'MarkerEdgeColor', c*0.6, ...
                'LineWidth', 0.3);
        end

        % Combine treatment + experiment handles
        hAll = [hTreat(:); hExp(:)];
        namesAll = [ ...
            "Treatment: " + treatNames(:); ...
            "Experiment: " + string(expsAll) ...
            ];

        lgd = legend(axLeg, hAll, cellstr(namesAll), "Location","northwest");
        lgd.Box = "off";
        % Place legend inside the legend tile (axLeg), using figure-normalized coords
        lgd.Units = "normalized";
        axLeg.Units = "normalized";

        axp = axLeg.Position;   % [x y w h] in figure-normalized units
        lgp = lgd.Position;     % current legend size [x y w h]

        padX = 0.005 * axp(3);   % 0.5% of legend-tile width
        padY = 0.02 * axp(4);   % 2% of legend-tile height

        % Left-align legend within the legend tile, top-align
        lgp(1) = axp(1) + padX;
        lgp(2) = axp(2) + axp(4) - lgp(4) - padY;

        lgd.Position = lgp;

    end

    if opts.saveFigs
        suffix = stem;

        if opts.byStage
            suffix = suffix + "_stagesWithinPanel";
        end

        outFig = fullfile(opts.outDir, ...
            opts.fileStem + "_" + suffix + "." + opts.figFormat);

        exportgraphics(fig, outFig, "Resolution", 300);
    end
end

function col = defaultColumn(n, cls)
    switch cls
        case {"double","single"}
            col = nan(n,1);
        case {"string"}
            col = strings(n,1);
        case {"cell"}
            col = cell(n,1);
        case {"categorical"}
            col = categorical(repmat(missing, n, 1));
        case {"datetime"}
            col = NaT(n,1);
        case {"logical"}
            col = false(n,1);
        otherwise
            % fallback: make it a string column of missing
            col = strings(n,1);
    end
end

% ==================== PLOT CALLS ====================
dye = opts.dyeName;

% Define metrics depending on mode
if opts.minimalDyeOnly
    metricsRaw = {
        meanVar,   dye + " Mean Intensity",      meanVar
        intDenVar, dye + " Integrated Density",  intDenVar
        };

    metricsNorm = {
        normMeanVar,   dye + " Mean (Normalised)",    normMeanVar
        normIntDenVar, dye + " IntDen (Normalised)",  normIntDenVar
        };
else
    metricsRaw = {
        meanVar,      dye + " Mean Intensity",        meanVar
        intDenVar,    dye + " Integrated Density",    intDenVar
        ratioVarName, "Ratio: GFP / " + dye,          ratioVarName
        };

    metricsNorm = {
        normMeanVar,   dye + " Mean (Normalised)",   normMeanVar
        normIntDenVar, dye + " IntDen (Normalised)", normIntDenVar
        normRatioVar,  "Ratio: GFP / " + dye + " (Normalised)", normRatioVar
        };
end

% ---- Raw plots ----
for mi = 1:size(metricsRaw,1)
    plotMetric(T, metricsRaw{mi,1}, metricsRaw{mi,2}, metricsRaw{mi,3}, opts);
end

% ---- Normalized plots (+ stats + stars) ----
if opts.doStats
    if opts.byStage
        fprintf("\nRunning stats on experiment medians (per Time × Stage)...\n");
    elseif opts.allowStatsWithoutStage
        fprintf("\nRunning stats on experiment medians (per Time; no stage)...\n");
    else
        fprintf("\nStats requested, but byStage=false and allowStatsWithoutStage=false. Skipping stats.\n");
    end
end

for mi = 1:size(metricsNorm,1)

    varName = metricsNorm{mi,1};
    yLabel  = metricsNorm{mi,2};
    stem    = metricsNorm{mi,3};

    % Fresh cache per metric
    statsCache = struct();

    if opts.doStats && (opts.byStage || opts.allowStatsWithoutStage)

        [Ssum, Spair, Sann] = computeStatsByTimeStage(T, string(varName), opts);

        % Save CSVs
        writetable(Ssum,  fullfile(opts.outDir, ...
            opts.fileStem + "_STATS_summary_" + string(varName) + ".csv"));

        writetable(Spair, fullfile(opts.outDir, ...
            opts.fileStem + "_STATS_pairwise_" + string(varName) + ".csv"));

        % Store ONLY this metric
        key = matlab.lang.makeValidName(string(varName));
        statsCache.(key).annot = Sann;

        % Attach to opts
        opts.statsCache = statsCache;
    end

    plotMetric(T, varName, yLabel, stem, opts);
end

% Optional: save parsed + ratio table for MATLAB downstream
writetable(T, fullfile(opts.outDir, opts.fileStem + "_parsed_withRatios_andNorm.csv"));

fprintf("Done. Outputs written to:\n  %s\n", opts.outDir);
end

function [Condition, Dose_mM, Time_min, Series] = parseH2O2Filename(fileNames)
% parseH2O2Filename
%
% Robustly parses condition, dose, time, and series from Fiji File names.
%
% Supports examples such as:
%   H2O2_100mM_60min_series1.tif
%   Fzo1D_H2O2_60min_series1.tif
%   20260121_H2O2_100mM_60min_series1.tif
%   Veh_series1.tif
%   Vehicle_30min_series2.tif

n = numel(fileNames);

Condition = strings(n,1);
Dose_mM   = nan(n,1);
Time_min  = nan(n,1);
Series    = nan(n,1);

for i = 1:n
    fn = string(fileNames(i));
    fn = erase(fn, '"');

    % Work with the file name only, even if a path was stored
    [~, baseName, ext] = fileparts(fn);
    fnClean = lower(baseName + ext);

    % Remove common extensions if present
    fnClean = regexprep(fnClean, "\.(tif|tiff|csv)$", "");

    % -------------------- Condition --------------------
    if contains(fnClean, "h2o2")
        Condition(i) = "H2O2";
    elseif contains(fnClean, "vehicle") || contains(fnClean, "veh")
        Condition(i) = "Veh";
    else
        Condition(i) = "Unparsed";
    end

    % -------------------- Dose --------------------
    if Condition(i) == "Veh"
        Dose_mM(i) = 0;

    elseif Condition(i) == "H2O2"
        % Prefer dose immediately after H2O2, e.g. H2O2_100mM
        dTok = regexp(fnClean, "h2o2[_\-\s]*(\d+)\s*mm", "tokens", "once");

        if isempty(dTok)
            % Fallback: any explicit mM token
            dTok = regexp(fnClean, "(\d+)\s*mm", "tokens", "once");
        end

        if ~isempty(dTok)
            Dose_mM(i) = str2double(dTok{1});
        else
            Dose_mM(i) = NaN;
        end
    end

    % -------------------- Time --------------------
    % Default to 0 if no explicit time is present
    Time_min(i) = 0;

    tTok = regexp(fnClean, "(\d+)\s*mins?", "tokens", "once");

    if ~isempty(tTok)
        Time_min(i) = str2double(tTok{1});
    else
        hTok = regexp(fnClean, "(\d+)\s*h\b", "tokens", "once");

        if ~isempty(hTok)
            Time_min(i) = 60 * str2double(hTok{1});
        end
    end

    % -------------------- Series --------------------
    sTok = regexp(fnClean, "series\s*(\d+)", "tokens", "once");

    if ~isempty(sTok)
        Series(i) = str2double(sTok{1});
    end
end

end

function timeLevels = orderTimeLevels(T, opts)
% Preferred time order for plots, with unexpected times appended.
% If opts.timePoints is provided, that order is respected.

detectedTimes = unique(T.Time_min(~isnan(T.Time_min)));
detectedTimes = detectedTimes(:)';

if ~isempty(opts.timePoints)
    prefTimes = opts.timePoints(:)';
else
    prefTimes = [0 20 30 60];
end

ordered = prefTimes(ismember(prefTimes, detectedTimes));
extras  = setdiff(detectedTimes, prefTimes);

timeLevels = [ordered, extras(:)'];

end

function plotSingleGrouped(T, varName, yLabel, stem, opts)
% plotSingleGrouped
%
% Single-axis grouped layout.
%
% No panel format. Data are arranged as:
%
%   Vehicle 5 mins | Vehicle 60 mins |
%   H2O2 10 mM 5 mins | H2O2 10 mM 60 mins |
%   H2O2 50 mM 5 mins | H2O2 50 mM 60 mins |
%   ...
%
% If byStage=true, this structure is repeated for each stage.

varName = string(varName);

if ~ismember(varName, string(T.Properties.VariableNames))
    warning("Column %s not found. Skipping plot.", varName);
    return;
end

valsAll = T.(varName);

if ~isnumeric(valsAll)
    warning("Column %s is not numeric. Skipping plot.", varName);
    return;
end

% Keep only finite values for this plot
plotKeep = isfinite(valsAll);

% Optional plot-only outlier removal.
% This does NOT affect the parsed table, normalization, counts, or stats.
if opts.removePlotOutliers
    finiteVals = valsAll(plotKeep);

    [lo, hi] = getPlotOutlierLimits(finiteVals, opts);

    outlierKeep = valsAll >= lo & valsAll <= hi;
    plotKeep = plotKeep & outlierKeep;

    fprintf("Plot outlier filter for %s: kept %d of %d finite values; limits = [%g, %g].\n", ...
        varName, nnz(plotKeep), nnz(isfinite(valsAll)), lo, hi);
end

T = T(plotKeep, :);

if isempty(T)
    warning("No finite values remain for %s after plot filtering. Skipping plot.", varName);
    return;
end

% Remove unused categories
if iscategorical(T.PlotGroup)
    T.PlotGroup = removecats(T.PlotGroup);
end

if opts.byStage && iscategorical(T.Stage)
    T.Stage = removecats(T.Stage);
end

timeLevels = orderTimeLevels(T, opts);
groupLevels = string(categories(T.PlotGroup));

% Remove truly unparsed groups from plotting unless they are the only group
if numel(groupLevels) > 1
    groupLevels = groupLevels(groupLevels ~= "Unparsed");
end

% Stage levels
if opts.byStage
    stagePref = ["Horsetail","Meiosis I","Meiosis II","Late Meiosis II","Spores"];

    stageCats = string(categories(T.Stage));
    stageCats = stageCats(:);

    prefPresent = stagePref(ismember(stagePref, stageCats));
    prefPresent = prefPresent(:);

    extras = setdiff(stageCats, stagePref);
    extras = string(extras);
    extras = extras(:);

    stageLevels = [prefPresent; extras];
    stageLevels = stageLevels(stageLevels ~= "DELETE");
else
    stageLevels = "ALL";
end

% Experiment shade factors
expsAll = categories(T.Experiment);
nEall = numel(expsAll);

if nEall <= 1
    aAll = opts.expShadeMin;
else
    aAll = linspace(opts.expShadeMin, opts.expShadeMax, nEall);
end

getA = @(expCat) aAll(find(strcmp(expsAll, expCat), 1, 'first'));

% Layout constants
nStages = numel(stageLevels);
nGroups = numel(groupLevels);
nTimes  = numel(timeLevels);

timeGap  = 0.50;
groupGap = 0.5;
stageGap = 0.5;

xTicks = [];
xTickLabels = strings(0);

xMap = nan(nStages, nGroups, nTimes);

xCursor = 1;

for si = 1:nStages
    for gi = 1:nGroups
        for ti = 1:nTimes
            xMap(si, gi, ti) = xCursor;

            switch opts.xTickLabelMode

                case "timeonly"
                    xTickLabels(end+1) = string(formatTimeLabel(timeLevels(ti), opts)); %#ok<AGROW>

                case "full"
                    if opts.byStage
                        xTickLabels(end+1) = sprintf("%s %s %s", ...
                            string(shortStageLabel(stageLevels(si))), ...
                            string(groupLevels(gi)), ...
                            string(formatTimeLabel(timeLevels(ti), opts))); %#ok<AGROW>
                    else
                        xTickLabels(end+1) = sprintf("%s %s", ...
                            string(groupLevels(gi)), ...
                            string(formatTimeLabel(timeLevels(ti), opts))); %#ok<AGROW>
                    end
            end

            xTicks(end+1) = xCursor; %#ok<AGROW>
            xCursor = xCursor + timeGap;
        end

        % Extra gap between treatment groups
        xCursor = xCursor + groupGap;
    end

    % Extra gap between stage blocks
    xCursor = xCursor + stageGap;
end

nX = numel(xTicks);

fig = figure("Color","w", ...
    "Units","centimeters", ...
    "Position",[2 2 opts.figSizeCm(1) opts.figSizeCm(2)], ...
    "PaperUnits","centimeters", ...
    "PaperSize", opts.figSizeCm, ...
    "PaperPosition", [0 0 opts.figSizeCm(1) opts.figSizeCm(2)], ...
    "Name", string(varName));

ax = axes(fig);
hold(ax, "on");

ax.Position = [0.28 0.18 0.80 0.65];

% Plotting aesthetics
dotSize  = opts.singlePointSize;
medSize  = opts.experimentMedianSize;
medianLW = opts.pooledMedianLineWidth;
medianHalfWidth = opts.pooledMedianHalfWidth;

% Main plotting loop
for si = 1:nStages

    if opts.byStage
        stName = stageLevels(si);
        subStage = T(T.Stage == stName, :);
    else
        stName = "ALL";
        subStage = T;
    end

    if isempty(subStage)
        continue;
    end

    for gi = 1:nGroups
        gname = groupLevels(gi);

        for ti = 1:nTimes
            tval = timeLevels(ti);
            xc = xMap(si, gi, ti);

            maskG = (string(subStage.PlotGroup) == gname) & ...
                    (subStage.Time_min == tval);

            valsG_all = subStage.(varName)(maskG);
            valsG_all = valsG_all(isfinite(valsG_all));

            if isempty(valsG_all)
                continue;
            end

            % Base color by treatment group
            base = [0.5 0.5 0.5];

            if isfield(opts, "colorMap") && ~isempty(opts.colorMap) && isKey(opts.colorMap, gname)
                base = opts.colorMap(gname);
            end

            % Raw points by experiment
            expsHere = categories(subStage.Experiment(maskG));

            for ei = 1:numel(expsHere)
                expName = expsHere{ei};

                maskE = maskG & (subStage.Experiment == expName);
                valsE = subStage.(varName)(maskE);
                valsE = valsE(isfinite(valsE));

                if isempty(valsE)
                    continue;
                end

                a = getA(expName);
                c_bg = base*(1-a) + [1 1 1]*a;

                xj = xc + 0.22*(rand(size(valsE)) - 0.5);

                s = scatter(ax, xj, valsE, dotSize, "filled", ...
                    "MarkerFaceColor", c_bg, ...
                    "MarkerEdgeColor", c_bg*0.6, ...
                    "LineWidth", 0.3);

                s.MarkerFaceAlpha = 0.65;
            end

            % Per-experiment medians
            for ei = 1:numel(expsHere)
                expName = expsHere{ei};

                maskE = maskG & (subStage.Experiment == expName);
                valsE = subStage.(varName)(maskE);
                valsE = valsE(isfinite(valsE));

                if isempty(valsE)
                    continue;
                end

                medVal = median(valsE, "omitnan");

                a = getA(expName);
                c_med = base*(1-a) + [1 1 1]*a;

                if numel(expsHere) == 1
                    xMed = xc;
                else
                    offsets = linspace(-0.12, 0.12, numel(expsHere));
                    xMed = xc + offsets(ei);
                end

                scatter(ax, xMed, medVal, medSize, "filled", ...
                    "MarkerFaceColor", c_med, ...
                    "MarkerEdgeColor", "k", ...
                    "LineWidth", 1);
            end

            % Pooled median line
            medPooled = median(valsG_all, "omitnan");

            plot(ax, [xc - medianHalfWidth, xc + medianHalfWidth], ...
                [medPooled, medPooled], ...
                "k-", ...
                "LineWidth", medianLW);

            % Optional stats annotation
            if isfield(opts, "statsCache") && opts.annotateStats

                thisDose = plotGroupToDose(gname);

                if thisDose ~= 0 && ~isnan(thisDose)

                    key = matlab.lang.makeValidName(string(varName));

                    if isfield(opts.statsCache, key) && isfield(opts.statsCache.(key), "annot")
                        Ann = opts.statsCache.(key).annot;

                        m = (Ann.Time_min == tval) & (Ann.Stage == string(stName));

                        if any(m)
                            row = find(m, 1, "first");

                            pHere = NaN;

                            if thisDose == 10
                                pHere = Ann.p10(row);
                            elseif thisDose == 50
                                pHere = Ann.p50(row);
                            elseif thisDose == 100
                                pHere = Ann.p100(row);
                            end

                            stars = pToStars(pHere);

                            if stars ~= ""
                                yl = ylim(ax);
                                yr = diff(yl);
                                yOff = opts.statsYOffsetFrac * yr;

                                yTop = max(valsG_all);

                                if ~isfinite(yTop)
                                    yTop = yl(2) - 2*yOff;
                                end

                                yStar = min(yTop + yOff, yl(2) - 0.5*yOff);

                                text(ax, xc, yStar, stars, ...
                                    "HorizontalAlignment","center", ...
                                    "VerticalAlignment","bottom", ...
                                    "FontSize", 12, ...
                                    "FontWeight","bold");
                            end
                        end
                    end
                end
            end
        end
    end

    % Stage separator
    if opts.byStage && si < nStages
        xSep = max(xMap(si,:,:), [], "all") + 0.5*stageGap;
        xline(ax, xSep, "-", ...
            "Color", [0.85 0.85 0.85], ...
            "LineWidth", 1);
    end
end

% Axis formatting
xlim(ax, [min(xTicks)-0.7, max(xTicks)+0.7]);

if opts.hideXAxisTicks
    xticks(ax, []);
    xticklabels(ax, {});
else
    xticks(ax, xTicks);
    xticklabels(ax, cellstr(xTickLabels));
    xtickangle(ax, opts.xTickAngle);
end

ylabel(ax, yLabel, ...
    "Interpreter", "none", ...
    "FontSize", opts.yLabelFontSize, ...
    "FontWeight", opts.yLabelFontWeight);

% Remove top x-axis and right y-axis
box(ax, "off");

% Remove background grid
if opts.showGrid
    grid(ax, "on");
else
    grid(ax, "off");
end

% Remove figure/axis title by default
if opts.showFigureTitle
    title(ax, yLabel, "Interpreter", "none", "FontWeight", "bold");
else
    title(ax, "");
end

% No legend by default
if opts.showLegend
    legend(ax, "show", "Location", "best");
end

set(ax, ...
    "FontSize", opts.xTickFontSize, ...
    "TickDir", "out", ...
    "LineWidth", opts.axisLineWidth);

ax.YAxis.FontSize = opts.yTickFontSize;
ax.XAxis.FontSize = opts.xTickFontSize;

% Make only x-tick labels bold
ax.XAxis.FontWeight = "bold";
ax.YAxis.FontWeight = "normal";

% Re-apply y-axis label formatting after y-axis tick formatting
ax.YLabel.FontSize = opts.yLabelFontSize;
ax.YLabel.FontWeight = opts.yLabelFontWeight;
ax.YLabel.Interpreter = "none";

% Slight y-padding
yl = ylim(ax);
yr = diff(yl);

if yr > 0
    ylim(ax, [yl(1) - 0.05*yr, yl(2) + 0.08*yr]);
end

% Save
if opts.saveFigs
    suffix = stem;

    if opts.byStage
        suffix = suffix + "_singleGrouped_byStage";
    else
        suffix = suffix + "_singleGrouped";
    end

    outFig = fullfile(opts.outDir, ...
        opts.fileStem + "_" + suffix + "." + opts.figFormat);

    exportgraphics(fig, outFig, "Resolution", 300);
end

end

function lab = formatTimeLabel(tval, opts)
% Converts numeric Time_min values into plot labels.
%
% opts.timeLabelMode:
%   "short"   -> 5 min, 60 min
%   "compact" -> 5', 60'

switch lower(opts.timeLabelMode)

    case "compact"
        if tval == 0
            % Extract the number from timeZeroLabel, e.g. "5 min" -> "5'"
            tok = regexp(opts.timeZeroLabel, "\d+", "match", "once");

            if isempty(tok)
                lab = opts.timeZeroLabel;
            else
                lab = tok + "'";
            end
        else
            lab = sprintf("%g'", tval);
        end

    case "short"
        if tval == 0
            lab = opts.timeZeroLabel;
        else
            lab = sprintf("%g min", tval);
        end
end

end

function lab = shortStageLabel(stageName)
% Short stage labels for compact grouped x-axis labels.

stageName = string(stageName);

switch stageName
    case "Horsetail"
        lab = "Ht";
    case "Meiosis I"
        lab = "M1";
    case "Meiosis II"
        lab = "M2";
    case "Late Meiosis II"
        lab = "LM2";
    case "Spores"
        lab = "S";
    otherwise
        lab = stageName;
end

end

function [lo, hi] = getPlotOutlierLimits(vals, opts)
% getPlotOutlierLimits
%
% Returns lower and upper limits for plot-only outlier removal.

vals = vals(isfinite(vals));

if isempty(vals)
    lo = -Inf;
    hi = Inf;
    return;
end

switch opts.plotOutlierMethod

    case "percentile"
        lo = prctile(vals, opts.plotOutlierPercentiles(1));
        hi = prctile(vals, opts.plotOutlierPercentiles(2));

    case "iqr"
        q1 = prctile(vals, 25);
        q3 = prctile(vals, 75);
        iq = q3 - q1;

        lo = q1 - opts.plotOutlierIQRMultiplier * iq;
        hi = q3 + opts.plotOutlierIQRMultiplier * iq;

    case "absolute"
        lo = opts.plotOutlierLimits(1);
        hi = opts.plotOutlierLimits(2);

    otherwise
        lo = -Inf;
        hi = Inf;
end

end

%-------------- STATS HELPERS --------------%
function stars = pToStars(p)
    if isnan(p) || isempty(p)
        stars = "";
    elseif p < 1e-4
        stars = "****";
    elseif p < 1e-3
        stars = "***";
    elseif p < 1e-2
        stars = "**";
    elseif p < 0.05
        stars = "*";
    else
        stars = "";
    end
end

function d = plotGroupToDose(pg)
    pg = string(pg);
    if pg == "Vehicle"
        d = 0;
        return;
    end
    tok = regexp(pg, "(\d+)\s*mM", "tokens", "once");
    if isempty(tok)
        d = NaN;
    else
        d = str2double(tok{1});
    end
end

function [summaryTbl, pairwiseTbl, annotTbl] = computeStatsByTimeStage(T, varName, opts)
% Computes experiment-median stats per Time×Stage for the variable varName.
% Returns:
%   summaryTbl: one row per Time×Stage
%   pairwiseTbl: Dunn-Sidak pairwise comparisons per Time×Stage
%   annotTbl: p-values Vehicle vs dose (for plot asterisks)

    % Build per-experiment medians (unit = Experiment)
    grpVars = ["Time_min","Stage","PlotGroup","Experiment"];
    G = groupsummary(T, grpVars, "median", varName);
    medCol = "median_" + varName;

    times  = unique(G.Time_min);
    stages = unique(G.Stage);

    summaryRows = table();
    pairRows    = table();
    annotRows   = table();

    for ti = 1:numel(times)
        tval = times(ti);

        for si = 1:numel(stages)
            st = stages(si);

            sub = G(G.Time_min==tval & G.Stage==st, :);
            if isempty(sub)
                continue;
            end

            % Values + groups
            y = sub.(medCol);
            g = sub.PlotGroup;

            % Drop nonfinite
            ok = isfinite(y) & ~isundefined(g);
            y = y(ok); g = g(ok);
            if numel(y) < 2 || numel(unique(g)) < 2
                continue;
            end

            % KW
            nGroups = numel(unique(g));

            if nGroups == 2

                gLevels = string(categories(categorical(g)));

                y1 = y(string(g) == gLevels(1));
                y2 = y(string(g) == gLevels(2));

                if numel(y1) >= 2 && numel(y2) >= 2
                    pKW = ranksum(y1, y2);
                else
                    pKW = NaN;
                end

                pw = table(tval, string(st), gLevels(1), gLevels(2), pKW, ...
                    'VariableNames', {'Time_min','Stage','Group1','Group2','p_adj'});

            else
                % original kruskalwallis + Dunn block

                pKW = NaN;
                stats = [];
                try
                    [pKW,~,stats] = kruskalwallis(y, g, 'off');
                catch
                    pKW = NaN;
                    stats = [];
                end

                % Dunn-Sidak pairwise (if possible)
                pw = table();
                if ~isempty(stats) && isstruct(stats)
                    try
                        c = multcompare(stats, 'CType','dunn-sidak', 'Display','off');
                        % c columns: [i j lower diff upper p]
                        cats = string(stats.gnames);
                        pw = table( ...
                            repmat(tval,size(c,1),1), repmat(string(st),size(c,1),1), ...
                            cats(c(:,1)), cats(c(:,2)), ...
                            c(:,6), ...
                            'VariableNames', {'Time_min','Stage','Group1','Group2','p_adj'} ...
                            );
                    catch
                        pw = table();
                    end
                end
            end
            % Spearman trend (dose vs value)
            dose = arrayfun(@plotGroupToDose, g);
            okd = isfinite(dose) & isfinite(y);
            rho = NaN; pTrend = NaN;
            if nnz(okd) >= 3 && numel(unique(dose(okd))) >= 2
                [rho,pTrend] = corr(dose(okd), y(okd), 'Type','Spearman', 'Rows','complete');
            end

            % Counts per treatment (how many experiments contributed)
            catsHere = categories(categorical(g));
            vehName = opts.vehicleLabel;

            nVeh  = sum(string(g) == vehName);
            n10   = sum(contains(string(g), "10 mM"));
            n50   = sum(contains(string(g), "50 mM"));
            n100  = sum(contains(string(g), "100 mM"));

            % Vehicle vs each dose p
            pV10  = NaN;
            pV50  = NaN;
            pV100 = NaN;

            if ~isempty(pw)
                pV10  = pickPairP(pw, vehName, "H2O2 10 mM");
                pV50  = pickPairP(pw, vehName, "H2O2 50 mM");
                pV100 = pickPairP(pw, vehName, "H2O2 100 mM");
            else
                % fallback, unadjusted ranksum if Dunn/multcompare not available
                pV10  = ranksumSafe(y, g, vehName, "H2O2 10 mM");
                pV50  = ranksumSafe(y, g, vehName, "H2O2 50 mM");
                pV100 = ranksumSafe(y, g, vehName, "H2O2 100 mM");
            end

            % Summary row
            summaryRows = [summaryRows; ...
                table(string(varName), tval, string(st), pKW, rho, pTrend, ...
                      nVeh, n10, n50, n100, pV10, pV50, pV100, ...
                      'VariableNames', ...
                      {'Metric','Time_min','Stage','KW_p','Spearman_rho','Spearman_p', ...
                       'nExp_Veh','nExp_10mM','nExp_50mM','nExp_100mM', ...
                       'p_Veh_vs_10','p_Veh_vs_50','p_Veh_vs_100'})]; %#ok<AGROW>

            % Pairwise rows
            if ~isempty(pw)
                pw.Metric = repmat(string(varName), height(pw), 1);
                pw = movevars(pw, "Metric", "Before", 1);
                pairRows = [pairRows; pw]; %#ok<AGROW>
            end

            % Annotation rows (for plotting stars)
            annotRows = [annotRows; ...
                table(string(varName), tval, string(st), ...
                      pV10, pV50, pV100, ...
                      'VariableNames', {'Metric','Time_min','Stage','p10','p50','p100'})]; %#ok<AGROW>
        end
    end

    summaryTbl  = summaryRows;
    pairwiseTbl = pairRows;
    annotTbl    = annotRows;

    if isempty(summaryTbl)
        summaryTbl = table( ...
            strings(0,1), zeros(0,1), strings(0,1), ...
            zeros(0,1), zeros(0,1), zeros(0,1), ...
            zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
            zeros(0,1), zeros(0,1), zeros(0,1), ...
            'VariableNames', ...
            {'Metric','Time_min','Stage','KW_p','Spearman_rho','Spearman_p', ...
            'nExp_Veh','nExp_10mM','nExp_50mM','nExp_100mM', ...
            'p_Veh_vs_10','p_Veh_vs_50','p_Veh_vs_100'});
    end

    if isempty(pairwiseTbl)
        pairwiseTbl = table( ...
            strings(0,1), zeros(0,1), strings(0,1), ...
            strings(0,1), strings(0,1), zeros(0,1), ...
            'VariableNames', ...
            {'Metric','Time_min','Stage','Group1','Group2','p_adj'});
    end

    if isempty(annotTbl)
        annotTbl = table( ...
            strings(0,1), zeros(0,1), strings(0,1), ...
            zeros(0,1), zeros(0,1), zeros(0,1), ...
            'VariableNames', ...
            {'Metric','Time_min','Stage','p10','p50','p100'});
    end
end

function p = pickPairP(pw, a, b)
    a = string(a); b = string(b);
    m = (pw.Group1==a & pw.Group2==b) | (pw.Group1==b & pw.Group2==a);
    if any(m)
        p = pw.p_adj(find(m,1,'first'));
    else
        p = NaN;
    end
end

function p = ranksumSafe(y, g, a, b)
    a = string(a); b = string(b);
    ya = y(g==a); yb = y(g==b);
    if numel(ya) < 2 || numel(yb) < 2
        p = NaN;
        return;
    end
    p = ranksum(ya, yb);
end

% To RUN
%
%analyze_MTxRos_fromScratch( ...
%"csv.file", ...
%outDir="outputDir", ...
%saveFigs=true, ...
%fileStem="H2O2Ac_MTxRos" ...
%byStage=false(true if analysis by stage)
%);

% To RUN Legacy Mode (Simple csvs with just Mean and IntDen
%analyse_H2O2perturbation( ...
%    "V:\path\to\your\legacy_dataset.csv", ...
%    outDir="V:\path\to\output_folder", ...
%    dyeName="MTxRos", ...
%    fileStem="H2O2Ac_MTxRos_legacy", ...
%    byStage=true, ...
%    minimalDyeOnly=true, ...
%    useProvidedNorm=true, ...
%    forceSingleTime=true, ...
%    singleTime_min=0, ...
%    forcePlotGroups=true, ...
%    treatedKey="H2O2", ...
%    treatedLabel="H2O2 10 mM", ...
%    doStats=true, ...
%    saveFigs=true ...
%);
