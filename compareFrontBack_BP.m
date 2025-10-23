function T = compareFrontBack_BP(distProfileGFP, distProfileRFP, outputCSV, outputDir, stageName)
%COMPAREFRONTBACK_BOXLOT Compare front vs back mitochondrial intensity for GFP and RFP
%   Uses cell-level intensity profiles and plots per-cell scatter + mean ± SD.
%   Skips stats if fewer than 3 valid data points. Saves figure and summary CSV.

    try
        %% Ensure input profiles are valid
        if isempty(distProfileGFP) || isempty(distProfileRFP)
            warning('Empty input data for %s. Skipping.', stageName);
            T = table();
            return;
        end

        %% Prepare output directory
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        %% === EXTRACT FRONT AND BACK INTENSITIES PER CELL ===
        fG = nanmean(distProfileGFP(:, 1:50), 2);
        bG = nanmean(distProfileGFP(:, 51:100), 2);
        fR = nanmean(distProfileRFP(:, 1:50), 2);
        bR = nanmean(distProfileRFP(:, 51:100), 2);

        valid = ~isnan(fG) & ~isnan(bG) & ~isnan(fR) & ~isnan(bR);
        fG = fG(valid); bG = bG(valid);
        fR = fR(valid); bR = bR(valid);

        %% === COMPILE DATA FOR PLOTTING ===
        frontGFP = fG(:)';
        backGFP  = bG(:)';
        frontRFP = fR(:)';
        backRFP  = bR(:)';

        %% === SCATTER PLOT WITH JITTER AND MEAN ± SD ===
        f = figure('Visible','off');
        hold on;

        jitterScale = 0.1;
        scatter(1 + jitterScale * randn(size(frontGFP)), frontGFP, 50, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(2 + jitterScale * randn(size(backGFP)), backGFP, 50, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(3 + jitterScale * randn(size(frontRFP)), frontRFP, 50, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(4 + jitterScale * randn(size(backRFP)), backRFP, 50, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.5);

        % Compute and overlay mean ± SD
        dataGroups = {frontGFP, backGFP, frontRFP, backRFP};
        for i = 1:4
            y = dataGroups{i};
            m = mean(y, 'omitnan');
            s = std(y, 'omitnan');
            errorbar(i, m, s, 'k', 'LineWidth', 1.5, 'CapSize', 10);
            plot(i, m, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        end

        % Axis formatting
        xlim([0.5, 4.5]);
        ylim([0, 1]);
        xticks(1:4);
        xticklabels({'Front GFP', 'Back GFP', 'Front RFP', 'Back RFP'});
        ylabel('Normalized Intensity');
        title(sprintf('Front vs Back Intensity - %s', stageName), 'Interpreter','none');
        grid on; box on;
        set(gca, 'FontSize', 12);

        %% === STATISTICS (skip if n < 3) ===
        if numel(frontGFP) >= 3 && numel(backGFP) >= 3
            [~, pGFP] = ttest(frontGFP, backGFP);
        else
            pGFP = NaN;
            disp(['Skipping GFP stats for ', stageName, ' (n < 3)']);
        end

        if numel(frontRFP) >= 3 && numel(backRFP) >= 3
            [~, pRFP] = ttest(frontRFP, backRFP);
        else
            pRFP = NaN;
            disp(['Skipping RFP stats for ', stageName, ' (n < 3)']);
        end

        % Annotate significance
        yMax = max([frontGFP, backGFP, frontRFP, backRFP], [], 'omitnan');
        yOffsetGFP = 0.05 * yMax;
        yOffsetRFP = 0.10 * yMax;

        if ~isnan(pGFP) && pGFP < 0.05
            plot([1 2], [yMax - yOffsetGFP, yMax - yOffsetGFP], 'k', 'LineWidth', 1.5);
            plot([1 1], [yMax - yOffsetGFP - 0.015*yMax, yMax - yOffsetGFP], 'k');
            plot([2 2], [yMax - yOffsetGFP - 0.015*yMax, yMax - yOffsetGFP], 'k');
            text(1.5, yMax - yOffsetGFP + 0.02*yMax, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
        end
        if ~isnan(pRFP) && pRFP < 0.05
            plot([3 4], [yMax - yOffsetRFP, yMax - yOffsetRFP], 'k', 'LineWidth', 1.5);
            plot([3 3], [yMax - yOffsetRFP - 0.015*yMax, yMax - yOffsetRFP], 'k');
            plot([4 4], [yMax - yOffsetRFP - 0.015*yMax, yMax - yOffsetRFP], 'k');
            text(3.5, yMax - yOffsetRFP + 0.02*yMax, '*', 'HorizontalAlignment', 'center', 'FontSize', 20, 'FontWeight', 'bold');
        end

        %% === SAVE FIGURE ===
        saveas(f, fullfile(outputDir, sprintf('%s_FrontBackScatter.png', stageName)));
        close(f);

        %% === CREATE OUTPUT TABLE ===
        numCells = length(frontGFP);
        T = table((1:numCells)', frontGFP', backGFP', (frontGFP ./ backGFP)', ...
                              frontRFP', backRFP', (frontRFP ./ backRFP)', ...
            'VariableNames', {'Cell', 'Front_GFP', 'Back_GFP', 'Ratio_GFP', ...
                                     'Front_RFP', 'Back_RFP', 'Ratio_RFP'});
        writetable(T, outputCSV);

        fprintf('\n%s: GFP p=%.4f | RFP p=%.4f | n=%d cells\n', stageName, pGFP, pRFP, numCells);

    catch ME
        warning('Error in %s: %s', stageName, ME.message);
        T = table();
    end
end



%Run the following two lines to apply mitoIntensityFS and compare front vs back
%[distProfileGFP, distProfileRFP] = mitoIntensityFS(pathGFP, pathRFP, outputDir, stageName); 
%compareFrontBack_BP(distProfileGFP, distProfileRFP, outputCSV, outputDir, stageName)