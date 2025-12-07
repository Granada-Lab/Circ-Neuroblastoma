function cwt5_plot_ridge_readout_distribution(reporter_colors, reporter_colors_dark)
% Carolin Ector, 02.11.2023
% Function plots distribution of continuous amplitudes and periods per cellline and per subtype

% Define parameters
all_celllines = {'IMR5', 'SKNAS', 'CHP212', 'SY5Y', 'GIMEN', 'Kelly', ...
                 'CLBGA', 'Lan5', 'NGP', 'SKNBE', 'SKNSH', 'SKNBE2'};
values = {'period', 'amplitude'};
xaxis_labels = {'Period (hours)', 'Amplitude (a.u.)'};
reporter_legend = {'Bmal1', 'Bmal1 fit', 'Per2', 'Per2 fit'};

disp('cwt5_plot_ridge_readout_distribution.m is executed')

% Load data
load('Results/sorted_periods_amplitudes.mat')

for v = 1:numel(values) % Loop through period and amplitude
    
    fprintf('\n=== Processing %s ===\n', values{v});
    
    % Load sorted continuous periods or amplitudes
    sorted_data = sorted_periods_amplitudes{v};
    
    % Find valid cell lines with data
    [valid_columns, valid_celllines] = find_valid_celllines(sorted_data, all_celllines);
    
    if isempty(valid_columns)
        fprintf('Warning: No valid data found for %s. Skipping.\n', values{v});
        continue;
    end
    
    % Filter data to only valid columns
    sorted_data_filtered = sorted_data(valid_columns);
    celllines = valid_celllines;
    
    fprintf('Found %d valid cell lines: %s\n', length(celllines), strjoin(celllines, ', '));
    
    % Create main figure
    fig = figure;
    fig.Position = [1, 1, 1440, 821];
    
    % Calculate subplot layout
    div = 4;
    n = ceil(numel(celllines) / div);
    
    % Initialize storage variables
    clear meandata
    meandata = cell(2, numel(celllines));
    
    % Process each cell line
    for c = 1:numel(celllines)
        fprintf('Processing cell line: %s\n', celllines{c});
        
        % Load cell line data
        cellline_data = sorted_data_filtered{c};
        legend_entries = {}; % Store legend entries for this subplot
        
        % Determine which reporters have data
        [reporter_start, reporter_end] = find_valid_reporters(cellline_data);
        
        if reporter_start > reporter_end
            fprintf('  Warning: No valid reporters for %s. Skipping.\n', celllines{c});
            continue;
        end
        
        % Process each reporter
        for a = reporter_start:reporter_end
            
            % Extract and validate reporter data
            [data_matrix, num_timeseries] = extract_reporter_data(cellline_data, a);
            
            if isempty(data_matrix)
                fprintf('    Warning: No data for reporter %d in %s. Skipping.\n', a, celllines{c});
                continue;
            end
            
            % Calculate mean values across replicates
            mean_vals = calculate_mean_values(data_matrix);
            
            if isempty(mean_vals)
                fprintf('    Warning: No valid mean values for reporter %d in %s. Skipping plot.\n', a, celllines{c});
                continue;
            end
            
            % Store for overlay plots
            meandata{a, c} = mean_vals;
            
            % Add legend entry with time-series count
            legend_entries{end+1} = sprintf('%s (%d)', reporter_legend{a}, num_timeseries);
            
            % Create subplot
            subplot(div, n, c);
            
            % Plot based on value type
            if v == 1 % Periods - use polar histogram
                plot_period_polar_histogram(mean_vals, reporter_colors{a});
            else % Amplitudes - use regular histogram
                plot_amplitude_histogram(mean_vals, reporter_colors{a}, reporter_colors_dark{a});
            end
            
        end % End reporter loop
        
        % Finalize subplot
        finalize_subplot(celllines{c}, v, legend_entries);
        
    end % End cell line loop
    
    % Add overall figure labels and save
    add_figure_labels(fig, xaxis_labels{v});
    figure_name = create_figure_name('histogram', values{v}, 'avgpercellline');
    if v == 1
        figure_name = create_figure_name('polarhistogram', values{v}, 'avgpercellline');
    end
    saveas(fig, figure_name, 'svg');
    
    % Create overlay plot
    create_overlay_plot(meandata, celllines, values{v}, xaxis_labels{v}, v);
    
end % End values loop

disp('cwt5_plot_ridge_readout_distribution.m is completed')

end % End main function

%% Helper Functions

function [valid_columns, valid_celllines] = find_valid_celllines(sorted_data, all_celllines)
    % Find columns with valid data
    valid_columns = [];
    valid_celllines = {};
    
    for c = 1:size(sorted_data, 2)
        has_data = false;
        
        if ~isempty(sorted_data{c})
            cell_data = sorted_data{c};
            
            % Check if any cell contains non-NaN data
            for i = 1:numel(cell_data)
                if ~isempty(cell_data{i}) && any(~isnan(cell_data{i}))
                    has_data = true;
                    break;
                end
            end
        end
        
        if has_data
            valid_columns(end+1) = c;
            if c <= length(all_celllines)
                valid_celllines{end+1} = all_celllines{c};
            else
                valid_celllines{end+1} = sprintf('CellLine_%d', c);
            end
        end
    end
end

function [reporter_start, reporter_end] = find_valid_reporters(cellline_data)
    % Determine which reporters have data
    reporter_start = 1;
    reporter_end = 2;
    
    num_rows = size(cellline_data, 1);
    for row_idx = 1:num_rows
        row_empty_check = cellfun(@isempty, cellline_data(row_idx, :));
        
        if all(row_empty_check)
            fprintf('  Reporter row %d is completely empty.\n', row_idx);
            if row_idx == 1
                reporter_start = 2;
            elseif row_idx == 2
                reporter_end = 1;
            end
        end
    end
end

function [data_matrix, num_timeseries] = extract_reporter_data(cellline_data, reporter_idx)
    % Extract and organize reporter data
    reporter_data = cellline_data(reporter_idx, :);
    
    valid_replicates = {};
    
    % Collect valid replicates
    for rep = 1:length(reporter_data)
        if ~isempty(reporter_data{rep}) && any(~isnan(reporter_data{rep}))
            valid_replicates{end+1} = reporter_data{rep}(:); % Store as column vector
        end
    end
    
    if isempty(valid_replicates)
        data_matrix = [];
        num_timeseries = 0;
        return;
    end
    
    % Create data matrix - each column is a replicate
    max_length = max(cellfun(@length, valid_replicates));
    data_matrix = nan(max_length, length(valid_replicates));
    
    for rep = 1:length(valid_replicates)
        rep_data = valid_replicates{rep};
        data_matrix(1:length(rep_data), rep) = rep_data;
    end
    
    num_timeseries = length(valid_replicates);
end

function mean_vals = calculate_mean_values(data_matrix)
    % Calculate mean values across replicates
    mean_vals = [];
    
    for f = 1:size(data_matrix, 1)
        mean_val = mean(data_matrix(f, :), 'omitnan');
        if ~isnan(mean_val)
            mean_vals(end+1) = mean_val;
        end
    end
end

function plot_period_polar_histogram(periods, color)
    % Plot periods in polar histogram
    valid_periods = periods(~isnan(periods) & periods >= 20 & periods <= 40);
    
    if length(valid_periods) < 2
        text(0.5, 0.5, sprintf('No valid periods\n(%d points)', length(valid_periods)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'red');
        return;
    end
    
    % Convert periods to angles (22-36 hours to 0-π radians)
    angles = (valid_periods - 22) / (36 - 22) * pi;
    angles = angles(angles >= 0 & angles <= pi);
    
    if isempty(angles)
        return;
    end
    
    % Create histogram
    num_bins = min(10, max(3, length(angles)));
    p = polarhistogram(angles, num_bins, 'Normalization', 'pdf');
    hold on;
    
    p.FaceAlpha = 0.4;
    p.FaceColor = color;
    p.EdgeColor = 'none';
end

function plot_amplitude_histogram(amplitudes, color, dark_color)
    % Plot amplitude histogram
    if length(amplitudes) < 2
        text(0.5, 0.5, sprintf('Insufficient data\n(%d points)', length(amplitudes)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'red');
        return;
    end
    
    try
        if length(amplitudes) >= 10
            h = histfit(amplitudes, min(10, length(amplitudes)));
            hold on
            h(1).FaceColor = color;
            h(1).FaceAlpha = 0.3;
            h(2).Color = dark_color;
        else
            h = histogram(amplitudes, min(5, length(amplitudes)));
            hold on
            h.FaceColor = color;
            h.FaceAlpha = 0.3;
        end
    catch ME
        fprintf('    Warning: Failed to plot histogram: %s\n', ME.message);
    end
end

function finalize_subplot(cellline_name, value_type, legend_entries)
    % Finalize subplot formatting
    ax = gca;
    title(cellline_name, 'FontSize', 10, 'FontName', 'Arial', 'Interpreter', 'none');
    
    if value_type == 1 % Periods
        if strcmp(ax.Type, 'polaraxes')
            ax.ThetaTick = linspace(0, 180, 8);
            ax.ThetaTickLabel = round(linspace(22, 36, numel(ax.ThetaTick)), 1);
            ax.ThetaLim = [0 180];
        end
    end
    
    set(ax, 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
    if value_type == 2
        set(ax, 'XLimitMethod', 'padded');
    end
    
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial');
    end
    
    hold off;
end

function add_figure_labels(fig, xlabel_text)
    % Add overall figure labels
    han1 = axes(fig, 'visible', 'off');
    han1.Title.Visible = 'on';
    han1.XLabel.Visible = 'on';
    han1.YLabel.Visible = 'on';
    xlabel(han1, xlabel_text, 'FontWeight', 'bold', 'FontSize', 15);
    ylabel(han1, 'Count', 'FontWeight', 'bold', 'FontSize', 15);
end

function figure_name = create_figure_name(plot_type, value_type, suffix)
    % Create standardized figure names
    figure_name = sprintf('Figures/cwt_figures/%s_%s_%s', plot_type, value_type, suffix);
end

function create_overlay_plot(meandata, celllines, value_name, xlabel_text, value_type)
    % Create overlay plot combining all cell lines
    if isempty(meandata) || all(cellfun(@isempty, meandata(:)))
        fprintf('Warning: No meandata available for overlay plots in %s\n', value_name);
        return;
    end
    
    fig = figure;
    fig.Position = [1, 1, 1440, 821];
    
    valid_celllines_overlay = {};
    parameters = [];
    
    for k = 1:numel(celllines)
        
        % Check if we have data for this cell line
        if all(cellfun(@isempty, meandata(:, k)))
            continue;
        end
        
        % Combine data from all reporters
        combined_data = [];
        for rep = 1:size(meandata, 1)
            if ~isempty(meandata{rep, k})
                combined_data = [combined_data; meandata{rep, k}(:)];
            end
        end
        
        if length(combined_data) < 2
            continue;
        end
        
        valid_celllines_overlay{end+1} = celllines{k};
        
        % Plot based on value type
        if value_type == 1 % Periods - polar histogram
            plot_overlay_periods(combined_data);
        else % Amplitudes - regular histogram
            plot_overlay_amplitudes(combined_data);
        end
        
        % Calculate parameters
        if isempty(parameters)
            parameters = nan(2, numel(celllines));
        end
        col_idx = find(strcmp(celllines, celllines{k}));
        if ~isempty(col_idx)
            parameters(1, col_idx) = median(combined_data, 'omitnan');
            if mean(combined_data, 'omitnan') ~= 0
                parameters(2, col_idx) = std(combined_data, 'omitnan') / mean(combined_data, 'omitnan');
            end
        end
        
    end
    
    if ~isempty(valid_celllines_overlay)
        % Finalize overlay plot
        finalize_overlay_plot(fig, value_type, value_name, xlabel_text, valid_celllines_overlay);
        
        % Save parameters if available
        if ~isempty(parameters) && any(~isnan(parameters(:)))
            save_table = array2table(parameters, 'VariableNames', celllines);
            fprintf('Parameters table created for %s\n', value_name);
        end
    else
        fprintf('Warning: No valid cell lines for overlay plot in %s\n', value_name);
        close(fig);
    end
end

function plot_overlay_periods(periods)
    % Plot periods in overlay polar histogram
    valid_periods = periods(~isnan(periods) & periods >= 20 & periods <= 40);
    
    if length(valid_periods) < 2
        return;
    end
    
    angles = (valid_periods - 22) / (36 - 22) * pi;
    angles = angles(angles >= 0 & angles <= pi);
    
    if isempty(angles)
        return;
    end
    
    num_bins = min(5, max(3, length(angles)));
    p = polarhistogram(angles, num_bins, 'Normalization', 'pdf');
    hold on;
    p.FaceAlpha = 0.4;
    p.EdgeColor = 'none';
end

function plot_overlay_amplitudes(amplitudes)
    % Plot amplitudes in overlay histogram
    try
        if length(amplitudes) >= 10
            h2 = histfit(amplitudes, min(10, length(amplitudes)));
            hold on
            h2(1).FaceAlpha = 0;
            h2(1).EdgeColor = 'none';
            h2(1).HandleVisibility = 'off';
            bar_color = h2(1).FaceColor;
            h2(2).Color = bar_color;
        else
            h2 = histogram(amplitudes, min(5, length(amplitudes)));
            hold on
            h2.FaceAlpha = 0.7;
            h2.HandleVisibility = 'on';
        end
    catch ME
        fprintf('  Warning: Failed to create overlay plot: %s\n', ME.message);
    end
end

function finalize_overlay_plot(fig, value_type, value_name, xlabel_text, celllines)
    % Finalize overlay plot formatting
    if value_type == 1 % Periods
        ax = gca;
        if strcmp(ax.Type, 'polaraxes')
            ax.ThetaTick = linspace(0, 180, 8);
            ax.ThetaTickLabel = round(linspace(22, 36, numel(ax.ThetaTick)), 1);
            ax.ThetaLim = [0 180];
            title('Half-Polar Histogram of Periods (22–36 Hours)');
        end
        figure_name = create_figure_name('polarhistogramfits', value_name, 'nBins5_all_NB_celllines');
    else
        figure_name = create_figure_name('histogramfits', value_name, 'all_NB_celllines');
    end
    
    ax = gca;
    legend(celllines, 'Location', 'eastoutside', 'FontSize', 12, 'FontName', 'Arial');
    set(ax, 'LineWidth', 2, 'FontSize', 15, 'FontName', 'Arial');
    
    han1 = axes(fig, 'visible', 'off');
    han1.Title.Visible = 'on';
    han1.XLabel.Visible = 'on';
    han1.YLabel.Visible = 'on';
    title(han1, 'all NB cell lines', 'FontWeight', 'bold', 'FontSize', 18);
    xlabel(han1, xlabel_text, 'FontWeight', 'bold', 'FontSize', 18);
    ylabel(han1, 'Count', 'FontWeight', 'bold', 'FontSize', 18);
    
    hold off
    saveas(fig, figure_name, 'svg');
end