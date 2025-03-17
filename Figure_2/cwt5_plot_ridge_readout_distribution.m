function cwt5_plot_ridge_readout_distribution(reporter_colors,reporter_colors_dark)

%Carolin Ector, 02.11.2023

%Function plots distribution of continuous amplitudes and periods per cellline and per subtype

%define parameters
celllines = {'IMR5','GIMEN','Kelly','CLBGA','Lan5','NGP','SKNBE','SKNSH'};
values = {'period';'amplitude'};
xaxisvalues = {'Period (hours)';'Amplitude (a.u.)'};
% weightsheets = {'weight_ridgelength_norm';'weight_ridgelength_unnorm'};
finer_subtype = {'TNBC-BL1';'TNBC-BL2';'TNBC-M';'LumA';'Epithelial';'Sarcoma'};
finer_sub = [4,3,4,3,1,1]; %for loop finer subtypes
reporterlegend = {'Bmal1','Bmal1 fit','Per2','Per2 fit'};

disp('cwt5_plot_ridge_readout_distribution.m is executed')

load('Results/sorted_periods_amplitudes.mat')

for v = 1:numel(values) %loop v values

    fig = figure;%('Visible','off');
    fig.Position = [1,1,1440,821];

    %calculate number of subplots per figure (including all cell models together)

    div = 3;
    n = ceil([numel(celllines)/div]);

    %load sorted continuous periods or amplitdes
    sorted_data = sorted_periods_amplitudes{v};
    sorted_data(:,[2,3,4]) = [];

    for c = 1:numel(celllines) % loop c celllines

        aa1 = 1;
        aa2 = 2;

        %load cell line data
        celllinedata = sorted_data{:,c};
        legend_entries1 = {}; % Collect cell line names for this subtype

        numRows = size(celllinedata, 1);
        for rowIdx = 1:numRows
            % Create a logical vector for the row: true if the cell is empty.
            rowEmptyCheck = cellfun(@isempty, celllinedata(rowIdx, :));

            % Check if all cells in the row are empty.
            if all(rowEmptyCheck)
                fprintf('In sorted_data{%d}, row %d is completely empty.\n', i, rowIdx);
                if rowIdx == 1
                    aa1 = 2;
                elseif rowIdx == 2
                    aa2 = 1;
                end
            end
        end

        for a = aa1:aa2 %loop a reporters

            %load data of specific reporter
            data = cell2mat(celllinedata(a,:));

            legend_entries1{end+1} = [reporterlegend{a},reporterlegend{a+1}];

            % %load weights (calculated from fitting a sigmoidal function done in "weighted_boxplot_circadian_values.m"
            % inputfile_weights = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
            % [weights] = table2array(readtable(inputfile_weights, 'Sheet',weightsheets{v}));
            % w = rmmissing(transpose(weights(:,c+1)));

            %% plot histogram per cell line model, overlay reporters

            for f = 1:size(data,1)
                meanval = mean(data(f,:),'omitnan');
                meantimeseries(f,c) = meanval;
            end

            subplot(div,n,c);
            x = meantimeseries(:,c);

            if v == 1 %plot periods in polarhistogram

                periods = x;

                % Convert periods to angles (Mapping 24–36 hours to 0–π radians)
                angles = (periods - 22) / (36 - 22) * pi;

                % Define number of bins
                numBins = 10;

                % Create half-polar histogram
                p = polarhistogram(angles, numBins, 'Normalization', 'pdf'); % Normalized histogram
                hold on;

                % Adjust visualization
                p.FaceAlpha = 0.4; % Make bars semi-transparent
                p.FaceColor = reporter_colors{a};
                p.EdgeColor = 'none'; % Remove bar edges

            else

                %plot bars
                h = histfit(x,3); hold on
                h(1).FaceColor = reporter_colors{a};
                h(1).FaceAlpha = 0.3;
                h(2).Color = reporter_colors_dark{a};

            end

            varstoclear = {'data','w','weighted_sum','count','x','x_values','y_values'};
            clear(varstoclear{:})

            meandata{a,c} = meantimeseries(:,c);
            % resultexcel = append('mean_by_cellline_',values{v},'_',reporters{a},'.xlsx');
            % outputsheet = celllines{c};
            % writematrix(meantimeseries(:,c),resultexcel,'sheet',outputsheet);

        end %loop a reporters

        if v == 1 %periods

            % Adjust axis ticks and labels
            ax = gca;
            ax.ThetaTick = linspace(0, 180, 8); % Set half-circle ticks (0° to 180°)
            ax.ThetaTickLabel = round(linspace(22, 36, numel(ax.ThetaTick)), 1); % Label in hours
            ax.ThetaLim = [0 180]; % Restrict to half-circle

            title('Half-Polar Histogram of Periods (22–36 Hours)');
            set(ax,'linewidth',1.5,'FontSize',12,'FontName','Arial');
            figurename1 = append('Figures/cwt_figures/polarhistogram_',values{v},'_avgpercellline');
        end

        hold off
        ax = gca;  % Get the current axes
        title(celllines{c},'FontSize',10,'FontName','Arial','Interpreter','none');

        if v == 2
            set(ax,'XLimitMethod','padded','linewidth',1.5,'FontSize',12,'FontName','Arial');
            figurename1 = append('Figures/cwt_figures/histogram_',values{v},'_avgpercellline');
        end
        
        legend(legend_entries1,'Location','northeast','FontSize',10,'FontName','Arial');

    end %loop c celllines

    han1=axes(fig,'visible','off');

    han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
    xlabel(han1,xaxisvalues{v},'FontWeight','bold','FontSize',15);
    ylabel(han1,'Count','FontWeight','bold','FontSize',15);

    hold off
    saveas(fig, figurename1, 'svg');

    %% overlay results all cell lines

    fig = figure;%('visible','off');
    fig.Position = [1,1,1440,821];

    for k = 1:numel(celllines) %loop k celllines

        x1 = mean(cell2mat(meandata(:,k)'),2,'omitnan');

        %plot histogtam
        % numBins = 10; % Increase this value for finer granularity
        % h2 = histfit(x1, numBins, 'kernel'); hold on

        if v == 1 %plot periods in polarhistogram

            periods =x1;

            % Convert periods to angles (Mapping 24–36 hours to 0–π radians)
            angles = (periods - 22) / (36 - 22) * pi;

            % Define number of bins
            numBins = 5;

            % Create half-polar histogram
            p = polarhistogram(angles, numBins, 'Normalization', 'pdf'); % Normalized histogram
            hold on;

            % Adjust visualization
            p.FaceAlpha = 0.4; % Make bars semi-transparent
            p.EdgeColor = 'none'; % Remove bar edges
            
        else %plot amplitudes in histogram

            h2 = histfit(x1,3); hold on
            h2(1).FaceAlpha = 0;
            h2(1).EdgeColor = 'none';
            h2(1).HandleVisibility = 'off';
            barColor = h2(1).FaceColor;
            h2(2).Color = barColor;

        end

        %calculate parameters:
        parameters(1,k) = median(x1,'omitnan');
        parameters(2,k) = std(x1,'omitnan')/mean(x1,'omitnan');

        varstoclear1 = {'x1'};
        clear(varstoclear1{:});

    end
  
    savetable = array2table(parameters,'variablenames',celllines);

    figurename2 = append('Figures/cwt_figures/histogramfits_',values{v},'_all_NB_celllines');

    if v == 1 %periods
        % Adjust axis ticks and labels
        ax = gca;
        ax.ThetaTick = linspace(0, 180, 8); % Set half-circle ticks (0° to 180°)
        ax.ThetaTickLabel = round(linspace(22, 36, numel(ax.ThetaTick)), 1); % Label in hours
        ax.ThetaLim = [0 180]; % Restrict to half-circle
        title('Half-Polar Histogram of Periods (22–36 Hours)');
        figurename2 = append('Figures/cwt_figures/polarhistogramfits_',values{v},'_nBins5_all_NB_celllines');
    end

    %add legend and modify appearance of the polarhistogram
    ax = gca;  % Get the current axes
    legend(celllines,'Location','eastoutside','FontSize',12,'FontName','Arial');
    set(ax,'linewidth',2,'FontSize',15,'FontName','Arial');

    han1=axes(fig,'visible','off');
    han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
    title(han1,'all NB cell lines','FontWeight','bold','FontSize',18);
    xlabel(han1,xaxisvalues{v},'FontWeight','bold','FontSize',18);
    ylabel(han1,'Count','FontWeight','bold','FontSize',18);

    hold off
    clear legend_entries

    saveas(fig, figurename2, 'svg');

    clear inputdata

end %loop v values

clear input
disp('cwt5_plot_ridge_readout_distribution.m is completed')

end %function
