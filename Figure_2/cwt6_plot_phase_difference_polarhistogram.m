function cwt6_plot_phase_difference_polarhistogram

%Carolin Ector, 31.10.2023

%Function plots weighted Bmal1-Per2-phase differences per cell line and per subtype in polarhistograms

load('Results/extracted_phase_differences.mat')

disp('cwt6_plot_phase_difference_polarhistogram.m is executed')

celllines = {'GIMEN', 'NGP', 'SKNSH'};
total_time = (0:0.16666667:137.7)'; %total recording time

%create figure for polarhistograms
fig1 = figure;%('Visible','off');
fig1.Position = [1,1,2560,1361];
hold all

for c = 1:numel(celllines)

    phdiff = phdiff_all{c};
    time = time_all{c};

    row_indices = [];
    g = 0;

    for h = 1:numel(phdiff)
        test = isnan(phdiff{h});
        if test == 1
            g = g+1;
            row_indices(g,:) = h;
        end
        clear test
    end

    if ~isempty(row_indices)
        phdiff(row_indices,:) = [];
        time(row_indices,:) = [];
    end

    u = 0;
    for i = 1:numel(phdiff)
        p = phdiff{i};
        t = time{i};
        p1 = isnan(p);
        if p1 ~= 1
            u = u+1;
            if t(1,:) ~= 0
                coltoadd = numel(0:0.1666666667:(t(1,:)-0.1));
                emptycols(1:coltoadd,:) = NaN;
                t = [emptycols;t];
                p = [emptycols;p];
                clear coltoadd
                clear emptycols
            end

            le_total = length(total_time);

            if length(t) < le_total
                t((end+1:le_total),:) = NaN;
                p((end+1:le_total),:) = NaN;
            end

            TF1 = isnan(t);

            x(:,u) = p;

            varstoclear6 = {'t','p'};
            clear(varstoclear6{:})
        else

        end
        clear p1
    end

    % Initialize variables to store the weighted sum and count
    sum = zeros(size(x, 1), 1);
    count = zeros(size(x, 1), 1);

    % Loop through each row in x
    for row = 1:size(x, 1)
        % Initialize variables to track the sum and count for the current row
        row_sum = 0;
        row_count = 0;

        % Loop through each time series
        for col = 1:size(x, 2)
            % Check if the value is not NaN
            if ~isnan(x(row, col))
                % Update the sum and count for the current row
                row_sum = row_sum + x(row, col);
                row_count = row_count + 1;
            end
        end

        % Check if there are at least two non-NaN values in the current row
        if row_count >= 3
            % Update the weighted sum and count for the final result
            sum(row) = row_sum;
            count(row) = row_count;
        end
    end

    % Calculate the weighted mean time-series
    meanphasespercellline(:,c) = sum ./ count;

    TF = isnan(meanphasespercellline(:,c));
    x2 = total_time;
    x2(TF,:) = NaN;

    varstoclear3 = {'x','phdiff2','subtitles','times','commonridge1','commonridge2','w','weights2','row_indices'};
    clear(varstoclear3{:});

    %% phase difference over time per cell line - weighted mean of all replicates (weight = common ridge length)
    subplot(1,3,c);

    plot(x2,meanphasespercellline(:,c));

    ylim([-2*pi,2*pi]);
    yticks(-2*pi:pi:2*pi);
    yticklabels({'-2π','-π','0','π','2π'});
    xticks(0:24:144);
    xticklabels({'0','24','48','72','96','120','144'});
    title(celllines{c},'FontSize',10,'FontName','Arial','Interpreter','none');
    ax=gca;
    set(ax,'XLimitMethod','padded','linewidth',1.5,'FontSize',10,'FontName','Arial');

    clear x2
    clear y
    clear w

end %celllines

% resultexcel_x = append('continuous_phdiff_by_cellline_mean.xlsx');
% table = array2table(meanphasespercellline,'VariableNames',celllines);
% writetable(table,resultexcel_x,'Sheet','weighted_mean_phdiff');

han1=axes(fig1,'visible','off');

han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
xlabel(han1,'Time (h)','FontWeight','bold','FontSize',11);
ylabel(han1,'Phase difference (radians)','FontWeight','bold','FontSize',11);

hold off

%save figure
figurename1 = append('Figures/cwt_figures/phase_difference_over_time_avgpercellline.svg');
saveas(fig1, figurename1);

%% polarhistogram per subtype, TNBC subtypes separate -> Figure 2E

fig = figure;%('Visible','off');
fig.Position = [1,1,2560,1361];
aa = 0;

for k = 1:numel(celllines) %loop k subtypes

    subplot(1,3,k);

    aa = aa+1;
    %load data for specific combination
    alpha2 = rmmissing(meanphasespercellline(:,k));

    %plot polarhistogtam
    polarhistogram(alpha2,6,'Normalization','probability','FaceAlpha',0.4); %v2
    hold all

    %calculate parameters:
    parameters(1,k) = circ_median(alpha2);
    parameters(2,k) = circ_var(alpha2);

    %add legend and modify appearance of the polarhistogram
    ax = gca;  % Get the current axes
    legend(celllines{k},'Location','eastoutside','FontSize',15,'FontName','Arial');
    set(ax,'linewidth',3,'FontSize',20,'FontName','Arial');
    ax.ThetaZeroLocation = 'top';  % Adjust as needed

    varstoclear4 = {'alpha2'};
    clear(varstoclear4{:});


end %loop celllines per subtype

savetable = array2table(parameters,'variablenames',celllines);

hold off
clear celllines

%save figure
figurename3 = append('Figures/cwt_figures/polarhistogram_phase_difference_percellline.svg');
saveas(fig, figurename3);

disp('cwt6_plot_phase_difference_polarhistogram.m is completed')

end %function
