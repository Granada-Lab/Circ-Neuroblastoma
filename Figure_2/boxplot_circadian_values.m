function boxplot_circadian_values

%Carolin Ector, 25.10.2023
%version 3, modified on 03.07.2024 to fit Neuroblastoma data

%Function sorts cell lines by unweighted circadian parameters and creates boxplots

%input: stored in "workspace_circadian_parameters.mat"
% pathtofolder: path to manuscript folder
% celllines_text: names of the cell lines being analysed, optimized for being used as x-labels
% colorbp: colors used in the graphs for the two circadian clock luciferase reporters, or the merged version (Bmal1+Per2)
colorbp = {[0.1216,0.4667,0.7059];[1,0.4980,0.0549]};
metric = {'ridgelength_unnorm'};
yaxisnames = {'Ridge length (hours)'};

for m = 1:1%numel(metric) %loop m metric

    fig = figure;
    fig.Position = [420,285,525,425];

    %load data
    file_bmal = 'extracted_circadian_parameters_by_replicate_Bmal1.xlsx';
    file_per = 'extracted_circadian_parameters_by_replicate_Per2.xlsx';

    [t_bmal] = readtable(file_bmal,'sheet',metric{m});
    [t_per] = readtable(file_per,'sheet',metric{m});
    bmal = table2array(t_bmal(:,2:end));
    per = table2array(t_per(:,2:end));
    datacombined = [bmal;per];

    celllines = t_bmal.Properties.VariableNames;
    celllines(:,1) =[];

    fig = figure;
    fig.Position = [420,285,525,425];

    % Replicate each data point by its weight
    
    for c = 1:size(datacombined,2)
        n_samples(:,c) = numel(rmmissing(datacombined(:,c)));
        xlabeltext(:,c) = append(celllines(:,c),' (',num2str(n_samples(:,c)),')');
    end

    %sort cell lines by their median lag or peak value (ascending)
    order = 'descend';
    Med = median(datacombined,1,'omitnan');
    nanColumns = any(isnan(Med), 1);
    bmal_clean = bmal(:, ~nanColumns); 
    per_clean = per(:, ~nanColumns);
    xlabeltext_clean = xlabeltext(:, ~nanColumns);
    Med_clean = Med(:, ~nanColumns);

    [~, sortIdx] = sort(Med_clean,order);
    bmalSorted = bmal_clean(:,sortIdx);
    perSorted = per_clean(:,sortIdx);
    celllinesSorted = xlabeltext_clean(:,sortIdx);
    dataSorted = [bmalSorted;perSorted];

    hold all

    %add datapoints as scatter
    datasize1 = size(bmalSorted);
    datasize2 = size(perSorted);
    x1 = repmat(1:datasize1(:,2),datasize1(:,1),1);
    x2 = repmat(1:datasize2(:,2),datasize2(:,1),1);
    s1 = scatter(x1,bmalSorted,67,colorbp{1},'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
    s2 = scatter(x2,perSorted,67,colorbp{2},'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
    set(s2,'MarkerEdgeColor','none');
    set(s1,'MarkerEdgeColor','none');

    %plot boxplot
    %     boxplot(dataDescend,'Labels',celllinesDescend,'PlotStyle','compact','Symbol','+r','MedianStyle','line'); hold all
    boxplot(dataSorted,'Labels',celllinesSorted);
    boxes = findobj(gca, 'Tag', 'Box');
    for j=1:length(boxes)
        set(boxes(j), 'Color', [1 1 1],'LineWidth',0.8); % Set edge color to green
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), [0.4 0.4 0.4], 'FaceAlpha', 0.2); % Set face color to red with alpha value of 0.3
    end

    hold off

    % Aesthetics
    ylabel(yaxisnames{m})
    hold off

    ax = gca;
    set(findobj(gca,'type','line'),'linew',1);
    numCategories = numel(celllinesSorted);
    set(ax, 'XTick', 1:numCategories, 'XTickLabel', celllinesSorted, 'XTickLabelRotation', 45);
    box on;
    grid on
    set(ax,'YLimitMethod', 'padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    figurename = 'boxplot_ridgelengthes.svg';
    saveas(fig,figurename);
    clear figurename

end %loop m metric
end %function