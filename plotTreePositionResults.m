

%accuracy k=1 across sessions

mice = {'Blossom', 'Bubbles', 'Buttercup', 'Daisy', 'Oksana', 'Poppy', 'Villanelle'};
for iMouse = 1:length(mice)
    mouse = mice{iMouse};
    positionpath =  ['Z:\luciede\NTNU\Codes\Wes\' mouse '\'];
    load([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat']);

    figure 
    for SI = [1 2 3 5 6 7 8 9 10]
        temp(:,SI) = TreeResults_CV_Compact.Scores(SI,:);
        temp(temp==0) = NaN;
        boxchart(temp)
        hold on
        clear temp
        title(mice{iMouse})
         xlabel('Session')
        ylabel('Scores')
        ylim([0 0.2])
    end


    figure
    for SI = [1 2 3 5 6 7 8 9 10]
        temp(:,SI) = TreeResults_CV_Compact.ScoresR(SI,:);
        temp(temp==0) = NaN;
        boxchart(temp)
        hold on
        clear temp
        title(mice{iMouse})
         xlabel('Session')
        ylabel('Scores Random')
    end


    figure
    for SI = [1 2 3 5 6 7 8 9 10]
        temp(:,SI) = TreeResults_CV_Compact.Depths(SI,:);
        temp(temp==0) = NaN;
        boxchart(temp)
        hold on
        clear temp
        title(mice{iMouse})
        xlabel('Session')
        ylabel('Tree depth')
        %ylim([0 0.2])
    end


    figure 

    for SI = [1 2 3 5 6 7 8 9 10]
        temp(:,SI) = TreeResults_CV_Compact.DepthsR(SI,:);
        temp(temp==0) = NaN;
        boxchart(temp)
        hold on
        clear temp
        title(mice{iMouse})
        xlabel('Session')
        ylabel('Tree depth Random')
    end
end


%% Pool the data across mice

for SI =  [1 2 3 5 6 7 8 9 10]

    for iMouse = 1:length(mice)
        mouse = mice{iMouse};
        positionpath =  ['Z:\luciede\NTNU\Codes\Wes\' mouse '\'];
        load([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat']);

        %Scores
        result1 = TreeResults_CV_Compact.Scores(SI,:);

        if iMouse > 1
            if size(result1,1) < size( Scores_mice{SI},2)
                result1(1,size( Scores_mice{SI},2)) = NaN;
            end

        end
        Scores_mice{SI}(iMouse,:)= result1;

        %Random scores
        result2 = TreeResults_CV_Compact.ScoresR(SI,:);

        if iMouse > 1
            if size(result2,1) < size( ScoresR_mice{SI},2)
                result2(1,size( ScoresR_mice{SI},2)) = NaN;
            end

        end
        ScoresR_mice{SI}(iMouse,:)= result2;

        %Depth
        result3 = TreeResults_CV_Compact.Depths(SI,:);

        if iMouse > 1
            if size(result3,1) < size(Depth_mice{SI},2)
                result3(1,size(Depth_mice{SI},2)) = NaN;
            end

        end
        Depth_mice{SI}(iMouse,:)= result3;

        %Random Depth
        result4 = TreeResults_CV_Compact.DepthsR(SI,:);

        if iMouse > 1
            if size(result4,1) < size(DepthR_mice{SI},2)
                result4(1,size(DepthR_mice{SI},2)) = NaN;
            end

        end
        DepthR_mice{SI}(iMouse,:)= result4;
    end

    clear TreeResults_CV_Compact result1 result2 result3 result4
end

%% Plot for Day 1, all mice all K

%Day 1, scores
figure
SI = 3;
for K = 1:11
    temp(:,K) = Scores_mice{SI}(:,K);
    temp(temp==0) = NaN;
    boxchart(temp, 'BoxFaceColor', 'black')
    hold on
    clear temp
end

%Day 1, scores random on same figure
SI = 3;
for K = 1:11
    temp(:,K) = ScoresR_mice{SI}(:,K);
    temp(temp==0) = NaN;
    boxchart(temp, 'BoxFaceColor', 'red')
    hold on
    clear temp
end
xlabel('K step')
ylabel('Error rate at Day 1')
title('Day 1')


%%%%%%same but as a line plot, not a boxplot%%%%
SI = 3;
for K = 1:11
    std_data(1,K) = std( Scores_mice{SI}(:,K), 'omitnan');
    std_random(1,K) = std( ScoresR_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(Scores_mice{SI}));
sem_random= std_random/sqrt(length(ScoresR_mice{SI}));


figure
y = nanmean(Scores_mice{SI});
y_sem = sem_data;
x = 1:length(y);
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on

y = nanmean(ScoresR_mice{SI});
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(Scores_mice{SI}), 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')
plot(nanmean(ScoresR_mice{SI}), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Error rate at Day 1')
title('Day 1')



%Are the scores significantly lower than the random, control data?
SI = 3;
for K = 1:11
    [h,p] = ttest2(Scores_mice{SI}(:,K),ScoresR_mice{SI}(:,K) );
    ttest_scores_D1(1,K) = h;
    ttest_scores_D1(2,K)=p;
    clear h p
end

%Are the scores significantly different across K, for the real data?
SI = 3;


[p,tbl,stats] = anova1(Scores_mice{SI});
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/anova_day1.csv')

%% Plot for Day 7, all mice all K

%Day 7, scores
figure
SI = 10;
for K = 1:11
    temp(:,K) = Scores_mice{SI}(:,K);
    temp(temp==0) = NaN;
    boxchart(temp, 'BoxFaceColor', 'black')
    hold on
    clear temp
end

%Day 7, scores random on same figure
SI = 10;
for K = 1:11
    temp(:,K) = ScoresR_mice{SI}(:,K);
    temp(temp==0) = NaN;
    boxchart(temp, 'BoxFaceColor', 'red')
    hold on
    clear temp
end
xlabel('K step')
ylabel('Error rate at Day 7')
title('Day 7')

%%%%%%same but as a line plot, not a boxplot%%%%
SI = 10;
for K = 1:11
    std_data(1,K) = std( Scores_mice{SI}(:,K), 'omitnan');
    std_random(1,K) = std( ScoresR_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(Scores_mice{SI}));
sem_random= std_random/sqrt(length(ScoresR_mice{SI}));


figure
y = nanmean(Scores_mice{SI});
x = 1:length(y);
y_sem = sem_data;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

y = nanmean(ScoresR_mice{SI});
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(Scores_mice{SI}), 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')
plot(nanmean(ScoresR_mice{SI}), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Error rate at Day 7')
title('Day 7')

%Are the scores significantly lower than the random, control data?
SI = 10;
for K = 1:11
    [h,p] = ttest2(Scores_mice{SI}(:,K),ScoresR_mice{SI}(:,K) );
    ttest_scores_D7(1,K) = h;
    ttest_scores_D7(2,K)=p;
    clear h p
end

%Are the scores significantly different across K, for the real data?
SI = 10;
[p,tbl,stats] = anova1(Scores_mice{SI});
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/anova_day7.csv')

%% Is it doing better at Day 1 or Day7?
%For each K, ttest between S3 and S10 
for K = 1:11
    [h,p] = ttest2(Scores_mice{3}(:,K),Scores_mice{10}(:,K), 'tail', 'right'); %Testing if error rate on day 1> error rate d7
    ttest_scores_D1vsD7(1,K) = h;
    ttest_scores_D1vsD7(2,K) = p;
end
%% Error rates all days on the same plot, real data

figure
for SI = [3 5 6 7 8 9 10]
    legends{SI} = num2str(SI);
for K = 1:11
    std_data(1,K) = std( Scores_mice{SI}(:,K), 'omitnan');
    std_random(1,K) = std( ScoresR_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(Scores_mice{SI}));
sem_random= std_random/sqrt(length(ScoresR_mice{SI}));

y = nanmean(Scores_mice{SI});
y_sem = sem_data;
x = 1:length(y);
%patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
plot(nanmean(Scores_mice{SI}), 'LineWidth', 2, 'marker', 'o')
hold on
ylim([0 0.5])
xlabel('K step')
ylabel('Error rate at Day 1')
title('Day 1')
clear y y_sem
legend =legends{SI};
end


%% Depth, Plot for Day 1, all mice all K

%Day 1, depth
figure
SI = 3;
for K = 1:11
    temp(:,K) = Depth_mice{SI}(:,K);
    temp(temp==0) = NaN;
    b = boxchart(temp, 'BoxFaceColor', 'black');
    b.MarkerColor = 'black';
    hold on
    clear temp b
end

%Day 1, scores random on same figure
SI = 3;
for K = 1:11
    temp(:,K) = DepthR_mice{SI}(:,K);
    temp(temp==0) = NaN;
    b =   boxchart(temp, 'BoxFaceColor', 'red')
    b.MarkerColor = 'red';
    hold on
    clear temp b
end
xlabel('K step')
ylabel('Max tree depth at Day 1')
title('Day 1')


%%%%%%same but as a line plot, not a boxplot%%%%
SI = 3;
for K = 1:11
    std_data(1,K) = std(Depth_mice{SI}(:,K), 'omitnan');
    std_random(1,K) = std( DepthR_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(Depth_mice{SI}));
sem_random= std_random/sqrt(length(DepthR_mice{SI}));

figure
y = nanmean(Depth_mice{SI});
y_sem = sem_data;
x = 1:length(y);
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

y = nanmean(DepthR_mice{SI});
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(Depth_mice{SI}), 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')
plot(nanmean(DepthR_mice{SI}), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 450])
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Max tree depth at Day 1')
title('Day 1')

%Are the max depth significantly lower than the random, control data?
SI = 3;
for K = 1:11
    [h,p] = ttest2(Depth_mice{SI}(:,K),DepthR_mice{SI}(:,K) );
    ttest_depth_D1(1,K) = h;
    ttest_depth_D1(2,K)=p;
    clear h p
end

%Are the scores significantly different across K, for the real data?
SI = 3;
[p,tbl,stats] = anova1(Depth_mice{SI});
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/anova_day1_depth.csv')


%% Depth, Plot for Day 7 all mice all K

%Day 7, depth
figure
SI = 10;
for K = 1:11
    temp(:,K) = Depth_mice{SI}(:,K);
    temp(temp==0) = NaN;
    b = boxchart(temp, 'BoxFaceColor', 'black');
    b.MarkerColor = 'black';
    hold on
    clear temp b
end

%Day 7, scores random on same figure
SI = 10;
for K = 1:11
    temp(:,K) = DepthR_mice{SI}(:,K);
    temp(temp==0) = NaN;
    b =   boxchart(temp, 'BoxFaceColor', 'red')
    b.MarkerColor = 'red';
    hold on
    clear temp b
end
xlabel('K step')
ylabel('Max tree depth at Day 7')
title('Day 7')



%%%%%%same but as a line plot, not a boxplot%%%%
SI = 10;
for K = 1:11
    std_data(1,K) = std(Depth_mice{SI}(:,K), 'omitnan');
    std_random(1,K) = std( DepthR_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(Depth_mice{SI}));
sem_random= std_random/sqrt(length(DepthR_mice{SI}));

figure
x = 1:11;
y = nanmean(Depth_mice{SI});
y_sem = sem_data;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

y = nanmean(DepthR_mice{SI});
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(Depth_mice{SI}), 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')
plot(nanmean(DepthR_mice{SI}), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 450])
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Max tree depth at Day 7')
title('Day 7')

%Are the max depth significantly lower than the random, control data?
SI = 10;
for K = 1:11
    [h,p] = ttest2(Depth_mice{SI}(:,K),DepthR_mice{SI}(:,K) );
    ttest_depth_D7(1,K) = h;
    ttest_depth_D7(2,K)=p;
    clear h p
end

%Are the scores significantly different across K, for the real data?
SI = 10;

clear p tbl stats
[p,tbl,stats] = anova1(Depth_mice{SI});
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/anova_day7_depth.csv')

%% Is it doing better at Day 1 or Day7? - depth
%For each K, ttest between S3 and S10 
for K = 1:11
    [h,p] = ttest2(Depth_mice{3}(:,K),Depth_mice{10}(:,K));
    ttest_depth_D1vsD7(1,K) = h;
    ttest_depth_D1vsD7(2,K) = p;
end
%% Find the sweet K spot (optimal K) for accuracy x efficiency 

%if optimal K not computed
mice = {'Blossom', 'Bubbles', 'Buttercup', 'Daisy', 'Oksana', 'Poppy', 'Villanelle'};
for iMouse = 1:length(mice)
    mouse = mice{iMouse};
    positionpath =  ['Z:\luciede\NTNU\Codes\Wes\' mouse '\'];
    load([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat']);

 %Loop through SI
 for SI = 1:10
     if  SI ==4
         optimal_K(SI) = NaN;
         continue
     end
     
     %normalise between 0 and 1
     score_norm = rescale(TreeResults_CV_Compact.Scores(SI,:), 0, 1);
    
     depth_norm = rescale(TreeResults_CV_Compact.Depths(SI,:), 0, 1);

     [min_value, K_index] = min(score_norm + depth_norm);

     optimal_K(SI) = K_index;

     clearvars -except TreeResults_CV_Compact SI optimal_K mice mouse positionpath iMouse opt_K_mice
 end

TreeResults_CV_Compact.optimal_K_OF = optimal_K;
save([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat'], 'TreeResults_CV_Compact')

opt_K_mice{iMouse} = optimal_K;
clearvars -except mice opt_K_mice iMouse
end

%Load optimal K if computed before 
mice = {'Blossom', 'Bubbles', 'Buttercup', 'Daisy', 'Oksana', 'Poppy', 'Villanelle'};
for iMouse = 1:length(mice)
    mouse = mice{iMouse};
    positionpath =  ['Z:\luciede\NTNU\Codes\Wes\' mouse '\'];
    load([positionpath 'TreeResults_S1S10_1515_separated_Position_CV_Compact.mat']);


opt_K_pos_mice{iMouse} = TreeResults_CV_Compact.optimal_K_OF;
clearvars -except mice  opt_K_behaviour_mice iMouse opt_K_pos_mice
end


% Plot the results of sweet K spot
figure
h = histogram(vertcat(opt_K_mice{:}), 'FaceColor', 'black')
xlabel('Optimal K step')
ylabel('Count')

%One plot per session
mat = vertcat(opt_K_mice{:});
for SI = [1 2 3 5 6 7 8 9 10]
    figure
    histogram( mat(:,SI))
    xlim([1 9])
    xlabel('Optimal K step')
    ylabel('Count')
    title(['Day ' num2str(SI)])

end

%% Compare OF to Day 1! is the addition of objects increasing decoder accuracy?

SI = 1;
for K = 1:11
    std_pos(1,K) = std( Scores_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_pos= std_pos/sqrt(length(Scores_mice{SI}));


figure
x = 1:11;
y = nanmean(Scores_mice{SI});
y_sem = sem_pos;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], [33/255 113/255 181/255] , 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

SI = 2 ;
for K = 1:11
    std_pos(1,K) = std( Scores_mice{SI}(:,K), 'omitnan');
end

%SEM
sem_pos= std_pos/sqrt(length(Scores_mice{SI}));

y = nanmean(Scores_mice{SI});
y_sem = sem_pos;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], [8/255 48/255 107/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

SI = 1;
plot(nanmean(Scores_mice{SI}), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
SI = 2;
plot(nanmean(Scores_mice{SI}), 'Color',  [8/255 48/255 107/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', [8/255 48/255 107/255])
ylim([0 0.5])
legend({'OF1', 'OF2'})
xlabel('K step')
ylabel('Error rate')
title('OF1 vs OF2')

for K = 1:11
    [h,p] = ttest2(Scores_mice{1}(:,K),Scores_mice{2}(:,K) );
    ttest_scores_1_2(1,K) = h;
    ttest_scores_1_2(2,K)=p;
    clear h p
end

%% Check for each mouse individually how the error rates evolves after introducing objects in the environment
for iMouse = 1:7

    figure
    x = 1:11;
    SI = 1;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
    hold on
    SI = 2;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [8/255 48/255 107/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', [8/255 48/255 107/255])
    ylim([0 0.5])
    legend({'OF1', 'OF2'})
    xlabel('K step')
    ylabel('Error rate')
    title(sprintf(' %s ', mice{iMouse}))


    figure
    x = 1:11;
    SI = 1;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
    hold on
    SI = 3;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [253/255 141/255 60/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [253/255 141/255 60/255])
    ylim([0 0.5])
    legend({'OF1', 'Day 1'})
    xlabel('K step')
    ylabel('Error rate')
    title(sprintf(' %s ', mice{iMouse}))


    figure
    x = 1:11;
    SI = 1;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
    hold on
    SI = 10;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [253/255 141/255 60/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [253/255 141/255 60/255])
    ylim([0 0.5])
    legend({'OF1', 'Day 7'})
    xlabel('K step')
    ylabel('Error rate')
    title(sprintf(' %s ', mice{iMouse}))

    figure
    x = 1:11;
    SI = 2;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
    hold on
    SI = 3;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [253/255 141/255 60/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', [253/255 141/255 60/255])
    ylim([0 0.5])
    legend({'OF2', 'Day 1'})
    xlabel('K step')
    ylabel('Error rate')
    title(sprintf(' %s ', mice{iMouse}))

    figure
    x = 1:11;
    SI = 2;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [33/255 113/255 181/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',  [33/255 113/255 181/255])
    hold on
    SI = 10;
    plot(Scores_mice{SI}(iMouse,:), 'Color',  [253/255 141/255 60/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', [253/255 141/255 60/255])
    ylim([0 0.5])
    legend({'OF2', 'Day 7'})
    xlabel('K step')
    ylabel('Error rate')
    title(sprintf(' %s ', mice{iMouse}))

end

%% Overall error rate at a single K (K1 and optimal Ks), across sessions, mice pooled

K = 1;

for SI = [1 2 3 5 6 7 8 9 10]
    to_plot_K1_data(:,SI) =Scores_mice{SI}(:,K);
    to_plot_K1_random(:,SI) =ScoresR_mice{SI}(:,K);
end
to_plot_K1_data(to_plot_K1_data==0) = NaN;
to_plot_K1_random(to_plot_K1_random==0) = NaN;


std_data= std( to_plot_K1_data, 'omitnan');
std_random = std(to_plot_K1_random, 'omitnan');
sem_data= std_data/sqrt(length(to_plot_K1_data));
sem_random= std_random/sqrt(length(to_plot_K1_random));

figure
y = nanmean(to_plot_K1_data);
y_sem = sem_data;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
x = 1:length(y);
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on

y = nanmean(to_plot_K1_random);
y_sem = sem_random;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

y = nanmean(to_plot_K1_data);
y(isnan(y)) = [];
plot(y, 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')

y = nanmean(to_plot_K1_random);
y(isnan(y)) = [];
plot(y, 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('Session')
xticklabels({'OF1', 'OF2','Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6', 'Day7'})
ylabel('Error rate across the experiment')
title('K1')



K = 7;
for SI = [1 2 3 5 6 7 8 9 10]
    to_plot_K7_data(:,SI) =Scores_mice{SI}(:,K);
    to_plot_K7_random(:,SI) =ScoresR_mice{SI}(:,K);
end
to_plot_K7_data(to_plot_K7_data==0) = NaN;
to_plot_K7_random(to_plot_K7_random==0) = NaN;

std_data= std( to_plot_K7_data, 'omitnan');
std_random = std(to_plot_K7_random, 'omitnan');
sem_data= std_data/sqrt(length(to_plot_K7_data));
sem_random= std_random/sqrt(length(to_plot_K7_random));

figure
y = nanmean(to_plot_K7_data);
y_sem = sem_data;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
x = 1:length(y);
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on

y = nanmean(to_plot_K7_random);
y_sem = sem_random;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

y = nanmean(to_plot_K7_data);
y(isnan(y)) = [];
plot(y, 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')

y = nanmean(to_plot_K7_random);
y(isnan(y)) = [];
plot(y, 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('Session')
xticklabels({'OF1', 'OF2','Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6', 'Day7'})
ylabel('Error rate across the experiment')
title('K7')



K = 8;
for SI = [1 2 3 5 6 7 8 9 10]
    to_plot_K8_data(:,SI) =Scores_mice{SI}(:,K);
    to_plot_K8_random(:,SI) =ScoresR_mice{SI}(:,K);
end
to_plot_K8_data(to_plot_K8_data==0) = NaN;
to_plot_K8_random(to_plot_K8_random==0) = NaN;


std_data= std( to_plot_K8_data, 'omitnan');
std_random = std(to_plot_K8_random, 'omitnan');
sem_data= std_data/sqrt(length(to_plot_K8_data));
sem_random= std_random/sqrt(length(to_plot_K8_random));

figure
y = nanmean(to_plot_K8_data);
y_sem = sem_data;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
x = 1:length(y);
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'red', 'FaceAlpha', 0.2, 'EdgeColor', 'none')
hold on

y = nanmean(to_plot_K8_random);
y_sem = sem_random;
y(isnan(y)) = [];
y_sem(isnan(y_sem)) = [];
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )],  'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

y = nanmean(to_plot_K8_data);
y(isnan(y)) = [];
plot(y, 'Color', 'red', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'red')

y = nanmean(to_plot_K8_random);
y(isnan(y)) = [];
plot(y, 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('Session')
xticklabels({'OF1', 'OF2','Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6', 'Day7'})
ylabel('Error rate across the experiment')
title('K8')

%anova across days
to_anova_K1_data =to_plot_K1_data;
to_anova_K1_data(:,4) =[];
[p,tbl,stats] = anova1(to_anova_K1_data);
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/K1_across_days.csv')

to_anova_K7_data =to_plot_K7_data;
to_anova_K7_data(:,4) =[];
[p,tbl,stats] = anova1(to_anova_K7_data);
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/K7_across_days.csv')

to_anova_K8_data =to_plot_K8_data;
to_anova_K8_data(:,4) =[];
[p,tbl,stats] = anova1(to_anova_K8_data);
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/K8_across_days.csv')


%% With all sessions averaged, at which K is real data as bad a shuffled data?


 to_concat_data = Scores_mice;
 to_concat_data{1} = [];
 to_concat_data{2} = [];
 to_plot_data = vertcat(to_concat_data{:});

 to_concat_random = ScoresR_mice;
 to_concat_random{1} = [];
 to_concat_random{2} = [];
 to_plot_random = vertcat(to_concat_random{:});


for K = 1:11
    std_data(1,K) = std( to_plot_data(:,K), 'omitnan');
    std_random(1,K) = std( to_plot_random(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(to_plot_data));
sem_random= std_random/sqrt(length(to_plot_random));


figure
y = nanmean(to_plot_data);
x = 1:length(y);
y_sem = sem_data;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], [140/255 81/255 10/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

y = nanmean(to_plot_random);
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(to_plot_data), 'Color',[140/255 81/255 10/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',[140/255 81/255 10/255])
plot(nanmean(to_plot_random), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
ylim([0 0.5])
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Error rate across sessions')


for K = 1:11
    [h,p] = ttest2(to_plot_data(:,K),to_plot_random(:,K) );
    ttest_score_alldays(1,K) = h;
    ttest_score_alldays(2,K)=p;
    clear h p
end

%Anova to compare between Ks
[p,tbl,stats] = anova1(to_plot_data);
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/ErrorRate_ave_days.csv')


%depth
clear to_concat_random to_concat_data to_plot_data to_plot_random


to_concat_data = Depth_mice;
to_concat_data{1} = [];
to_concat_data{2} = [];
to_plot_data = vertcat(to_concat_data{:});

to_concat_random = DepthR_mice;
to_concat_random{1} = [];
to_concat_random{2} = [];
to_plot_random = vertcat(to_concat_random{:});


for K = 1:11
    std_data(1,K) = std( to_plot_data(:,K), 'omitnan');
    std_random(1,K) = std( to_plot_random(:,K), 'omitnan');
end

%SEM
sem_data= std_data/sqrt(length(to_plot_data));
sem_random= std_random/sqrt(length(to_plot_random));


figure
y = nanmean(to_plot_data);
x = 1:length(y);
y_sem = sem_data;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], [140/255 81/255 10/255], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on

y = nanmean(to_plot_random);
y_sem = sem_random;
patch([x fliplr(x)], [y-y_sem fliplr(y+y_sem )], 'black', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(nanmean(to_plot_data), 'Color', [140/255 81/255 10/255], 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor',[140/255 81/255 10/255])
plot(nanmean(to_plot_random), 'Color', 'black', 'LineWidth', 2, 'marker', 'o', 'MarkerFaceColor', 'black')
legend({'Data', 'Shuffled'})
xlabel('K step')
ylabel('Max tree depth across sessions')


for K = 1:11
    [h,p] = ttest2(to_plot_data(:,K),to_plot_random(:,K) );
    ttest_depth_alldays(1,K) = h;
    ttest_depth_alldays(2,K)=p;
    clear h p
end

%Anova to compare between Ks
[p,tbl,stats] = anova1(to_plot_data);
results = multcompare(stats);
tbl = array2table(results,"VariableNames",   ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])
writetable(tbl, 'Z:\luciede\NTNU\Miniscope-rig\2025\CV Trees\Position/MaxDepth_ave_days.csv')

%% Changing visualisation - plotting only optimal K = K8

%error rate d1 vs d7
figure
SI = 3;
K =8;
temp(:,1) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'magenta')
hold on
clear temp


SI = 10;
K=8;
temp(:,2) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'cyan')
hold on
clear temp

xlabel('Session')
xticklabels({'Day 1', 'Day 7'})
ylabel('Error rate at Day 1 vs Day 7')

%depth d1 vs d7
figure
SI = 3;
K =8;
temp(:,1) = Depth_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'magenta')
hold on
clear temp


SI = 10;
K=8;
temp(:,2) = Depth_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'cyan')
hold on
clear temp

xlabel('Session')
xticklabels({'Day 1', 'Day 7'})
ylabel('Max depth at Day 1 vs Day 7')

%error rate Of2 vs D1 vs d7
figure

SI = 2;
K =8;
temp(:,1) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'green')
hold on
clear temp

SI = 3;
K =8;
temp(:,2) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'magenta')
hold on
clear temp


SI = 10;
K=8;
temp(:,3) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'cyan')
hold on
clear temp

xlabel('Session')
xticklabels({'OF2', 'Day 1', 'Day 7'})
ylabel('Error rate')

%error rate of1 vs Of2 vs d1 vs d7
figure

SI = 1;
K =8;
temp(:,1) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'yellow')
hold on
clear temp

SI = 2;
K =8;
temp(:,2) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'green')
hold on
clear temp

SI = 3;
K =8;
temp(:,3) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'magenta')
hold on
clear temp


SI = 10;
K=8;
temp(:,4) = Scores_mice{SI}(:,K);
temp(temp==0) = NaN;
boxchart(temp, 'BoxFaceColor', 'cyan')
hold on
clear temp

xlabel('Session')
xticklabels({'OF1', 'OF2', 'Day 1', 'Day 7'})
ylabel('Error rate')