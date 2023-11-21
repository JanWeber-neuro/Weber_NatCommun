function plot_fig_1c(path)

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_1c');

%% Plot Reaction Time

avg_RT_0  = tbl.Var1(2:end);
avg_RT_25 = tbl.Var2(2:end);
avg_RT_75 = tbl.Var3(2:end);

% create cell array
plotRT = {avg_RT_0, avg_RT_25, avg_RT_75}';

% colors use for plotting
cb = {'#008450', '#EFB700', '#B81D13'};

figure(1); 
h = rm_raincloud(plotRT, [1 0 0], 0, 'ks');

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
xticks([0 0.15 0.25]);
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('RT [sec.]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);

%% Plot Accuracy

avg_ACC_0  = tbl.Var5(2:end);
avg_ACC_25 = tbl.Var6(2:end);
avg_ACC_75 = tbl.Var7(2:end);

% create cell array
plotAcc = {avg_ACC_0, avg_ACC_25, avg_ACC_75}';

figure(2); 
h = rm_raincloud(plotAcc, [1 0 0], 0, 'ks');

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('Accuracy');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);

%% Plot IQR

avg_iqr_0  = tbl.Var9(2:end);
avg_iqr_25 = tbl.Var10(2:end);
avg_iqr_75 = tbl.Var11(2:end);

% create cell array
plotIQR = {avg_iqr_0, avg_iqr_25, avg_iqr_75}';

% color to use for plotting
cb = {'#008450', '#EFB700', '#B81D13'};

figure(3); 
h = rm_raincloud(plotIQR, [1 0 0], 0, 'ks', 0.009);

% change density plot color
h.p{1,1}.FaceColor = cb{1};
h.p{2,1}.FaceColor = cb{2};
h.p{3,1}.FaceColor = cb{3};
h.p{1,1}.EdgeColor = cb{1};
h.p{2,1}.EdgeColor = cb{2};
h.p{3,1}.EdgeColor = cb{3};
h.p{1,1}.FaceAlpha = 0.8;
h.p{2,1}.FaceAlpha = 0.8;
h.p{3,1}.FaceAlpha = 0.8;

% change color for indiv. dots
h.s{1,1}.MarkerFaceColor = cb{1};
h.s{2,1}.MarkerFaceColor = cb{2};
h.s{3,1}.MarkerFaceColor = cb{3};
h.s{1,1}.MarkerFaceAlpha = 0.8;
h.s{2,1}.MarkerFaceAlpha = 0.8;
h.s{3,1}.MarkerFaceAlpha = 0.8;
h.s{1,1}.SizeData = 60;
h.s{2,1}.SizeData = 60;
h.s{3,1}.SizeData = 60;

% change mean dots
h.m(1).MarkerFaceColor = cb{1};
h.m(2).MarkerFaceColor = cb{2};
h.m(3).MarkerFaceColor = cb{3};
h.m(1).SizeData = 100;
h.m(2).SizeData = 100;
h.m(3).SizeData = 100;

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('Likelihood of Stop'); % x and y axis are flipped in the function
xlabel('IQR (RT)');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [425 276 290 202]);
