function plot_fig_2ab_PeakMetrics_HFA(path)

% grand average structure
ga_0  = struct();
ga_25 = struct();
ga_75 = struct();

% colors
cb = {'#008450', '#EFB700', '#B81D13'};

%% Plot peak latency frontal

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2a upper left');

ga_0.pfc.latency  = table2array(tbl(2:end,1));
ga_25.pfc.latency = table2array(tbl(2:end,2));
ga_75.pfc.latency = table2array(tbl(2:end,3));

figure; 
h = rm_raincloud({ga_0.pfc.latency, ga_25.pfc.latency, ...
    ga_75.pfc.latency}', [1 0 0], 0, 'ks', 0.03);

% admin
for Icond = 1:3
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('Time to HLL [s]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [441 334 266 126]);
xlim([-0.5 0.3]);

view([0 90]);

clear tbl

%% Plot peak latency motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2b upper left');

ga_0.motor.latency  = table2array(tbl(2:end,1));
ga_25.motor.latency = table2array(tbl(2:end,2));
ga_75.motor.latency = table2array(tbl(2:end,3));

figure; 
h = rm_raincloud({ga_0.motor.latency, ga_25.motor.latency, ...
    ga_75.motor.latency}', [1 0 0], 0, 'ks', 0.03);

% admin
for Icond = 1:3
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('Time to HLL [s]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [441 334 266 126]);
xlim([-0.5 0.3]);

view([0 90]);

%% Plot peak amplitude frontal

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2a lower right');

ga_0.pfc.amplitude  = table2array(tbl(2:end,1));
ga_25.pfc.amplitude = table2array(tbl(2:end,2));
ga_75.pfc.amplitude = table2array(tbl(2:end,3));

figure; 
h = rm_raincloud({ga_0.pfc.amplitude, ga_25.pfc.amplitude, ...
    ga_75.pfc.amplitude}', [1 0 0], 0, 'ks', 0.9);

% admin
for Icond = 1:3
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('HFA Peak [z]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [544 272 186 189]);

%% Plot peak amplitude motor cortex

% load data 
tbl = readtable(fullfile(path.data, 'SourceData.xlsx'), 'sheet', 'Fig_2b lower right');

ga_0.motor.amplitude  = table2array(tbl([2:8 10:end],1));
ga_25.motor.amplitude = table2array(tbl([2:8 10:end],2));
ga_75.motor.amplitude = table2array(tbl([2:8 10:end],3));

figure; 
h = rm_raincloud({ga_0.motor.amplitude, ga_25.motor.amplitude, ...
    ga_75.motor.amplitude}', [1 0 0], 0, 'ks', 1);

% admin
for Icond = 1:3
    
    % change density plot color
    h.p{Icond,1}.FaceColor        = cb{Icond};
    h.p{Icond,1}.EdgeColor        = cb{Icond};
    h.p{Icond,1}.FaceAlpha        = 0.8;

    % change color for indiv. dots
    h.s{Icond,1}.MarkerFaceColor  = cb{Icond};
    h.s{Icond,1}.MarkerFaceAlpha  = 0.8;
    h.s{Icond,1}.SizeData         = 40;

    % change mean dots
    h.m(Icond).MarkerFaceColor    = cb{Icond};
    h.m(Icond).SizeData           = 80;
    
end

h.l(1).Color    = [0.7 0.7 0.7];
h.l(2).Color    = [0.7 0.7 0.7];

% general figure settings
tmpxticks = yticks; % due to rotation is different
yticks(tmpxticks);
xlim([3 15]);
yticklabels(fliplr({'0%', '25%', '75%'}));
ylabel('LL Stop'); % x and y axis are flipped in the function
xlabel('HFA Peak [z]');
set(gca, 'linewidth', 3, 'fontsize', 13);
box off
set(gcf, 'Position', [544 272 186 189]);
