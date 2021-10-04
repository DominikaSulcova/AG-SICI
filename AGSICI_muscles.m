%% AG-SICI: MUSCULAR CONTRIBUTION TO TEPs
% Written by Dominika for GABA-AD project (2021)
% 
% ----- extracts and quantifies muscular activity -----
% 1) prepares the data
%    - loads preprocessed data from letswave
%       - only muslce-related ICs kept from second ICA
%       - hilbert transformed
%       - baseline corrected [-0.2 -0.005]s
%    - crops the data in predefined window
%    - 
% 2) 

%% parameters
clear all; clc

% dataset
subject = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        a = 1;
        for p = 1:length(position)
            for c = 1:length(current)
               colours(a, :) = uisetcolor; 
               a = a + 1;
            end
        end
    case 'NO'
        if exist('colours.mat') > 0
            load('colours.mat')
        else
            disp('No colour scheme found in this directory!')    
        end
end
save('colours.mat', 'colours'); 
clear a c p answer

% load a random header
load('avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 06 across reversed stim_140.lw6', '-mat')

% visualization 
figure_counter = 1;
xstep = header.xstep; 
time_window = [-0.05, 0.3]; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

% create output folder
output_folder = [pwd '\AG-SICI_muscles'];
if ~exist(output_folder) 
    mkdir(pwd, 'AG-SICI_muscles')
end    

% identify matlab file with results
results_path = uigetdir(pwd, 'Choose the directory with the global output file');
results_folder = [results_path '\AG-SICI_plus.mat'];

%% 1) prepare the data
% ----- adjustable parameters -----
prefix_m = 'icfilt ica art-sup raw P1'; 
prefix_IC = 'icfilt_N45 ica2 avg avgchan bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
% ----- adjustable parameters -----

% load from lw
% AGSICI_muscle_activity = struct;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % define subject
                if subject(s) < 10
                    subj = ['0' num2str(subject(s))];
                else
                    subj = num2str(subject(s));
                end
                
                % load and crop - muscular contraction 
                name = [prefix_m ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
                load(name)
                AGSICI_muscle_activity.muscle.data(p, c, i, s, :, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));
                
                % load and crop - N45 IC 
                name = [prefix_IC ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
                load(name)
                AGSICI_muscle_activity.N45.data(p, c, i, s, :, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));            
            end 
        end
    end
end
clear p c i s name prefix_m prefix_IC subj

%% 2) calculate GFP
% calculate GFP for each subject/condition
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)  
                % calculate GFP
                AGSICI_muscle_activity.muscles.GFP(p, c, i, s, :) = std(squeeze(AGSICI_muscle_activity.muscles.data(p, c, i, s, :, :)), 1);
                AGSICI_muscle_activity.N45.GFP(p, c, i, s, :) = std(squeeze(AGSICI_muscle_activity.N45.data(p, c, i, s, :, :)), 1);
            end
        end
    end
end
clear p c i s 

% calculate mean GFP per condition
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for t = 1:size(AGSICI_muscle_activity.muscles.GFP, 5)    
                % muscle activity
                AGSICI_muscle_activity.muscles.GFP_mean(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.muscles.GFP(p, c, i, :, t)));
                
                % muscle activity
                AGSICI_muscle_activity.N45.GFP_mean(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.N45.GFP(p, c, i, :, t)));
            end
        end
    end
end
clear p c i t

% append new variables to the general MATLAB file
save(results_folder, 'AGSICI_muscle_activity', '-append');

%% 3) plot mean GFP
% ----- adjustable parameters -----
TOI_muscles = [0.010, 0.050];
TOI_N45 = [0.030, 0.060];
% ----- adjustable parameters -----

% plot mean GFP per condition
for p = 1:length(position)
    for c = 1:length(current)
        % launch the figure
        fig = figure(figure_counter);        
        sgtitle(sprintf('%s STS - %s current', position{p}, current{c}), 'FontSize', 16, 'FontWeight', 'bold')
        
        % plot GFP of muscular activity, mark TOI
        subplot(2, 1, 1)
        hold on
        rectangle('Position', [TOI_muscles(1), 0.01, TOI_muscles(2) - TOI_muscles(1), 5.5-0.01], 'FaceColor', [0.98, 0.83, 0.83], 'EdgeColor', 'none')
        plot_TEP(x, squeeze(AGSICI_muscle_activity.muscles.GFP_mean(p, c, :, :)), ...
            'limit', [0 5.5], 'legend', {'100 %rMT' '120 %rMT' '140 %rMT'})
        title('GFP - muscular activity')
        
        % plot GFP of N45 component, mark TOI
        subplot(2, 1, 2)
        hold on
        rectangle('Position', [TOI_N45(1), 0.01, TOI_N45(2) - TOI_N45(1), 1-0.01], 'FaceColor', [0.98, 0.83, 0.83], 'EdgeColor', 'none')
        plot_TEP(x, squeeze(AGSICI_muscle_activity.N45.GFP_mean(p, c, :, :)), 'limit', [0 1])
        title('GFP - N45 component')    
        
        % save figure
        figure_name = sprintf('AGSICI_muscles_GFP_%s_%s', position{p}, current{c});
        savefig([output_folder '\' figure_name '.fig'])
        saveas(fig, [output_folder '\' figure_name '.png'])
        
        % update counter
        figure_counter = figure_counter + 1;
    end
end
clear p c fig

%% 4) extract mean GFP values, export for R 
% ----- adjustable parameters -----
TOI_muscles = [0.010, 0.050];
TOI_N45 = [0.030, 0.060];
% ----- adjustable parameters -----

% set cropping limits for TOIs
x_start_muscles = (TOI_muscles(1) - time_window(1))/xstep;
x_end_muscles = (TOI_muscles(2) - time_window(1))/xstep;
x_start_N45 = (TOI_N45(1) - time_window(1))/xstep;
x_end_N45 = (TOI_N45(2) - time_window(1))/xstep;

% calculate mean GFP over TOI for each subject/condition 
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)  
                AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, s) = mean(squeeze(AGSICI_muscle_activity.muscles.GFP(p, c, i, s, x_start_muscles:x_end_muscles)));
                AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, s) = mean(squeeze(AGSICI_muscle_activity.N45.GFP(p, c, i, s, x_start_N45:x_end_N45)));
            end
        end
    end
end
clear p c i s x_start_muscles x_end_muscles x_start_N45 x_end_N45

% export as a long-format table
AGSICI_muscles = table;
AGSICI_muscles.subject = zeros(0); AGSICI_muscles.orientation = {}; AGSICI_muscles.intensity = zeros(0); 
AGSICI_muscles.GFP_muscle = zeros(0); AGSICI_muscles.GFP_N45 = zeros(0); 
row_cnt = 1;
for s = 1:length(subject)
    for p = 1:length(position)
        for c = 1:length(current)
            for i = 1:length(intensity)
                % fill in the row
                AGSICI_muscles.subject(row_cnt) = subject(s);             
                AGSICI_muscles.orientation(row_cnt) = {[position{p} '-' current{c}]};
                AGSICI_muscles.intensity(row_cnt) = str2double(intensity{i}(end-2:end));
                AGSICI_muscles.GFP_muscle(row_cnt) = AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, s);
                AGSICI_muscles.GFP_N45(row_cnt) = AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, s);
                
                % update row counter
                row_cnt = row_cnt + 1;
            end
        end
    end
end
writetable(AGSICI_muscles, [output_folder '\AGSICI_muscles.csv'])
clear p c s i row_cnt

% append new variables to the general MATLAB file
save(results_folder, 'AGSICI_muscle_activity', 'AGSICI_muscles', '-append');

%% 5) plot mean GFP over TOIs
% ----- muscular activity -----
% launch the figure
fig = figure(figure_counter); 
hold on

% plot the data
counter = 1;
for p = 1:length(position)
    for c = 1:length(current)
        % calculate mean and SEM
        for i = 1:length(intensity)            
            y(i) = mean(AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, :));
            SEM(i) = std(AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, :)) / sqrt(length(subject));
        end

        % plot
        perr(counter) = errorbar(1:length(intensity), y, SEM);

        % adjust parameters
        perr(counter).Color = colours(counter, :);
        perr(counter).LineWidth = 1.5;
        perr(counter).Marker = 'o';
        perr(counter).MarkerFaceColor = colours(counter, :);
        perr(counter).MarkerSize = 10;

        % update counter
        counter = counter + 1;
    end
end
clear p c i counter

% add features, adjust parameters
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
title(sprintf('Muscular activity (mean GFP): %d - %dms', TOI_muscles(1) * 1000, TOI_muscles(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('GFP (\muV \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);

% save the figure       
savefig([output_folder '\AGSICI_muscles_all.fig'])
saveas(fig, [pwd '\AGSICI_muscles_all.png'])

% update the counter
figure_counter = figure_counter + 1;  
clear fig leg perr y SEM

% ----- N45 -----
% launch the figure
fig = figure(figure_counter); 
hold on

% plot the data
counter = 1;
for p = 1:length(position)
    for c = 1:length(current)
        % calculate mean and SEM
        for i = 1:length(intensity)            
            y(i) = mean(AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, :));
            SEM(i) = std(AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, :)) / sqrt(length(subject));
        end

        % plot
        perr(counter) = errorbar(1:length(intensity), y, SEM);

        % adjust parameters
        perr(counter).Color = colours(counter, :);
        perr(counter).LineWidth = 1.5;
        perr(counter).Marker = 'o';
        perr(counter).MarkerFaceColor = colours(counter, :);
        perr(counter).MarkerSize = 10;

        % update counter
        counter = counter + 1;
    end
end
clear p c i counter

% add features, adjust parameters
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
title(sprintf('Component N45 (mean GFP): %d - %dms', TOI_N45(1) * 1000, TOI_N45(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('GFP (\muV \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);

% save the figure       
savefig([output_folder '\AGSICI_N45_all.fig'])
saveas(fig, [pwd '\AGSICI_N45_all.png'])

% update the counter
figure_counter = figure_counter + 1;  
clear fig leg perr y SEM

%% 6) overall correlation
% extract data
data_corr_muscle = [];
data_corr_N45 = [];
for p = 1:length(position)
    for c = 1:length(current)        
        for i = 1:length(intensity)
            data_corr_muscle = [data_corr_muscle; squeeze(AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, :))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, :))];
        end
    end
end
data_corr = [data_corr_muscle data_corr_N45];
clear p c i data_corr_muscle data_corr_N45

% choose colours
marker_col = [];
for a = 1:4
    for b = 1:57
        marker_col = [marker_col; colours(a, :)];
    end
end
clear a b

% ----- linear correlation -----
% prepare linear model: y ~ 1 + x
data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {'muscular activity' 'N45'});

% plot data + regression line
fig = figure(figure_counter);
hold on
plot_corr(data_model, data_corr, marker_col, 'Pearson')
title('Linear correlation', 'FontWeight', 'bold', 'FontSize', 16)

% save the figure       
savefig([output_folder '\AGSICI_muscles_corr_all.fig'])
saveas(fig, [output_folder '\AGSICI_muscles_corr_all.png'])

% update the counter
figure_counter = figure_counter + 1;  

% ----- non-linear correlation -----
% rank the data
for a = 1:size(data_corr, 2)
    [temp, data_corr_ranked(:, a)]  = ismember(data_corr(:, a), unique(data_corr(:, a)));
end
clear a temp

% prepare linear model: y ~ 1 + x
data_model_ranked = fitlm(data_corr_ranked(:, 1), data_corr_ranked(:, 2), 'VarNames', {'muscular activity' 'N45'});

% plot data + regression line
fig = figure(figure_counter);
hold on
plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')
title('Non-linear correlation: ranked data', 'FontWeight', 'bold', 'FontSize', 16)

% save the figure       
savefig([output_folder '\AGSICI_muscles_corr_all_ranked.fig'])
saveas(fig, [output_folder '\AGSICI_muscles_corr_all_ranked.png'])

% update the counter
figure_counter = figure_counter + 1;  

clear data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 7) correlation per stimulation orientation
for p = 1:length(position)
    for c = 1:length(current)       
        % extract data
        data_corr_muscle = [];
        data_corr_N45 = []; 
        for i = 1:length(intensity)
            data_corr_muscle = [data_corr_muscle; squeeze(AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, :))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, :))];
        end
        data_corr = [data_corr_muscle data_corr_N45];
        clear i data_corr_muscle data_corr_N45
        
        % choose colours
        marker_col = [];
        for a = 1:size(data_corr, 1)
            marker_col = [marker_col; colours((p-1)*2 + c, :)];
        end
        clear a 
        
        % ----- linear correlation -----
        % prepare linear model: y ~ 1 + x
        data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {'muscular activity' 'N45'});

        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_corr(data_model, data_corr, marker_col, 'Pearson')
        title(sprintf('Linear correlation: %s - %s', position{p}, current{c}), 'FontWeight', 'bold', 'FontSize', 16)

        % save the figure       
        savefig([output_folder '\AGSICI_muscles_corr_' position{p} '-' current{c} '.fig'])
        saveas(fig, [output_folder '\AGSICI_muscles_corr_' position{p} '-' current{c} '.png'])

        % update the counter
        figure_counter = figure_counter + 1;  

        % ----- non-linear correlation -----
        % clculate correlation coeficient and p
        [cor_coef, cor_p] = corr(data_corr, 'Type', 'Spearman');

        % rank the data
        for a = 1:size(data_corr, 2)
            [temp, data_corr_ranked(:, a)]  = ismember(data_corr(:, a), unique(data_corr(:, a)));
        end
        clear a

        % prepare linear model: y ~ 1 + x
        data_model_ranked = fitlm(data_corr_ranked(:, 1), data_corr_ranked(:, 2), 'VarNames', {'muscular activity' 'N45'});

        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')
        title(sprintf('Non-linear correlation: %s - %s', position{p}, current{c}), 'FontWeight', 'bold', 'FontSize', 16)

        % save the figure       
        savefig([output_folder '\AGSICI_muscles_corr_ranked_' position{p} '-' current{c} '.fig'])
        saveas(fig, [output_folder '\AGSICI_muscles_corr_ranked_' position{p} '-' current{c} '.png'])

        % update the counter
        figure_counter = figure_counter + 1;  
    end
end
clear p c data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 8) correlation per stimulation intensity
% choose colours
marker_col = [];
for a = 1:4
    for b = 1:length(subject)
        marker_col = [marker_col; colours(a, :)];
    end
end
clear a b

% plot
for i = 1:length(intensity)         
    % extract data
    data_corr_muscle = [];
    data_corr_N45 = []; 
    for p = 1:length(position)
        for c = 1:length(current)   
            data_corr_muscle = [data_corr_muscle; squeeze(AGSICI_muscle_activity.muscles.GFP_TOI(p, c, i, :))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_TOI(p, c, i, :))];
        end
    end
    data_corr = [data_corr_muscle data_corr_N45];
    clear p c data_corr_muscle data_corr_N45
    
    % ----- linear correlation -----
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {'muscular activity' 'N45'});

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model, data_corr, marker_col, 'Pearson')
    title(sprintf('Linear correlation: %s %%rMT', intensity{i}(end-2:end)), 'FontWeight', 'bold', 'FontSize', 16)

    % save the figure       
    savefig([output_folder '\AGSICI_muscles_corr_int_' intensity{i}(end-2:end) '.fig'])
    saveas(fig, [output_folder '\AGSICI_muscles_corr_int_' intensity{i}(end-2:end) '.png'])

    % update the counter
    figure_counter = figure_counter + 1;  

    % ----- non-linear correlation -----
    % clculate correlation coeficient and p
    [cor_coef, cor_p] = corr(data_corr, 'Type', 'Spearman');

    % rank the data
    for a = 1:size(data_corr, 2)
        [temp, data_corr_ranked(:, a)]  = ismember(data_corr(:, a), unique(data_corr(:, a)));
    end
    clear a

    % prepare linear model: y ~ 1 + x
    data_model_ranked = fitlm(data_corr_ranked(:, 1), data_corr_ranked(:, 2), 'VarNames', {'muscular activity' 'N45'});

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')
    title(sprintf('Non-linear correlation, ranked: %s %%rMT', intensity{i}(end-2:end)), 'FontWeight', 'bold', 'FontSize', 16)

    % save the figure       
    savefig([output_folder '\AGSICI_muscles_corr_int_' intensity{i}(end-2:end) '_ranked.fig'])
    saveas(fig, [output_folder '\AGSICI_muscles_corr_int_' intensity{i}(end-2:end) '_ranked.png'])

    % update the counter
    figure_counter = figure_counter + 1;  
end

clear data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 9) tonic muscular activity
% ----- adjustable parameters -----
prefix_m_tonic = 'dc muscle ica visual crop but fft-notchfilt prefilt prea P1'; 
x_end_tonic = (-0.005 - header.xstart)/xstep;
% ----- adjustable parameters -----

% extract RMS of baseline signal
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % define subject
                if subject(s) < 10
                    subj = ['0' num2str(subject(s))];
                else
                    subj = num2str(subject(s));
                end
                
                % load and crop baseline data
                name = [prefix_m_tonic ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
                load(name)
                data = squeeze(data(:, 1:32, :, :, :, 1:x_end_tonic));
                
                % calculate RMS for each trial
                for t = 1:size(data, 1)
                    for e = 1:size(data, 2)
                        rms(t, e) = sqrt(mean(data(t, e, :).^2));
                    end
                end
                
                % average RMS across trials and electrodes
                AGSICI_muscle_activity.tonic(p, c, i, s) = mean(rms, 'all');                               
            end
        end
    end
end
clear p c i s subj data t e rms

% append new variables to the general MATLAB file
save(results_folder, 'AGSICI_muscle_activity', '-append');
%% functions
function plot_TEP(x, data_visual, varargin)
% check whether to plot labels 
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'peak_latency'));
    if ~isempty(a)
        latency = varargin{a + 1};
    else
        latency = false;
    end
else
    latency = false;
end

% check whether to plot legend
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'legend'));
    if ~isempty(a)
        legend_on = varargin{a + 1};
    else
        legend_on = {};
    end
else
    legend_on = {};
end

% check for colours
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'colours'));
    if ~isempty(a)
        col = varargin{a + 1};
    else
        col = [0.24 0.49 0.99; 0.72 0.27 1; 0.87 0.16 0.40; 0.99 0.18 0.18];
    end
else
    col = [0.24 0.49 0.99; 0.72 0.27 1; 0.87 0.16 0.40; 0.99 0.18 0.18];
end

% check for limits
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'limit'));
    if ~isempty(a)
        lim = varargin{a + 1};
    end
end

% set figure limits
hold on
if ~exist('lim')
    % identify dataset with largest + and - amplitudes
    for b = 1:size(data_visual, 1)
        data_max(b) = max(data_visual(b, :));
        data_min(b) = min(data_visual(b, :));
    end
    data_lim = [find(data_max == max(data_max)) find(data_min == min(data_min))];

    % find limits
    plot(x, data_visual(data_lim(1), :), x, data_visual(data_lim(2), :))
    yl = get(gca, 'ylim');
    lim = [yl(1), yl(2) + 0.15*(yl(2) - yl(1))];
end
ylim(lim)
xlim([x(1) x(end)])

% plot background objects
rectangle('Position', [0, lim(1), 0.01, lim(2) - lim(1)], ...
    'FaceColor', [0.75, 0.75, 0.75], 'EdgeColor', 'none')
line([x(1) x(end)], [0 0], 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1)

% plot data
for c = 1:size(data_visual, 1)
    P(c) = plot(x, data_visual(c, :), 'Color', col(c, :), 'LineWidth', 2);
end

% mark the TMS stimulus
line([0 0], lim, 'Color', [0 0 0], 'LineWidth', 2.5, 'LineStyle', '--')

% plot legend if required
if length(legend_on) > 0
    if length(legend_on) >= 4
        n_col = 2;
    else
        n_col = 1;
    end
    lgd = legend(P, legend_on, ...
        'Location', 'northeast', 'NumColumns', n_col);
    lgd.FontSize = 10;
end

% mark the peak latency, if required
if latency
    % identify the peak
    peak_pos = find(abs(data_visual) == max(abs(data_visual)));
    peak_x = x(peak_pos);
    peak_y = data_visual(peak_pos);

    % mark the peak
    plot(peak_x, peak_y, 'o', 'MarkerFaceColor', [0.888 0.196 0.028], 'MarkerSize', 12, 'MarkerEdgeColor', 'none')
    line([peak_x peak_x], [lim(1), peak_y], 'Color', [0.888 0.196 0.028], 'LineWidth', 2, 'LineStyle', ':')

    % add annotation
    text(0.025, lim(2) + 0.10*(lim(2) - lim(1)), sprintf('peak latency: %3.0f ms', ...
        round(peak_x *1000)), 'Color', [0.888 0.196 0.028], 'FontSize', 14)
end

% set other parameters
set(gca, 'fontsize', 12)
xlabel('time (s)')
ylabel('amplitude (\muV)')

hold off 
end
function plot_corr(data_model, data_corr, marker_col, corr_type)
% calculate correlation coefficient and p
[cor_coef, cor_p] = corr(data_corr, 'Type', corr_type);

% plot correlation
plot_cor = plotAdded(data_model);

% adjust parameters    
set(gca, 'FontSize', 14)
xlabel('muscular activity (mean GFP)'); ylabel('N45 (mean GFP)');
plot_cor(2).Color = [0 0 0]; plot_cor(2).LineWidth = 4; 
plot_cor(3).Color = [0 0 0]; plot_cor(3).LineWidth = 2;
legend off

% add annotations
if data_model.Coefficients.Estimate(2) > 0
    text_pos = [0.95 0.85 0.75];
else
    text_pos = [0.25 0.15 0.05];
end
T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(1, 2), cor_p(1, 2)), 'Units', 'Normalized');
set(T(1), 'fontsize', 14, 'fontangle', 'italic'); 
set(T(2), 'fontsize', 14); 
set(T(3), 'fontsize', 14, 'fontweight', 'bold'); 

% replot markers
for c = 1:size(data_corr, 1)
    scatter(data_corr(c, 1), data_corr(c, 2), 50, marker_col(c, :), 'filled');
    hold on
end
end