%% AG-SICI: TMS-EVOKED POTENTIALS - visualization
% Written by Dominika for AG-SICI project (2021 -22)
% 
% Colection of scripts to output the results of the study in plots 
% intended for publication (figures finalized in Illustrator)
% 
% Output:
%   --> figures are saved in a folder 'AG-SICI_P1_figures
%   --> variables are saved in a global MATLAB output file 'AG-SICI_P1.mat'
% 
% 1) PREPARE TEP DATA
%       - load individual data, remove possible outliers
%       - average across selected channels
%       - merge per categor & save for letswave
%       - compute average + SEM 
%       - load RAGU output
% 
% 2) TEP BUTTERFLY + GFP - individual conditions
%       - colourcoded
%       - optionally highlight selected channel
% 
% 3) GFP PEAKS - per position 
%       - calculate mean GFP across all intensities
%       - identify local maxima, extract latencies
%       - encode into 'AGSICI_TEP_default' structure
%       - plot GFP + butterfly + peak topographies
%
% 4) TANOVA P-VALUE + EV
%       - extract p-values and % of explained variance from RAGU output
%           --> for ALL comparisons
%       - plot, highlight significant intervals
% 
% 5) GFP P-VALUE + EV 
%       - extract p-values and % of explained variance from RAGU output
%           --> for selected comparison
%       - plot, highlight significant intervals
% 
% 6) GFP BOXPLOT
%       - compute individual GFP, average across chosen TOI
%       - plot a boxplot for chosen comparison
%           --> orientation / interaction
% 
% 7) GFP MEAN VALUES
%       - extract group GFP values from chosen TOI 
%       - encode to 'AGSICI_GFP_summary_WO' (WO = without outlier)
% 
% 8) GFP RADAR PLOT + CHANGE BOXPLOT
%       - extract group GFP values from chosen TOIs --> split P25 and N45
%       - plot radar plots per TOI - position on axes, group by intensity
%       - calculate GFP change between intensities
%       - plot boxplots per TOI
% 
% 9) MUSCLES OUTLIERS
%       - load unfiltered data
%       - compute individual GFP values of muscular artifact [3 10]ms
%       - plot boxplots, identify outliers
%       - plot boxplots WO
% 
% 10) MUSCLES MEAN VALUES - per position
%       - encode to 'AGSICI_muscles_summary_WO'
% 
% 11) MUSCLES LINEPLOT
%       - compute average muscle GFP + SEM
%       - plot all conditions in one figure
% 
% 12) MUSCLES SORT BY CONTRACTION SIZE - all conditions
%       - pool all subjects + all conditions 
%       - sort by GFP of muscular contraction
%       - split in 4 categories
%       - extract mean values, encode to 'AGSICI_muscles_summary'
%       - save corresponding TEP data for letswave
%           --> for topography extraction
%       - plot GFP timecourse
% 
% 13) MUSCLES ARTIFACT
%       - average across all subjects and conditions
%       - save corresponding TEP data for letswave
%           --> for topography extraction
%       - plot timecourse
% 
% 14) MUSCLES SORT BY CONTRACTION SIZE - all conditions
%       - without outliers
%       - same as (11) but separately for across and along datasets

%% parameters
clear all; clc

% ----------- dataset -----------
study = 'P1';
subject = [1, 3:18, 20, 21];
position = {'along_normal' 'along_reversed' 'across_normal' 'across_reversed'};
intensity = {'100' '120' '140'};
electrodes = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2',...
    'CP1','CP2','FC5','FC6','CP5','CP6','P5', 'P6', 'C1','C2'};
% --------------------------------
% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
addpath(folder_git)
figure_counter = 1;

% input/output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_figures = [folder_results '\AG-SICI_' study '_figures_SSP'];
output_file = [folder_results '\AG-SICI_' study '_SSP.mat'];

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for p = 1:length(position)
            for i = 1:length(intensity)
                colours((p-1)*length(intensity) + i, :) = uisetcolor; 
            end
        end
        save([folder_git '\AGSICI_' study '_colours.mat'], 'colours')
    case 'NO'
        if exist([folder_git '\AGSICI_' study '_colours.mat']) > 0
            load([folder_git '\AGSICI_' study '_colours.mat'])
        else
            disp('No colour scheme found in git directory!')    
        end
end
clear answer p i

%% 1) PREPARE TEP DATA
% ----- section input -----
outliers = [];
target = {'CP5' 'P3' 'P5' 'P7'};
% -------------------------
% load individual data, filter out outliers
load(output_file, 'data_individual')
load(output_file, 'header')
load(output_file, 'AGSICI_data')
idx = find(~ismember(1:size(AGSICI_data, 4), outliers));
counter = 1; 
for a = 1:size(AGSICI_data, 1)
    for b = 1:size(AGSICI_data, 2)
        for c = 1:size(AGSICI_data, 3)
            data_individual(counter, c, :, :, :) = AGSICI_data(a, b, c, idx, 1:length(electrodes), :);
        end
        counter = counter + 1;
    end
end
clear a b c AGSICI_data counter

% identify channels to pool
channels = [];
for t = 1:length(target)
    for e = 1:length(electrodes)
        if strcmp(electrodes{e}, target{t})
            channels(end+1) = e;
        end
    end
end

% create the target channel & update electrodes 
data_individual(:, :, :, size(data_individual, 4) + 1, :) = mean(data_individual(:, :, :, channels, :), 4);
electrodes{end + 1} = 'target';

% save data per category for letswave
load([folder_git '\header_example.mat'])
for p = 1:length(position)
    for i = 1:length(intensity)
        % subset data
        for s = 1:size(data_individual, 3)
            for e = 1:size(data_individual, 4)
                for t = 1:size(data_individual, 5)
                    data(s, e, 1, 1, 1, t) = data_individual(p, i, s, e, t);
                end
            end
        end
        
        % decide name
        name = sprintf('merged_subj AGSICI %s TEP %s %s', study, position{p}, intensity{i});
        
        % save for lw
        savelw(data, header, name)        
    end
end

% compute average and SEM
data_mean = squeeze(mean(data_individual, 3));
data_SEM =  squeeze(std(data_individual, 0, 3)/sqrt(size(data_individual, 3)));

% load RAGU data
load(output_file, 'AGSICI_microstates')
data_MS = AGSICI_microstates;
load([folder_results '\AG-SICI_P1_microstates\AG-SICI_microstates_SSP.mat'])
data_MS = rd;
clear AGSICI_microstates rd
save(output_file, 'data_MS', '-append');

% encode in a default structure
AGSICI_TEP_default = struct;
AGSICI_TEP_default.peak = {'P25' 'N45' 'N75' 'N100' 'P180'};
AGSICI_TEP_default.center = [];
AGSICI_TEP_default.span = [];
AGSICI_TEP_default.eoi = [];

% append to the outcome MATLAB file
save(output_file, 'AGSICI_TEP_default', '-append');
clear labels AG_peaks t e a b c counter AGSICI_muscle_activity outliers idx t e p i s t channels data header name

%% 2) TEP BUTTERFLY + GFP - individual conditions
% ----- section input -----
channel = 'target';
x_lim = [-50, 300];
x_step = 0.5;
y_lim = {[-5 5], [0 3]};
% -------------------------
% calculate mean GFP (exclude target channel)
for p = 1:length(position)
    for i = 1:length(intensity)
        data_gfp(p, i, :) = std(squeeze(data_mean(p, i, 1:end-1, :)), 1);  
    end 
end
clear p i

% prepare x axis
x = x_lim(1):x_step:x_lim(2);

% identify channel to highlight
for e = 1:length(electrodes)
    if strcmp(electrodes{e}, channel)
        channel_n = e;
    end
end

% loop through conditions to plot
for p = 1:length(position)
    for i = 1:length(intensity)
        % plot butterfly plot
        data_visual = squeeze(data_mean(p, i, :, :));
        fig = plot_TEP(figure_counter, x, data_visual, y_lim{1}, 'channel', channel_n, 'colour', colours((p-1)*length(intensity) + i, :));
        
        % name figure & save
        fig_name = sprintf('AGSICI_TEP_butterfly_%s_%s', position{p}, intensity{i});
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
        
        % plot GFP
        data_visual = squeeze(data_gfp(p, i, :));
        fig = plot_GFP(figure_counter, x, data_visual', y_lim{2}, 'colour', colours((p-1)*length(intensity) + i, :));
        
        % name & save
        fig_name = sprintf('AGSICI_TEP_GFP_%s_%s', position{p}, intensity{i});
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
    end
end
clear p i e fig fig_name data_visual channel x x_lim x_step y_lim channel_n

%% 3) GFP PEAKS - per position 
% ----- section input -----
labeled = 'off';
max_peaks = 6;
x_lim = [-50, 300];
x_step = 0.5;
y_lim = {[-5 5], [0 3]};
% -------------------------
% plot separately for different coil positions
for p = 1:length(position)
    % average across intensities & subjects
    data = squeeze(mean(data_individual(p, :, :, :, :), [2, 3]));
    
    % calculate GFP
    gfp = std(data, 1); 
    
    % ----- PEAK IDENTIFICATION -----
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % extract peak latencies
    h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
    AG_peaks = GFP_peaks(gfp, x_lim, x_step, labeled, 'max_peaks', max_peaks);    
    AGSICI_TEP_default.center{p} = AG_peaks;

    % choose data & header for topoplots 
    for e = 1:size(data, 1)
        for i = 1:size(data, 2)
            data_topoplot(1, e, 1, 1, 1, i) = data(e, i);
        end
    end
    load([folder_git '\header_example.lw6'], '-mat')
    header.datasize = size(data_topoplot);
    clear e i

    % add topoplots
    for t = 1:length(AG_peaks)
        % plot the topoplot
        h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);
        topo_plot(header, data_topoplot, AG_peaks(t)/1000, x_lim(1)/1000, [-2, 2])

        % shift down
        pos = get(h_axis(1 + t), 'Position');
        pos(2) = pos(2) - 0.05;
        set(h_axis(1 + t), 'Position', pos);

        % add timing
        text(-0.3, -0.8, sprintf('%1.0f ms', AG_peaks(t)), 'Color', [1 0 0], 'FontSize', 14)
    end
    clear t
    hold off

    % save figure
    figure_name = sprintf('AGSICI_TEP_GFP_%s_peaks', position{p});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'])

    % update figure counteer
    figure_counter = figure_counter + 1;
    
    % ----- BUTTERFLY + GFP PLOTS -----   
    % preapre x axis
    x = x_lim(1):x_step:x_lim(2);
    
    % plot butterfly plot
    fig = plot_TEP(figure_counter, x, data, y_lim{1}, 'colour', colours((p-1)*length(intensity) + 2, :));

    % name figure & save
    fig_name = sprintf('AGSICI_TEP_butterfly_%s', position{p});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
    figure_counter = figure_counter + 1;
    
    % plot GFP
    fig = plot_GFP(figure_counter, x, gfp, y_lim{2}, 'colour', colours((p-1)*length(intensity) + 2, :));

    % name & save
    fig_name = sprintf('AGSICI_TEP_GFP_%s', position{p});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
    figure_counter = figure_counter + 1;
end
clear labeled max_peaks p data gfp fig fig_name h_axis data_topoplot header pos figure_name x x_lim x_step y_lim

%% 4) TANOVA P-VALUE + EV
% ----- section input -----
x_lim = [10, 300];
x_step = 0.5;
% -------------------------
% define explored comparisons
comparison = {data_MS.strF1, data_MS.strF2, [data_MS.strF1 '_' data_MS.strF2]};

% preapre x axis
x = x_lim(1):x_step:x_lim(2);

% plot TANOVA results for all comparisons
for c = 1:numel(comparison)
    % ----- TANOVA p value -----
    % load data
    x_start = (x_lim(1) + 50)/x_step;
    x_end = (x_lim(2) + 50)/x_step;
    data_visual = squeeze(data_MS.PTanova(1, 1 + c, x_start:x_end, 1));      

    % plot the p-value timecourse
    fig = plot_p(x, data_visual, figure_counter);                   

    % name and save figure
    fig_name = ['TANOVA_' comparison{c}];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;       
    
    % ----- explained variance -----
    % load data
    data_visual = squeeze(data_MS. TExpVar{1}(1, 1 + c, x_start:x_end));    

    % plot EV timecourse
    fig = plot_EV(x, data_visual, figure_counter);

    % name and save figure
    fig_name = ['TANOVA_EV_' comparison{c}];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;       
end
clear comparison c data_visual fig fig_name x x_start x_end x_lim x_step

%% 5) GFP P-VALUE + EV
% ----- section input -----
factor = 'interaction';         % orientation - intensity - interaction
x_lim = [10, 300];
x_step = 0.5;
% -------------------------
% identify factor
switch factor
    case 'orientation'
        factor_n = 2;
    case 'intensity'
        factor_n = 3;
    case 'interaction'
        factor_n = 4;
end

% preapre x axis
x = x_lim(1):x_step:x_lim(2);

% ----- p value -----
% load data
x_start = (x_lim(1) + 50)/x_step;
x_end = (x_lim(2) + 50)/x_step;
data_visual = squeeze(data_MS.GFPPTanova(1, factor_n, x_start:x_end, 1))';      

% plot the p-value timecourse
fig = plot_p(x, data_visual, figure_counter);                   

% name and save figure
fig_name = ['GFP_' factor];
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

% update figure counter
figure_counter = figure_counter + 1 ;       

% ----- explained variance -----
% load data
data_visual = squeeze(data_MS.GFPExpVar{1}(1, factor_n, x_start:x_end));    

% plot EV timecourse
fig = plot_EV(x, data_visual, figure_counter);

% name and save figure
fig_name = ['GFP_EV_' factor];
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

% update figure counter
figure_counter = figure_counter + 1 ;

clear factor factor_n data_visual fig x_lim x_step x x_start x_end fig fig_name

%% 6) GFP BOXPLOT
% ----- section input -----
factor = 'interaction';         % orientation - interaction
toi = [43, 56];
x_step = 0.5;
y_lim = [0, 4.5];
% -------------------------
% determine y info for plotting
y_info = {'GFP (\muV)', y_lim};

% switch between factors
switch factor
    case 'orientation'
        % define x info for plotting
        x_info = [factor position];

        % calculate individual GFP per condition (exclude target channel)
        data_visual = [];
        for s = 1:size(data_individual, 3)            
            for p = 1:length(position) 
                data_visual(s, p, :) = std(squeeze(mean(data_individual(p, :, s, 1:end-1, :), 2)), 1);  
            end
        end

        % average across toi
        x_start = (toi(1) + 50)/x_step;
        x_end = (toi(2) + 50)/x_step;
        data_visual = mean(data_visual(:, :, x_start:x_end), 3);
        
        % plot boxplot
        fig = plot_scatter(figure_counter, data_visual, x_info, y_info, 'colour', colours([2, 5, 8, 11], :));
        
        % name and save figure
        fig_name = sprintf('GFP_boxplot_%s', factor);
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;
        
    case 'interaction'              % plots a separate plot per intensity
        % loop through positions
         for i = 1:length(intensity)         
            % define x info for plotting
            x_info = ['orientation' position];

            % calculate individual GFP (exclude target channel)
            data_visual = [];
            for s = 1:size(data_individual, 3)            
                for p = 1:length(position) 
                    data_visual(s, p, :) = std(squeeze(data_individual(p, i, s, 1:end-1, :)), 1);  
                end
            end

            % average across toi
            x_start = (toi(1) + 50)/x_step;
            x_end = (toi(2) + 50)/x_step;
            data_visual = mean(data_visual(:, :, x_start:x_end), 3);

            % plot boxplot
            fig = plot_scatter(figure_counter, data_visual, x_info, y_info, 'colour', colours((i-1) + [1, 4, 7, 10], :));

            % name and save figure
            fig_name = sprintf('GFP_boxplot_%s_%s_%s', factor, intensity{i}, 'p0.01');
            savefig([folder_figures '\' fig_name '.fig'])
            saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

            % update figure counter
            figure_counter = figure_counter + 1 ;
         end
end

clear factor toi x_step y_lim x_start x_end x_info y_info data_visual fig fig_name p s i

%% 7) GFP MEAN VALUES
% ----- section input -----
toi = [43, 56];
x_step = 0.5;
% -------------------------
% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% create output structure
AGSICI_GFP_summary = struct;
for p = 1:length(position)
    for i = 1:length(intensity)
        % calculate individual GFP per condition (exclude target channel)
        data = [];
        for s = 1:size(data_individual, 3)            
            data(s, :) = std(squeeze(data_individual(p, i, s, 1:end-1, :)), 1);  
        end
        data = mean(data(:, x_start:x_end), 2);

        % group median
        statement = ['AGSICI_GFP_summary.' position{p} '_' intensity{i} '.median = median(data);'];
        eval(statement)

        % group median
        statement = ['AGSICI_GFP_summary.' position{p} '_' intensity{i} '.mean = mean(data);'];
        eval(statement)

        % group median
        statement = ['AGSICI_GFP_summary.' position{p} '_' intensity{i} '.SD = std(data);'];
        eval(statement)

        % group median
        statement = ['AGSICI_GFP_summary.' position{p} '_' intensity{i} '.SEM = std(data)/sqrt(length(data));'];
        eval(statement)
    end
end

% append to general matlab file
save(output_file, 'AGSICI_GFP_summary', '-append');
clear toi x_step x_start x_end p i s data statement

%% 8) GFP RADAR PLOT + CHANGE BOXPLOT
% ----- section input -----
toi = {[43, 56]};
x_step = 0.5;
y_lim_rad = [0.8, 2.3]; 
y_lim_box = [-1.2, 2.2];
% -------------------------
% plot radar plot
for t = 1:length(toi)
    % calculate individual GFP (exclude target channel)
    data = [];
    for i = 1:length(intensity)          
        for p = 1:length(position) 
            for s = 1:size(data_individual, 3)   
                data(i, p, s, :) = std(squeeze(data_individual(p, i, s, 1:end-1, :)), 1);  
            end
        end
    end

    % average across subjects and toi
    x_start = (toi{t}(1) + 50)/x_step;
    x_end = (toi{t}(2) + 50)/x_step;
    data_visual = mean(data(:, :, :, x_start:x_end), [3, 4]);

    % reorder positions
    data_visual = data_visual(:, [3, 2, 4, 1]);

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % plot the radar plot
    spider_plot_R2019b(data_visual, 'AxesLabels', position([3, 2, 4, 1]), ...
        'AxesLimits', [repelem(y_lim_rad(1), size(data_visual, 2)); repelem(y_lim_rad(2), size(data_visual, 2))],...
        'AxesAngular', 'off')

    % add legend
    legend('100% rMT', '120% rMT', '140% rMT', 'Location', 'southoutside');

    % name and save figure
    fig_name = sprintf('GFP_radarplot_%s', AGSICI_TEP_default.peak{t});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;
end

% determine info for boxplot
x_info = ['orientation' position];
y_info = {'\Delta GFP (\muV)', y_lim_box};

% plot change boxplot
for t = 1:length(toi)
    % calculate individual GFP (exclude target channel)
    data = [];
    for i = 1:length(intensity)          
        for p = 1:length(position) 
            for s = 1:size(data_individual, 3)   
                data(i, p, s, :) = std(squeeze(data_individual(p, i, s, 1:end-1, :)), 1);  
            end
        end
    end
    
    % average across toi
    x_start = (toi{t}(1) + 50)/x_step;
    x_end = (toi{t}(2) + 50)/x_step;
    data = mean(data(:, :, :, x_start:x_end), 4);
    
    % for each intensity change
    for a = 1:length(intensity) - 1
        % calculate the change 
        data_visual = [];
        for p = 1:length(position) 
            for s = 1:size(data_individual, 3)
                data_visual(s, p) = data(a+1, p, s) - data(a, p, s); 
            end
        end
        
        % plot a boxplot
        fig = plot_scatter(figure_counter, data_visual, x_info, y_info, 'colour', colours([2, 5, 8, 11], :), 'zero_line', 'on');
        
        % name and save figure
        fig_name = sprintf('GFP_boxplot_change_%s_%s', AGSICI_TEP_default.peak{t}, [num2str(intensity{a+1}) '-' num2str(intensity{a})]);
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;
    end
end
clear toi x_step y_lim t data data_visual i p s x_start x_end fig a 

%% 9) MUSCLES OUTLIERS
% ----- section input -----
toi = [3, 10];
x_step = 0.5;
% -------------------------
% load individual data
load(output_file, 'AGSICI_muscle_activity')
counter = 1; 
for a = 1:size(AGSICI_muscle_activity.contraction.data, 1)
    for b = 1:size(AGSICI_muscle_activity.contraction.data, 2)
        for c = 1:size(AGSICI_muscle_activity.contraction.data, 3)
            data_muscles(counter, c, :, :, :) = AGSICI_muscle_activity.contraction.data(a, b, c, :, :, :);   
            data_muscles_gfp(counter, c, :, :) = AGSICI_muscle_activity.contraction.GFP(a, b, c, :, :);
        end
        counter = counter + 1;
    end
end
clear a b c AGSICI_muscle_activity

%% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% average individual gfp across TOI timepoints and all intensities
for s = 1:size(data_muscles_gfp, 3)
    for p = 1:size(data_muscles_gfp, 1)
        data_visual(s, p) = mean(data_muscles_gfp(p, :, s, x_start:x_end), [2, 4]);
    end
end
clear s p

% plot GFP per condition, identify outliers
[fig, outliers] = plot_box(figure_counter, data_visual, colours([3, 6, 9, 12], :));
outliers = sort(outliers(1, :));

% name and save figure
fig_name = 'MUSCLES_outliers';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

% update the counter
figure_counter = figure_counter + 1;  

% boxplot without outliers (average across intensities)
idx = find(~ismember(1:size(data_muscles_gfp, 3), outliers));
fig = plot_box(figure_counter, data_visual(idx, :), colours([2, 5, 8, 11], :));

% name and save figure
fig_name = 'MUSCLES_WO';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

% update the counter
figure_counter = figure_counter + 1; 

clear toi x_step counter x_start x_end data_visual fig fig_name idx

%% 10) MUSCLES MEAN VALUES - per position
% ----- section input -----
toi = [3, 10];
x_step = 0.5;
% -------------------------
% identify outliers
idx = find(~ismember(1:size(data_muscles_gfp, 3), outliers));

% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% create output structure
for p = 1:length(position)
    % subset data, average across intensities
    data = squeeze(mean(data_muscles_gfp(p, :, idx, x_start:x_end), [2, 4]));
    
    % group median
    statement = ['AGSICI_muscles_summary_WO.' position{p} '.median = median(data);'];
    eval(statement)
    
    % group median
    statement = ['AGSICI_muscles_summary_WO.' position{p} '.mean = mean(data);'];
    eval(statement)
    
    % group median
    statement = ['AGSICI_muscles_summary_WO.' position{p} '.SD = std(data);'];
    eval(statement)
    
    % group median
    statement = ['AGSICI_muscles_summary_WO.' position{p} '.SEM = std(data)/sqrt(length(data));'];
    eval(statement)
end

% append to general output file
save(output_file, 'AGSICI_muscles_summary_WO', '-append');
clear toi x_step idx x_start x_end p data statement

%% 11) MUSCLES LINEPLOT
% ----- section input -----
toi = [3 10];
x_step = 0.5;
% -------------------------
% identify outliers
idx = find(~ismember(1:size(data_muscles_gfp, 3), outliers));

% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% launch the figure
fig = figure(figure_counter); 
hold on

% plot the data
for p = 1:length(position)
    % calculate mean and SEM
    for i = 1:length(intensity)            
        y(i) = mean(data_muscles_gfp(p, i, idx, x_start:x_end), [3, 4]);
        SEM(i) = std(squeeze(mean(data_muscles_gfp(p, i, idx, x_start:x_end), 4))) / sqrt(length(idx));
    end

    % plot
    perr(p) = errorbar(1:length(y), y, SEM);

    % adjust parameters
    perr(p).Color = colours((p-1)*length(intensity) + 2, :);
    perr(p).LineWidth = 1.5;
    perr(p).Marker = 'o';
    perr(p).MarkerFaceColor = colours((p-1)*length(intensity) + 2, :);
    perr(p).MarkerSize = 10;
end

% add features, adjust parameters
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
xlabel('stimulation intensity (% rMT)'); 
ylabel('GFP (\muV \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
      
% name and save figure
fig_name = 'MUSCLES_lineplot_all_WO';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

% update the counter
figure_counter = figure_counter + 1;  

clear p i y SEM toi x_step x_start x_end fig perr fig_name idx

%% 12) MUSCLES SORT BY CONTRACTION SIZE - all conditions
% ----- section input -----
toi = [3 10];
x_lim = [-5, 30];
x_step = 0.5;
y_lim = [-10, 600];
categories = 4;
% -------------------------
% identify outliers
idx = find(~ismember(1:size(data_muscles_gfp, 3), outliers));

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for a = 1:2
            for c = 1:categories
                colours_muscles((a-1)*categories + c, :) = uisetcolor; 
            end
        end
        save([folder_git '\AGSICI_' study 'colours_muscles.mat'], 'colours_muscles')
    case 'NO'
        if exist([folder_git '\AGSICI_' study '_colours_muscles.mat']) > 0
            load([folder_git '\AGSICI_' study '_colours_muscles.mat'])
        else
            disp('No colour scheme found in git directory!')    
        end
end
clear answer a c

% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% create a table
muscles_table = table;
counter = 1;
for p = 1:length(position)
    for i = 1:length(intensity)
        for s = 1:size(data_muscles_gfp, 3) 
            % fill in the line
            muscles_table.subject(counter) = subject(s);
            muscles_table.subject_n(counter) = s;
            muscles_table.position(counter) = p;
            muscles_table.intensity(counter) = i;
            muscles_table.contraction(counter) = mean(data_muscles_gfp(p, i, s, x_start:x_end));
            
            % update the counter
            counter = counter + 1;
        end
    end
end
AGSICI_muscles_summary.data = muscles_table;

% sort by contraction size
muscles_table = sortrows(muscles_table, 5); 
for c = 1:categories
    muscles_table.category((c-1)*(height(muscles_table)/categories) + [1:height(muscles_table)/categories]) = c;
end

% loop through categories
for c = 1:categories
    % subset the data
    muscles_subset = muscles_table(muscles_table.category == c, 2:5);
    
    % count position ratio
    for p = 1:length(position)
        AGSICI_muscles_summary.all(c).position(p) = height(muscles_subset(muscles_subset.position == p, 1));
    end
    
    % count intensity ratio
    for i = 1:length(intensity)
        AGSICI_muscles_summary.all(c).intensity(i) = height(muscles_subset(muscles_subset.intensity == i, 1));
    end
    
    % extract mean values
    AGSICI_muscles_summary.all(c).mean = mean(muscles_subset.contraction);
    AGSICI_muscles_summary.all(c).SD = std(muscles_subset.contraction);
    AGSICI_muscles_summary.all(c).SEM = std(muscles_subset.contraction)/sqrt(size(data_muscles, 3));
    
    % save muscle data for letswave
    clear data header
    for a = 1:height(muscles_subset)
        data(a, :, 1, 1, 1, :) = data_muscles(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :, :);
    end
    load([folder_git '\header_example.mat'])
    name = sprintf('merged_subj AGSICI %s muscles all cat%d', study, c);
    savelw(data, header, name) 
    
    % save TEP data for letswave
    clear data header
    for a = 1:height(muscles_subset)
        data(a, :, 1, 1, 1, :) = data_individual(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :, :);
    end
    load([folder_git '\header_example.mat'])
    name = sprintf('merged_subj AGSICI %s muscles TEP all cat%d', study, c);
    savelw(data, header, name) 
    
    % prepare GFP for visualization
    for a = 1:height(muscles_subset)
        data_visual(a, :) = data_muscles_gfp(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :);
    end
    x = x_lim(1):x_step:x_lim(2);
    y = mean(data_visual(:, (x_lim(1) + 50)/x_step:(x_lim(2) + 50)/x_step));
    SEM = std(data_visual(:, (x_lim(1) + 50)/x_step:(x_lim(2) + 50)/x_step))/sqrt(size(data_visual, 1));
    
    % plot GFP    
    fig = plot_muscle(figure_counter, x, y, SEM, {'time (ms)' x_lim}, {'GFP (\muV)' y_lim}, 'toi', toi, 'colour', colours_muscles(c, :))
    
    % name and save figure
    fig_name = sprintf('MUSCLES_contraction_all_cat%d', c);
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

    % update the counter
    figure_counter = figure_counter + 1;    
end

% append to general output file
save(output_file, 'AGSICI_muscles_summary', '-append');
clear toi x_lim y_lim x_step x_start x_end p i s counter categories c a data header name data_visual x y SEM fig fig_name muscles_subset

%% 13) MUSCLES ARTIFACT
% ----- section input -----
channel = 'CP5';
x_lim = [-5, 30];
x_step = 0.5;
y_lim = [-50, 150];
interpolation = [-2 2];
% -------------------------
% avereage data across categories
data_visual = squeeze(mean(data_muscles(:, :, :, :, (x_lim(1) + 50)/x_step:(x_lim(2) + 50)/x_step), [1 2 3]));

% preapre x axis
x = x_lim(1):x_step:x_lim(2);

% identify channel to highlight
for e = 1:length(electrodes)
    if strcmp(electrodes{e}, channel)
        channel_n = e;
    end
end

% plot butterfly plot + highlight chosen channel
fig = plot_TEP(figure_counter, x, data_visual, y_lim, 'channel', channel_n, 'interpolation', interpolation);
        
% name figure & save
fig_name = sprintf('MUSCLES_butterfly_%s', channel);
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
figure_counter = figure_counter + 1;
        
clear channel interpolation x_lim x_step y_lim data_visual x e channel_n fig fig_name

%% 14) MUSCLES SORT BY CONTRACTION SIZE - per position
% ----- section input -----
orientation = {'along' 'across'};
toi = [3 10];
x_lim = [-5, 30];
x_step = 0.5;
y_lim = [-10, 350];
categories = 4;
% -------------------------
% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for a = 1:2
            for c = 1:categories
                colours_muscles((a-1)*categories + c, :) = uisetcolor; 
            end
        end
        save([folder_git '\AGSICI_' study 'colours_muscles.mat'], 'colours_muscles')
    case 'NO'
        if exist([folder_git '\AGSICI_' study '_colours_muscles.mat']) > 0
            load([folder_git '\AGSICI_' study '_colours_muscles.mat'])
        else
            disp('No colour scheme found in git directory!')    
        end
end
clear answer a c

% load the data table 
if exist('muscles_table')~=1
    load(output_file, 'AGSICI_muscles_summary');
end

% deal with potential outliers
muscles_table = AGSICI_muscles_summary.data(~ismember(AGSICI_muscles_summary.data.subject_n, outliers), :);
AGSICI_muscles_summary_WO.data = muscles_table;

% loop through orientations
for o = 1:length(orientation)    
    % subset and sort by contraction size
    muscles_table = AGSICI_muscles_summary_WO.data(AGSICI_muscles_summary_WO.data.position == (o - 1)*2 + 1 | ...
        AGSICI_muscles_summary_WO.data.position == (o - 1)*2 + 2, :);
    muscles_table = sortrows(muscles_table, 5); 
    for c = 1:categories
        muscles_table.category((c-1)*(height(muscles_table)/categories) + [1:height(muscles_table)/categories]) = c;
    end

    % loop through categories
    for c = 1:categories
        % subset the data
        muscles_subset = muscles_table(muscles_table.category == c, 2:5);

        % count position ratio
        for p = 1:length(position)
            statement = ['AGSICI_muscles_summary_WO.' orientation{o} '(c).position(p) = height(muscles_subset(muscles_subset.position == p, 1));'];
            eval(statement)
        end

        % count intensity ratio
        for i = 1:length(intensity)
            statement = ['AGSICI_muscles_summary_WO.' orientation{o} '(c).intensity(i) = height(muscles_subset(muscles_subset.intensity == i, 1));'];
            eval(statement)
        end

        % extract mean values
        statement = ['AGSICI_muscles_summary_WO.' orientation{o} '(c).mean = mean(muscles_subset.contraction);'];
        eval(statement)
        statement = ['AGSICI_muscles_summary_WO.' orientation{o} '(c).SD = std(muscles_subset.contraction);'];
        eval(statement)
        statement = ['AGSICI_muscles_summary_WO.' orientation{o} '(c).SEM = std(muscles_subset.contraction)/sqrt(size(data_muscles, 3) - length(outliers));'];
        eval(statement)
        
        % save muscle data for letswave
        clear data header
        for a = 1:height(muscles_subset)
            data(a, :, 1, 1, 1, :) = data_muscles(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :, :);
        end
        load([folder_git '\header_example.mat'])
        name = sprintf('merged_subj AGSICI %s muscles WO %s cat%d', study, orientation{o}, c);
        savelw(data, header, name) 

        % save TEP data for letswave
        clear data header
        for a = 1:height(muscles_subset)
            data(a, :, 1, 1, 1, :) = data_individual(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :, :);
        end
        load([folder_git '\header_example.mat'])
        name = sprintf('merged_subj AGSICI %s muscles TEP WO %s cat%d', study, orientation{o}, c);
        savelw(data, header, name) 

        % prepare GFP for visualization
        for a = 1:height(muscles_subset)
            data_visual(a, :) = data_muscles_gfp(muscles_subset.position(a), muscles_subset.intensity(a), muscles_subset.subject_n(a), :);
        end
        x = x_lim(1):x_step:x_lim(2);
        y = mean(data_visual(:, (x_lim(1) + 50)/x_step:(x_lim(2) + 50)/x_step));
        SEM = std(data_visual(:, (x_lim(1) + 50)/x_step:(x_lim(2) + 50)/x_step))/sqrt(size(data_visual, 1));

        % plot GFP    
        fig = plot_muscle(figure_counter, x, y, SEM, {'time (ms)' x_lim}, {'GFP (\muV)' y_lim}, ...
            'toi', toi, 'colour', colours_muscles((o-1)*categories + c, :))

        % name and save figure
        fig_name = sprintf('MUSCLES_contraction_WO_%s_cat%d', orientation{o}, c);
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')

        % update the counter
        figure_counter = figure_counter + 1;    
    end
end

% append to the output file
save(output_file, 'AGSICI_muscles_summary_WO', '-append');
clear orientation toi x_lim y_lim x_step x_start x_end p i s o counter categories c a data header name ...
    data_visual x y SEM fig fig_name muscles_subset statement

%% functions
function fig = plot_TEP(figure_counter, x, data_visual, y_lim, varargin)
    % deal with optional arguments
    if ~isempty(varargin)
        % channel to highlight
        a = find(strcmpi(varargin, 'channel'));
        if ~isempty(a)
            channel_n = varargin{a + 1};
        else
            channel_n = [];
        end
        
        % line colour
        b = find(strcmpi(varargin, 'colour'));
        if ~isempty(b)
            colour = varargin{b + 1};
        else
            colour = [0.65 0.65 0.65];
        end
        
        % interpolated interval
        c = find(strcmpi(varargin, 'interpolation'));
        if ~isempty(c)
            int = varargin{c + 1};
        else
            int = [-5 10];
        end
    end    
    
    % launch the figure
    fig = figure(figure_counter);
    set(gcf, 'units','centimeters','position',[10 10 20 10], 'color', 'w');
    hold on
    
    % determine axis properties
    xl = [x(1) - 25, x(end) + 25];
    xlim(xl);
    xlabel('time (ms)');  
    ylim(y_lim);
    ylabel('amplitude (\muV)')

    % loop through channels to plot
    for a = 1:size(data_visual, 1)     
        P(a) = plot(x, data_visual(a, :), 'Color', colour, 'LineWidth', 1);
    end

    % highlight channel
    if ~isempty(channel_n)
        P(end+1) =  plot(x, data_visual(channel_n, :), 'Color', [0 0 0], 'LineWidth', 3);
    end

    % shade interpolated interval 
    rectangle('Position', [int(1), y_lim(1), int(2) - int(1), y_lim(2) - y_lim(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')

    % TMS stimulus
    line([0, 0], y_lim, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % other parameters
    xlim([x(1) - length(x)*(x(2) - x(1))*0.05, x(end) + length(x)*(x(2) - x(1))*0.05])
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
end
function fig = plot_GFP(figure_counter, x, data_visual, y_lim, varargin)
    % deal with optional arguments
    if ~isempty(varargin)       
        % fill colour
        a = find(strcmpi(varargin, 'colour'));
        if ~isempty(a)
            colour = varargin{a + 1};
        else
            colour = [0.7 0.7 0.7];
        end
    end    
    
    % launch the figure
    fig = figure(figure_counter);
    set(gcf, 'units','centimeters','position',[10 10 20 7], 'color', 'w');
    hold on
    
    % determine axis properties
    xl = [x(1) - 25, x(end) + 25];
    xlim(xl);
    xlabel('time (ms)');  
    ylim(y_lim);
    ylabel('amplitude (\muV)')
    
    % shade GFP
    F = fill([x fliplr(x)],[data_visual zeros(1, length(x))], colour, 'linestyle', 'none');

    % plot GFP line
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);

    % shade interpolated interval 
    rectangle('Position', [-5, y_lim(1), 15, y_lim(2) - y_lim(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')

    % TMS stimulus
    line([0, 0], y_lim, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % other parameters
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
end
function savelw(data, header, name)
    % modify header
    header.name = name;
    header.datasize = size(data);
    header.xstart = -0.05;
    for e = 1:size(data, 1)
        header.events(e).code = 'AG_TEP';
        header.events(e).latency = 0;
        header.events(e).epoch = e;
    end
    
    % save 
    save([name '.mat'], 'data')
    save([name '.lw6'], 'header')
end
function peak_x = GFP_peaks(y, time_window, xstep, labeled, varargin)
    % check whether to plot labels (default)
    if ~isempty(varargin)
        a = find(strcmpi(varargin, 'max_peaks'));
        if ~isempty(a)
            max_peaks = varargin{a + 1};
        end
    end
    
    % preapre x axis
    x = time_window(1):xstep:time_window(2);

    % launch the figure  
    plot(x, y)
    yl = get(gca, 'ylim');
    cla

    % plot interpolated part
    hold on
    xlim(time_window)
    rectangle('Position', [-5, yl(1), 15, yl(2) - yl(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')

    % plot data, mark TMS stimulus
    plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 3)

    % find peaks 
    [pks, locs] = findpeaks(y, 'MinPeakDistance', 5, 'MinPeakProminence', 0.005);
    for a = 1:length(locs)
        if time_window(1) + locs(a)*xstep <= 15
            idx(a) = false;
        elseif time_window(1) + locs(a)*xstep > 220
            idx(a) = false;        
        else
            idx(a) = true;
        end
    end
    pks = pks(idx); locs = locs(idx);
    if length(pks) > max_peaks
        pks = pks(1:max_peaks); 
        locs = locs(1:max_peaks);
    end

    % calculate peak coordinations
    for a = 1:length(locs)
        peak_x(a) = time_window(1) + locs(a)*xstep;
        peak_y(a) = pks(a);
    end
    peak_y = double(peak_y);

    % plot peaks
    for a = 1:length(locs)
        plot(peak_x(a), peak_y(a), ...
            'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
        line([peak_x(a), peak_x(a)], [yl(1), peak_y(a)], 'LineStyle', ':', 'Color', [1, 0, 0], 'LineWidth', 1.5)

        % label with latency (ms)
        if strcmp(labeled, 'on') 
            text(peak_x(a), peak_y(a) + 0.15, sprintf('%1.0fms', peak_x(a)*1000), 'Color', [1 0 0], 'FontSize', 14)
        end
    end

    % add parameters
    set(gca, 'Layer', 'Top')
    set(gca, 'fontsize', 14)
    ylim(yl)
    xlabel('time (ms)')
    ylabel('GFP (\muV)')
end
function topo_plot(header, data, x_pos, x_start, map_lims)
    varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

    % fetch data to display
    x_visual = ceil((x_pos - x_start)/header.xstep);
    vector = data(1, :, 1, 1, 1, x_visual);

    %fetch chanlocs
    chanlocs = header.chanlocs;

    %parse data and chanlocs 
    i=1;
    for chanpos=1:size(chanlocs,2);
        vector2(i)=double(vector(chanpos));
        chanlocs2(i)=chanlocs(chanpos);
        i=i+1;
    end;

    topoplot(vector2,chanlocs2,varargin{:});
    set(gcf,'color',[1 1 1]);
end
function fig = plot_p(x, data_visual, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % define intervals of significance
    signif_05 = data_visual < 0.05;
    signif_01 = data_visual < 0.01;
    
    % shade intervals of significance
    I(1) = area(x, signif_05);
    I(1).FaceColor = [1 0.73 0.73];
    I(1).EdgeColor = 'none';
    I(2) = area(x, signif_01);
    I(2).FaceColor = [1 0.44 0.44];
    I(2).EdgeColor = 'none';
    
    % plot p
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);
    
    % other parameters  
    xlabel('time (ms)')
    ylabel('p-value')
    set(gca, 'FontSize', 18)
    xlim([x(1), x(end)])    
    set(gca, 'layer', 'top');
    
    % change figure size
    fig.Position = [500 500 750 300];
    
    hold off
end
function fig = plot_EV(x, data_visual, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % shade area of explained variance
    A = area(x, data_visual*100);
    A.FaceColor = [0.65 0.65 0.65];
    A.EdgeColor = 'none';

    % other parameters    
    xlabel('time (ms)')
    ylabel('% explained variance')
    set(gca, 'FontSize', 18)
    xlim([x(1), x(end)])   
    set(gca, 'layer', 'top');

    % change figure size
    fig.Position = [500 500 750 300];
    
    hold off
end
function fig = plot_scatter(figure_counter, data_visual, x_info, y_info, varargin)
    % deal with optional arguments
    if ~isempty(varargin)        
        % line colour
        a = find(strcmpi(varargin, 'colour'));
        if ~isempty(a)
            colour = varargin{a + 1};
        else
            for b = 1:size(data_visual, 2)
                colour(b, :) = [0.7 0.7 0.7];
            end
        end
        
        % line colour
        c = find(strcmpi(varargin, 'zero_line'));
        if ~isempty(c)
            zero_line = varargin{c + 1};
        else
            zero_line = 'off';
        end
    end 
    
    % prepare data
    data_group = [];
    data_y = [];
    for d = 1:size(data_visual, 2)
        data_group(end + 1:end + size(data_visual, 1), 1) = d;
        data_y = [data_y; data_visual(:, d)];
    end
    
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % plot the markers
    beeswarm(data_group, data_y)
    
    % add boxplot
    for a = 1:size(data_visual, 2)
        boxchart(a * ones(size(data_visual, 1), 1), data_visual(:, a), 'BoxFaceColor', colour(a, :))
    end
    
    % add zero line if required
    if strcmp(zero_line, 'on')
        line([0.25, size(data_visual, 2) + 0.25], [0, 0], 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle', ':')
    end
    
    % font
    set(gca, 'Fontsize', 18)
    
    % x label
    xlabel(x_info{1})
    set(gca, 'xtick', 1:size(data_visual, 2), 'xticklabel', x_info(2:end))
    
    % y label
    ylabel(y_info{1})
    ylim(y_info{2})    
end
function beeswarm(x,y,varargin)
% function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',20/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',1)
addOptional(p,'MarkerEdgeColor', 'black')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .9; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    cla
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold on
    xlim(xl)
    ylim(yl)
end

% unsort so that output matches the original y data
x(sid)=x;
end
function fig = plot_muscle(figure_counter, x, y, SEM, x_info, y_info, varargin)
    % deal with optional arguments
    if ~isempty(varargin)        
        % colour
        a = find(strcmpi(varargin, 'colour'));
        if ~isempty(a)
            colour = varargin{a + 1};
        else
            colour = [0 0 0];
        end
        
        % toi
        b = find(strcmpi(varargin, 'toi'));
        if ~isempty(b)
            toi = varargin{b + 1};
        else
            toi = [];
        end
    end 
    
    % launch the figure
    fig = figure(figure_counter);
    set(gcf, 'units','centimeters','position',[10 10 14 10], 'color', 'w');
    hold on
    
    % shade TOI
    if ~isempty(toi)
        rectangle('Position', [toi(1), y_info{2}(1), toi(2), y_info{2}(2) - y_info{2}(1)], 'FaceColor', [0.98 0.98 0.43], 'EdgeColor', 'none')
    end
    
    % shade SEM
    F = fill([x fliplr(x)],[y + SEM fliplr(y - SEM)], colour, 'FaceAlpha', 0.25, 'linestyle', 'none');
    
    % plot data
    P = plot(x, y, 'Color', colour, 'LineWidth', 2);
    
    % mark TMS stimulus
    line([0, 0], y_info{2}, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % x params
    xlim([x_info{2}(1) - 2, x_info{2}(2) + 2])
    xlabel(x_info{1})
    
    % y params
    ylim(y_info{2})
    ylabel(y_info{1})
    
    % set other parameters
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')  
end
function [fig, varargout] = plot_box(figure_counter, data_visual, colour)
    % launch the figure
    fig = figure(figure_counter); 
    ax = gca;       
    hold on

    % plot the data
    boxplot(data_visual, 'color', colour)

    % axis labels
    ax.XTickLabel = '';
    label_array = {'along' 'along' 'across' 'across'; 'normal' 'reversed' 'normal' 'reversed'}; 
    for i = 1:length(label_array)
        text(i, ax.YLim(1), sprintf('%s\n%s', label_array{:, i}), 'FontSize', 14, ...
            'horizontalalignment', 'center', 'verticalalignment', 'top');    
    end
    ylabel('GFP (\muV)')
    set(gca, 'Fontsize', 14)

    % plot the markers
    for b = 1:size(data_visual, 2)
        scat(b) = scatter(repelem(b, size(data_visual, 1)), data_visual(:, b),...
            75, colour(b, :), 'filled');
    end

    % mark outliers
    h_out = flipud(findobj(gcf,'tag','Outliers'));
    for h = 1:length(h_out)
        x_out =  get(h_out(h), 'XData');
        y_out =  get(h_out(h), 'YData');
        for i = 1:length(x_out)
            if ~(isnan(x_out(i)))
                outliers(h, i) = find(data_visual(:, h) == y_out(i));
                text(x_out(i) + 0.1, double(y_out(i)), sprintf('%d', outliers(h, i)))
            end
        end
    end
    
    
end

