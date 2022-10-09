%% AG SICI: TEP ANALYSIS + VISUALIZATION
% Written by Dominika for the AG-SICI project, part 2 (2022)
% 
% Colection of scripts to analyse preprocessed TMS-EEG data
% - custom functions are included in the script 
% - letswave6 functions are being called from the directory (add path!)
% 
% Output:
%   --> variables are saved in a global MATLAB output file 'AG-SICI_P2.mat'
%   --> figures are saved in a folder 'AG-SICI_P2_figures'
% 
% 1) PREPARE DATA
%       - load RAGU output for both protocols
%       - load individual data, remove possible outliers
%       - average across selected target channels
%       - merge per category & save for letswave
%       - compute average + SEM 
% 
% 2) TEP BUTTERFLY - individual conditions
%       - colourcoded
%       - optionally highlight selected channel
% 
% 3) GFP - all conditions
%       - average all conditions, compute GFP
%       - plot GFP, identify peak latencies, add topoplots
%       - plot GFP alone --> for MS figure
% 
% 4) GFP - single pulse TEPs
%       - calculate individual GFP
%       - extract mean GFP over predefined TOIs (P25 + N75)
%       - plot boxplot per TOI 
%       - plot mean topography per TOI
% 
% 5) EXTRACT GFP AMPLITUDE - single pulse TEPs
%       - extract peak/mean amplitude and latency from individual GFPs
%       --> 'AGSICI_TEP_GFP'
%       - plot boxplot per GFP peak
%       - plot mean topography per GFP peak
% 
% 6) PLOT GFP & PERFORM RANOVA - single pulse 
%       - test for normality and sphericity
%       - run repeated measures ANOVA 
%       - in case of significant results, run post-hoc paired t-tests
% 
% 7) IDENTIFY EOIs
%       - 3 eois with the highest voltage based on MS class maps 
%       --> 'AGSICI_TEP_default'
% 
% 8) EXTRACT TEP AMPLITUDE 
%       - extract peak/mean amplitude and latency from individual TEPs
%       based on identified EOIs and timing from the grand average GFP
%       --> 'AGSICI_TEP'
%       - saves values for R
%       --> 'AGSICI_TEP_values'
% 
% 9) CALCULATE MEAN VALUES & PERFORM RANOVA
% 
% 10) PLOT ABSOLUTE PEAK VALUES - boxplot
%       - user can choose to plot different comparisons, either only
%       single pulse TEPs ('single'), or TS + paired pulse TEPs ('paired')
%       - plot for each peak
% 
% 11) PLOT SICI AMPLITUDE - boxplot
%       - calculate SICI as difference between paired pulse TEP amd TS
%       TEP amplitude
%       - plot for each peak
% 
% 12) PLOT SICI TIMECOURSE 
%       - calculate SICI as point-by-point subtraction of TS TEPs from
%       paired pulse TEPs
%       - plot butterfly plots per paired pulse condition
%       - plot GFP per paired pulse condition
%       - plot GFP + SEM from a chosen channel 
% 
% 13) EXTRACT SICI GFP AMPLITUDE
%       - extract peak/mean amplitude and latency from individual SICI GFPs
%       --> 'AGSICI_SICI_GFP'
% 
% 14) PLOT SICI GFP & PERFORM RANOVA
%       - plot peak topoplots based on individual GFP latency
%       - plot boxplots per GFP peak
%       - test for normality and sphericity
%       - run repeated measures ANOVA 
%       - in case of significant results, run post-hoc paired t-tests


%% parameters
clear all; clc

% ----------- dataset -----------
% dataset
study = 'P2';
subject = [1:18];
protocol = {'spTMS', 'ppTMS'};
intensity = {'CS1', 'CS2', 'CS3'};
prefix = 'avg bl_correct icfilt ica visual crop notch bandpass sp_filter prea'; 
time_window = [-0.05, 0.3];
outliers = [8, 11, 15];
% --------------------------------

% define conditions
condition{1} = [protocol{1} ' TS'];
for p = 1:length(protocol)
    for i = 1:length(intensity)
        condition{(p-1)*length(intensity) + i + 1} = [protocol{p} ' ' intensity{i}];
    end
end
clear p i

% choose the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');
addpath(folder_git)

% choose the input folder 
folder_input = uigetdir(pwd, 'Coose the input folder');

% choose the results folders + outputs
folder_results = uigetdir(pwd, 'Coose the results folder');
output_file = [folder_results '\AG-SICI_' study '.mat'];
folder_figures = [folder_results '\AG-SICI_' study '_figures'];

% load the info file
load(output_file, 'AGSICI_info')

% load a header
% load([folder_input '\' prefix ' ' study ' 01 ' condition{1} '.lw6'], '-mat')
load([folder_git '\header_example.mat'])

% visualization parameters
figure_counter = 1;
shade = 0.2;
electrodes = {header.chanlocs.labels};
x_step = header.xstep;
x_datastart = header.xstart;
x = [time_window(1):x_step:time_window(2)];
x_start = (time_window(1) - x_datastart)/x_step;
x_end = (time_window(2) - x_datastart)/x_step;

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for a = 1:length(intensity) + 1
           colours(a, :) = uisetcolor; 
        end
    case 'NO'
        if exist(['AGSICI_' study '_colours.mat']) > 0
            load(['AGSICI_' study '_colours.mat'])
        else
            disp('No colour scheme found in this directory!')    
        end
end
save([folder_git '\AGSICI_' study '_colours.mat'], 'colours'); 
clear a answer

%% 1) PREPARE DATA
% ----- section input -----
target = {'CP5' 'P3' 'P5' 'P7'};
% -------------------------
% load RAGU data
load([folder_results '\AG-SICI_' study '_microstates\AGSICI_' study '_MS_all.mat'])
data_MS = rd;
clear rd 

% load individual data, filter out outliers
idx = find(~ismember(subject, outliers));
for c = 1:length(condition)
    for i = 1:length(idx)
        % identify subject
        if subject(idx(i)) < 10
           subj = ['0' num2str(subject(idx(i)))];
        else
           subj = num2str(subject(idx(i))); 
        end

        % load the data, crop
        load([folder_input '\' prefix ' ' study ' ' subj ' ' condition{c} '.mat'])
        data_individual(c, i, :, :) = squeeze(data(1, :, 1, 1, x_start:x_end));
    end
end
clear idx c i subj data

% identify channels to pool
channels = [];
for t = 1:length(target)
    for e = 1:length(electrodes)
        if strcmp(electrodes{e}, target{t})
            channels(end+1) = e;
        end
    end
end
clear t e 

% create the target channel & update electrodes 
data_individual(:, :, size(data_individual, 3) + 1, :) = mean(data_individual(:, :, channels, :), 3);
electrodes{end + 1} = 'target';
clear channels

% save data per category for letswave
for c = 1:length(condition)
    % subset data
    for s = 1:size(data_individual, 2)
        for e = 1:size(data_individual, 3)
            for t = 1:size(data_individual, 4)
                data(s, e, 1, 1, 1, t) = data_individual(c, s, e, t);
            end
        end
    end
       
    % save for lw
    name = sprintf('merged WO AGSICI %s %s', study, condition{c});
    savelw(data, header, name, 'target')        
end
clear c s e t data name

% compute average and SEM
data_mean = squeeze(mean(data_individual, 2));
data_SEM =  squeeze(std(data_individual, 0, 2)/sqrt(size(data_individual, 2)));

% append to the outcome MATLAB file
save(output_file, 'data_MS', 'data_individual', 'data_mean', 'data_SEM', '-append');
clear target 

%% 2) TEP BUTTERFLY - individual conditions
% ----- section input -----
channel = 'Cz';
y_lim = {[-4 4] [0 2.2]};
% -------------------------
% identify channel to highlight
channel_n = find(strcmp(electrodes, channel));

% colour scheme
colours_all = [colours; colours(2:4, :)];

% loop through conditions to plot
for c = 1:length(condition)
    % plot butterfly plot
    data_visual = squeeze(data_mean(c, :, :));
    fig = plot_TEP(figure_counter, x*1000, data_visual, y_lim{1}, 'channel', channel_n, 'colour', colours_all(c, :));
        
    % name figure & save
    fig_name = sprintf('AGSICI_TEP_butterfly_%s', condition{c});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name '.svg'], 'svg')
    figure_counter = figure_counter + 1;
end
clear channel channel_n y_lim colours_all c data_visual fig fig_name

%% 3) GFP 
% ----- section input -----
labeled = 'off';
max_peaks = 6;
y_lim = [0 2.5];
% -------------------------
% ---- all conditions -----
% average across subjects & conditions
data = squeeze(mean(data_individual, [1, 2]));
    
% calculate GFP
gfp = std(data, 1); 
    
% launch the figure
fig = figure(figure_counter);
hold on
    
% extract peak latencies
h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
AG_peaks = GFP_peaks(gfp, time_window*1000, x_step*1000, labeled, 'max_peaks', max_peaks);    

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
    topo_plot(header, data_topoplot, AG_peaks(t)/1000, time_window(1), [-2, 2])

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
figure_name = 'AGSICI_TEP_GFP_all';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])
figure_counter = figure_counter + 1; 

% plot GFP
fig = plot_GFP(figure_counter, x*1000, gfp, y_lim, 'colour', colours(1, :));

% name & save
figure_name = 'AGSICI_TEP_GFP_all';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
figure_counter = figure_counter + 1;

% --- paired-pulse TEPs ---
% average across subjects & conditions
data = squeeze(mean(data_individual([1, 5, 6, 7], :, :, :), [1, 2]));
    
% calculate GFP
gfp = std(data, 1); 

% plot GFP
fig = plot_GFP(figure_counter, x*1000, gfp, y_lim, 'colour', colours(1, :));

% name & save
figure_name = 'AGSICI_TEP_GFP_paired';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
figure_counter = figure_counter + 1;

clear labeled max_peaks data gfp h_axis fig figure_name pos

%% 4) GFP - single pulse TEPs
% ----- section input -----
y_lim = [0 1.9];
peaks = {'P25' 'N75'};
TOI = [15 40; 55 80];
%-------------------------
% calculate GFP for each subject
for c = 1:size(data_individual, 1)
    for s = 1:size(data_individual, 2)
        data_GFP_individual(c, s, :) = std(squeeze(data_individual(c, s, 1:32, :)), 1);
    end
end

% calculate GFP for each condition
for c = 1:size(data_individual, 1)
    data_GFP(c, :) = std(squeeze(mean(data_individual(c, :, 1:32, :), 2)), 1);
end

% extract mean GFP from TOIs 
for c = 1:length(condition)
     for s = 1:size(data_individual, 2)
         for t = 1:size(TOI, 1)
             % calculate TOI limits
             lim = (TOI(t, :) - time_window(1)*1000)/(x_step*1000);
             
             % extract mean GFP
             AGSICI_GFP(c, s, t) = mean(data_GFP_individual(c, s, lim(1):lim(2)));
         end
     end
end

% plot boxplot per TOI and single-pulse condition
for t = 1:size(TOI, 1)    
    % plot GFP
    fig = plot_box(squeeze(AGSICI_GFP([2, 3, 4, 1], :, t))', 'GFP', condition([2, 3, 4, 1]), colours([2, 3, 4, 1], :), figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('AGSICI_TEP_GFP_%s', peaks{t});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end

% plot GFP per condition
colours_all = colours([1, 2, 3, 4, 2, 3, 4], :);
for c = 1:length(condition)
    % plot GFP
    fig = plot_GFP(figure_counter, x*1000, data_GFP(c, :), y_lim, 'colour', colours_all(c, :));

    % name & save
    figure_name = sprintf('AGSICI_TEP_GFP_%s', condition{c});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
    figure_counter = figure_counter + 1;
end

% plot TOI topoplots per condition
cmap = blue2red(48);
for c = 1:length(condition)
    for t = 1:size(TOI, 1)
        % calculate TOI limits
        lim = (TOI(t, :) - time_window(1)*1000)/(x_step*1000);

        % subset the data
        data_topoplot = double(squeeze(mean(data_individual(c, :, 1:32, lim(1):lim(2)), [2, 4])));
        chanlocs = header.chanlocs;

        % plot GFP
        fig = figure(figure_counter);
        topoplot(data_topoplot, chanlocs, 'maplimits', [-1 1], 'shading', 'interp', 'whitebk', 'on',...
            'colormap', cmap, 'style', 'map', 'electrodes', 'off')
        set(gcf,'color',[1 1 1]);

        % name & save
        figure_name = sprintf('AGSICI_TEP_%s_%s', peaks{t}, condition{c});
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
    end
end

% append to the outcome MATLAB file
save(output_file, 'data_GFP_individual', 'data_GFP', 'AGSICI_GFP', '-append');

clear c s t fig y_lim colours_all figure_name TOI lim peaks cmap data_topoplot chanlocs

%% 5) EXTRACT GFP AMPLITUDE - single pulse TEPs
% ----- section input -----
GFP_peak = {'P25' 'N75'};
TOI = [15 40; 55 80];
percent = 20;                           % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                      % y limits for topoplots 
%-------------------------
% calculate extraction parameters
for k = 1:length(GFP_peak)
    GFP_span(k) = (TOI(k, 2) - TOI(k, 1) + 1)/1000;
    GFP_center(k) = TOI(k, 1)/1000 + GFP_span(k)/2;
end

% visualization style
col_fig = [0.9294    0.1412    0.1412; 0.7176    0.2745    1.0000];
cmap = blue2red(120);

% loop through subjects
for s = 1:size(data_individual, 2)
    % setup names   
    idx = find(~ismember(subject, outliers));
    figure_title = sprintf('Subject n. %d', idx(s));    
       
    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;
        
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
    % loop through peaks
    for k = 1:length(GFP_peak)               
        % choose data
        for c = 1:4
            % data for timeseries visualization
            data_visual(c, :) = data_GFP_individual(c, s, :); 

            % data for topoplot
            for e = 1:32
                data_topoplot(c, e, 1, 1, 1, :) = squeeze(data_individual(c, s, e, :));
            end
        end
                       
        % define default TOI 
        center = GFP_center(k);
        span = GFP_span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;   
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 6, 1:18);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, AGSICI_TEP_default.peak{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 6, 1:18)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % extract amplitude and latency
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, x_step, time_window(1), 1);

            % update the figure            
            subplot(4, 6, 1:18)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
            end
            plot(x, data_visual(1, :), 'Color', [0 0 0], 'LineWidth', 2.5)
            plot(x, data_visual(2:4, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 2.5)
            plot(lat_peak, y_max, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')
            
            % add the topoplot  
            for a = 1:4
                subplot(4, 6, 18 + a);
                header.chanlocs = header.chanlocs(1:32);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims, cmap) 
            end

            % ask for approval
            answer = questdlg('Do you want to proceed?', GFP_peak{k},...
                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');

            % switch action
            switch answer
                case 'Yes, extract outcome values.'
                    % close the figure
                    close(fig_1)

                    % exit the while loop
                    finish = 1;

                case 'No, I will adjust TOI manually.'
                    % assign previous center and span
                    choose_center = center;  
                    choose_span = 2 * span;  

                    % identify the limits for visualisation of current peak
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / x_step);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / x_step);
                    choose_x = (choose_center - choose_span/2) : x_step : (choose_center + choose_span/2);

                    % prepare data and header for visualization
                    choose_data = data_visual(:, choose_x1 : choose_x2);
                    choose_header = header;
                    choose_header.datasize(6) = length(choose_data);  
                    choose_header.xstart = choose_center - choose_span/2;

                    % check if vector size matches
                    if size(choose_data, 2) ~= length(choose_x)
                        diff = size(choose_data, 2) - length(choose_x);
                        if diff > 0
                            choose_data = choose_data(:, 1:end - diff);
                        elseif diff < 0
                            choose_x = choose_x(1:end + diff);
                        end
                    end

                    % launch the choosing figure                 
                    choose_figure_name = ['Choose manually peak ' AGSICI_TEP_default.peak{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2, 'Color', [0.0549    0.5216    0.8118])
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    for a = 1
                        choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims, cmap);
                    end           

                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    for a = 1
                        subplot(4, 4, a) 
                        cla(choose_axesHandles(2))
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims, cmap);
                    end
                    hold off

                    % update the central latency
                    center = pos_x;

                    % close the choosing figure
                    pause(2)
                    close(choose_fig)

                    % close the the main figure
                    close(fig_1)
                end
        end
        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
            choose_data choose_center choose_axesHandles answer choose_figure_name

        % record outcome variables
        AGSICI_TEP_GFP.latency(:, s, k) = lat_peak; 
        AGSICI_TEP_GFP.amplitude_peak(:, s, k) = y_max; 
        AGSICI_TEP_GFP.amplitude_mean(:, s, k) = y_mean; 


        % set up the main figure
        figure(fig)
        subplot(length(GFP_peak), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', GFP_peak{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('GFP (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = {'120% rMT' '60% rMT' '80% rMT' '100% rMT'};
        for a = 1:4
            subplot(length(GFP_peak), 7, 3 + a + 7*(axis_counter-1))
            topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims, cmap);   
            if axis_counter == 1
                 set(get(gca, 'title'), 'string', topoplot_titles{a}, 'FontSize', 14);
            end
        end
            
        % update axis counter
        axis_counter = axis_counter + 1;
    end                                       

    % finalize the summary figure
    figure(fig)  
    subplot(length(GFP_peak), 7, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if subject(idx(s)) < 10
        subj = ['0' num2str(subject(idx(s)))];
    else
        subj = num2str(subject(idx(s)));
    end
    figure_name = ['AGSICI_' study '_amplitude_' subj];
    savefig([folder_figures '\TEP-GFP amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\TEP-GFP amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
end   
    
% append to the general MATLAB file
save(output_file, 'AGSICI_TEP_GFP', '-append');     

clear k s c a e TOI GFP_peaks GFP_span GFP_center percent map_lims idx figure_title fig axis_counter center span ...
    data_visual data_topoplot fig_1 yl  y_mean y_max lat_peak col_fig1 col_fig pos finish topoplot_titles ...
    yl_sgtitle figure_name cmap

%% 6) PLOT GFP & PERFORM RANOVA - single-pulse 
% ----- section input -----
GFP_peak = {'P25' 'N75'};
modify = 'off';
%-------------------------
% modify output structure if needed
if strcmp(modify, 'on')
    struct_temp = AGSICI_TEP_GFP;
    clear AGSICI_TEP_GFP
    AGSICI_TEP_GFP.latency.data = struct_temp.latency;
    AGSICI_TEP_GFP.amplitude_peak.data = struct_temp.amplitude_peak;
    AGSICI_TEP_GFP.amplitude_mean.data = struct_temp.amplitude_mean;
    modify = 'off';
    clear struct_temp 
end

% plot boxplot per peak 
for k = 1:length(GFP_peak)    
    % plot GFP
    fig = plot_box(squeeze(AGSICI_TEP_GFP.amplitude_mean.data([2, 3, 4, 1], :, k))', 'GFP', condition([2, 3, 4, 1]), colours([2, 3, 4, 1], :), figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('AGSICI_TEP_GFP_individualized_%s', GFP_peak{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end

% calculate mean data + plot peak topoplots
cmap = blue2red(48);
for c = 1:4
    for k = 1:length(GFP_peak) 
        % extract individual values
        for s = 1:size(data_individual, 2)
            peak_x(c, k, s) = (AGSICI_TEP_GFP.latency.data(c, s, k) - time_window(1))/x_step;
            peak_amplitude(c, k, s) = AGSICI_TEP_GFP.amplitude_peak.data(c, s, k);
            peak_latency(c, k, s) = AGSICI_TEP_GFP.latency.data(c, s, k) * 1000;
        end
        
        % compute mean values
        AGSICI_TEP_GFP.amplitude_peak.mean(c, k) = mean(peak_amplitude(c, k, :));
        AGSICI_TEP_GFP.amplitude_peak.SD(c, k) = std(peak_amplitude(c, k, :));
        AGSICI_TEP_GFP.amplitude_peak.SEM(c, k) = std(peak_amplitude(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_TEP_GFP.amplitude_peak.extremes(c, k, :) = [min(peak_amplitude(c, k, :)) max(peak_amplitude(c, k, :))];
        AGSICI_TEP_GFP.latency.mean(c, k) = mean(peak_latency(c, k, :));
        AGSICI_TEP_GFP.latency.SD(c, k) = std(peak_latency(c, k, :));
        AGSICI_TEP_GFP.latency.SEM(c, k) = std(peak_latency(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_TEP_GFP.latency.extremes(c, k, :) = [min(peak_latency(c, k, :)) max(peak_latency(c, k, :))];
        
        % extract peak topographies normalized by peak GFP
        data_topoplot = [];
        for s = 1:size(data_individual, 2)
            data_topoplot(s, :) = data_individual(c, s, 1:32, ceil(peak_x(c, k, s)))/peak_amplitude(c, k, s);
        end
        
        % prepare data 
        data_topoplot = double(mean(data_topoplot, 1));
        chanlocs = header.chanlocs;

        % plot peak topoplot
        fig = figure(figure_counter);
        topoplot(data_topoplot, chanlocs, 'maplimits', [-1 1], 'shading', 'interp', 'whitebk', 'on',...
            'colormap', cmap, 'style', 'map', 'electrodes', 'off')
        set(gcf,'color',[1 1 1]);

        % name & save
        figure_name = sprintf('AGSICI_TEP_topo_%s_%s', GFP_peak{k}, condition{c});
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
    end
end

% run RANOVA for amplitude and latency
idx = find(~ismember(subject, outliers));
condition_single = {'spTS' 'spCS1' 'spCS2' 'spCS3'};
for k = 1:length(GFP_peak)
    % ----- amplitude ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = [2, 3, 4, 1]
        statement = ['data_table.' condition_single{c} ' = AGSICI_TEP_GFP.amplitude_peak.data(c, :, k)'';'];
        eval(statement)
    end
    
    % within-subjects design
    wd = table(condition_single([2, 3, 4, 1])','VariableNames',{'stimulus'});
    wd.stimulus = categorical(wd.stimulus); 
    
    % rm-ANOVA
    rm = fitrm(data_table, 'spCS1-spTS ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' GFP_peak{k} '_amplitude = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)
    
    % extract p
    if rm.mauchly.pValue > 0.05
        statement = ['p = ranova_' GFP_peak{k} '_amplitude.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' GFP_peak{k} '_amplitude.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' GFP_peak{k} '_amplitude  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
    
    % ----- latency ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = [2, 3, 4, 1]
        statement = ['data_table.' condition_single{c} ' = AGSICI_TEP_GFP.latency.data(c, :, k)'';'];
        eval(statement)
    end
    
    % rm-ANOVA
    rm = fitrm(data_table, 'spCS1-spTS ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' GFP_peak{k} '_latency = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)    
    
        % extract p
    if rm.mauchly.pValue > 0.05
        statement = ['p = ranova_' GFP_peak{k} '_latency.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' GFP_peak{k} '_latency.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' GFP_peak{k} '_latency  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
end

% append to the general MATLAB file
save(output_file, 'AGSICI_TEP_GFP', '-append');  

clear c k s fig figure_name GFP_peak cmap peak_x peak_amplitude peak_latency data_topoplot chanlocs ...
    idx data_table condition_single statement wd rm p

%% 7) IDENTIFY EOIs 
% launch the default structure
eoi_n = [3,3,3,3,1];
AGSICI_TEP_default = struct;
AGSICI_TEP_default.peak = {'P20' 'P30' 'N75' 'N100' 'P180'};
AGSICI_TEP_default.center = [0.023 0.031 0.073 0.116 0.188];
AGSICI_TEP_default.span = [0.014 0.03 0.03 0.04 0.06];

% extract voltage from MS maps
MS_voltage = data_MS.MSMaps([1, 2, 3, 3, 5], :);

% choose EOIs based on the voltage
for k = 1:length(AGSICI_TEP_default.peak)
    % index max/min values
    if strcmp(AGSICI_TEP_default.peak{k}(1), 'P')
        [eoi_val, eoi_i] = maxk(MS_voltage(k, :), eoi_n(k));
    else
        [eoi_val, eoi_i] = mink(MS_voltage(k, :), eoi_n(k));
    end
    
    % identify EOI labels
    AGSICI_TEP_default.eoi{k} = electrodes(eoi_i);
end
clear k eoi_val eoi_i eoi_n

% save to the global MATLAB file
save(output_file, 'AGSICI_TEP_default', '-append');

%% 8) EXTRACT TEP AMPLITUDE
% ----- section input -----
percent = 20;                           % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                      % y limits for topoplots 
% -------------------------
% linestyle
col_fig = [1.0000    0.4118    0.1608; 0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.7176    0.2745    1.0000; 0.3647    0.2078    0.9882];

% loop through subjects
for s = 1:size(data_individual, 2)
    % setup names   
    idx = find(~ismember(subject, outliers));
    figure_title = sprintf('Subject n. %d', idx(s));    
       
    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;
        
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
    % loop through peaks
    for k = 1:length(AGSICI_TEP_default.peak)   
        % identify peak polarity
        if strcmp(AGSICI_TEP_default.peak{k}(1), 'P')
            polarity = 1;
        else
            polarity = -1;
        end
            
        % identify EOIs
        eoi = [];
        for e = 1:length(AGSICI_TEP_default.eoi{k})
            eoi(e) = find(strcmp(electrodes, AGSICI_TEP_default.eoi{k}{e}));
        end
            
        % choose data
        for c = 1:length(condition)
            % data for timeseries visualization
            data_visual(c, :) = squeeze(mean(data_individual(c, s, eoi, :), 3)); 

            % data for topoplot
            for e = 1:32
                data_topoplot(c, e, 1, 1, 1, :) = squeeze(data_individual(c, s, e, :));
            end
        end
                       
        % define default TOI 
        center = AGSICI_TEP_default.center(k);
        span = AGSICI_TEP_default.span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;                
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 6, 1:18);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, AGSICI_TEP_default.peak{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 6, 1:18)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % extract amplitude and latency
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, x_step, time_window(1), polarity);

            % update the figure            
            subplot(4, 6, 1:18)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
            end
            plot(x, data_visual(1, :), 'Color', [0 0 0], 'LineWidth', 2.5)
            plot(x, data_visual(2:4, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 2.5)
            plot(x, data_visual(5:7, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 2.5)
            plot(lat_peak, y_max, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')

            % add the topoplot   
            for a = 1:6
                subplot(4, 6, 18 + a);
                header.chanlocs = header.chanlocs(1:32);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims) 
            end

            % ask for approval
            answer = questdlg('Do you want to proceed?', AGSICI_TEP_default.peak{k},...
                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');

            % switch action
            switch answer
                case 'Yes, extract outcome values.'
                    % close the figure
                    close(fig_1)

                    % exit the while loop
                    finish = 1;

                case 'No, I will adjust TOI manually.'
                    % assign previous center and span
                    choose_center = center;  
                    choose_span = 2 * span;  

                    % identify the limits for visualisation of current peak
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / x_step);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / x_step);
                    choose_x = (choose_center - choose_span/2) : x_step : (choose_center + choose_span/2);

                    % prepare data and header for visualization
                    choose_data = data_visual(:, choose_x1 : choose_x2);
                    choose_header = header;
                    choose_header.datasize(6) = length(choose_data);  
                    choose_header.xstart = choose_center - choose_span/2;

                    % check if vector size matches
                    if size(choose_data, 2) ~= length(choose_x)
                        diff = size(choose_data, 2) - length(choose_x);
                        if diff > 0
                            choose_data = choose_data(:, 1:end - diff);
                        elseif diff < 0
                            choose_x = choose_x(1:end + diff);
                        end
                    end

                    % launch the choosing figure                 
                    choose_figure_name = ['Choose manually peak ' AGSICI_TEP_default.peak{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2, 'Color', [0.0549    0.5216    0.8118])
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    for a = 1
                        choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims);
                    end           

                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    for a = 1
                        subplot(4, 4, a) 
                        cla(choose_axesHandles(2))
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims);
                    end
                    hold off

                    % update the central latency
                    center = pos_x;

                    % close the choosing figure
                    pause(2)
                    close(choose_fig)

                    % close the the main figure
                    close(fig_1)
                end
        end
        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
            choose_data choose_center choose_axesHandles answer choose_figure_name

        % record outcome variables
        AGSICI_TEP.latency(:, s, k) = lat_peak; 
        AGSICI_TEP.amplitude_peak(:, s, k) = y_max; 
        AGSICI_TEP.amplitude_mean(:, s, k) = y_mean; 


        % set up the main figure
        figure(fig)
        subplot(length(AGSICI_TEP_default.peak), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', AGSICI_TEP_default.peak{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = {'120% rMT' '60% rMT' '80% rMT' '100% rMT'};
        for a = 1:4
            subplot(length(AGSICI_TEP_default.peak), 7, 3 + a + 7*(axis_counter-1))
            topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims);   
            if axis_counter == 1
                 set(get(gca, 'title'), 'string', topoplot_titles{a}, 'FontSize', 14);
            end
        end
            
        % update axis counter
        axis_counter = axis_counter + 1;
    end                                       

    % finalize the summary figure
    figure(fig)  
    subplot(length(AGSICI_TEP_default.peak), 7, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if subject(idx(s)) < 10
        subj = ['0' num2str(subject(idx(s)))];
    else
        subj = num2str(subject(idx(s)));
    end
    figure_name = ['AGSICI_' study '_amplitude_' subj];
    savefig([folder_figures '\TEP amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\TEP amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
end   
    
% append to the general MATLAB file
save(output_file, 'AGSICI_TEP', '-append');
clear s k e c a idx figure_title data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak ...
    col_fig1 col_fig pos finish axis_counter topoplot_titles eoi yl_sgtitle figure_name percent polarity map_lims

% save for R
idx = find(~ismember(subject, outliers));
AGSICI_TEP_values = table;
row_counter = height(AGSICI_TEP_values) + 1;
for s = 1:size(AGSICI_TEP.amplitude_peak, 2)    
    for c = 1:length(condition)  
        for k = 1:length(AGSICI_TEP_default.peak) 
            %fill in the table
            AGSICI_TEP_values.subject(row_counter) = idx(s);
            AGSICI_TEP_values.condition(row_counter) = condition(c);
            AGSICI_TEP_values.peak(row_counter) = AGSICI_TEP_default.peak(k);
            AGSICI_TEP_values.amplitude_peak(row_counter) = AGSICI_TEP.amplitude_peak(c, s, k);
            AGSICI_TEP_values.amplitude_mean(row_counter) = AGSICI_TEP.amplitude_mean(c, s, k);
            AGSICI_TEP_values.latency(row_counter) = AGSICI_TEP.latency(c, s, k);

            % update the counter
            row_counter = row_counter + 1;
        end
    end
end
writetable(AGSICI_TEP_values, [folder_results '\AGSICI_' study '_TEP.csv'])

%% 9) CALCULATE MEAN VALUES & PERFORM RANOVA
% ----- section input -----
TEP_peak = {'P20' 'P30' 'N75'};
modify = 'on';
% -------------------------
% modify output structure if needed
if strcmp(modify, 'on')
    struct_temp = AGSICI_TEP;
    clear AGSICI_TEP
    AGSICI_TEP.latency.data = struct_temp.latency;
    AGSICI_TEP.amplitude_peak.data = struct_temp.amplitude_peak;
    AGSICI_TEP.amplitude_mean.data = struct_temp.amplitude_mean;
    modify = 'off';
    clear struct_temp 
end

% calculate mean data 
for c = 1:length(condition)
    for k = 1:5
        % extract individual values
        for s = 1:size(data_individual, 2)
            peak_x(c, k, s) = (AGSICI_TEP.latency.data(c, s, k) - time_window(1))/x_step;
            peak_amplitude(c, k, s) = AGSICI_TEP.amplitude_peak.data(c, s, k);
            peak_latency(c, k, s) = AGSICI_TEP.latency.data(c, s, k) * 1000;
        end
        
        % compute mean values
        AGSICI_TEP.amplitude_peak.mean(c, k) = mean(peak_amplitude(c, k, :));
        AGSICI_TEP.amplitude_peak.SD(c, k) = std(peak_amplitude(c, k, :));
        AGSICI_TEP.amplitude_peak.SEM(c, k) = std(peak_amplitude(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_TEP.amplitude_peak.extremes(c, k, :) = [min(peak_amplitude(c, k, :)) max(peak_amplitude(c, k, :))];
        AGSICI_TEP.latency.mean(c, k) = mean(peak_latency(c, k, :));
        AGSICI_TEP.latency.SD(c, k) = std(peak_latency(c, k, :));
        AGSICI_TEP.latency.SEM(c, k) = std(peak_latency(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_TEP.latency.extremes(c, k, :) = [min(peak_latency(c, k, :)) max(peak_latency(c, k, :))];
    end
end

% run RANOVA for amplitude and latency
idx = find(~ismember(subject, outliers));
names = {'spTS' 'CS1_TS' 'CS2_TS' 'CS3_TS'};
for k = 1:length(TEP_peak)
    % ----- amplitude ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = 1:length(names)
        if c == 1
            statement = ['data_table.' names{c} ' = AGSICI_TEP.amplitude_peak.data(c, :, k)'';'];
            eval(statement)
        else
            statement = ['data_table.' names{c} ' = AGSICI_TEP.amplitude_peak.data(3 + c, :, k)'';'];
            eval(statement)
        end
    end
    
    % within-subjects design
    wd = table(names','VariableNames',{'stimulus'});
    wd.stimulus = categorical(wd.stimulus); 
    
    % rm-ANOVA
    rm = fitrm(data_table, 'spTS-CS3_TS ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' TEP_peak{k} '_amplitude = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)
    
    % extract p
    if rm.mauchly.pValue > 0.07
        statement = ['p = ranova_' TEP_peak{k} '_amplitude.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' TEP_peak{k} '_amplitude.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' TEP_peak{k} '_amplitude  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
    
    % ----- latency ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = 1:length(names)
        if c == 1
            statement = ['data_table.' names{c} ' = AGSICI_TEP.latency.data(c, :, k)'';'];
            eval(statement)
        else
            statement = ['data_table.' names{c} ' = AGSICI_TEP.latency.data(3 + c, :, k)'';'];
            eval(statement)
        end
    end
    
    % rm-ANOVA
    rm = fitrm(data_table, 'spTS-CS3_TS ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' TEP_peak{k} '_latency = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)    
    
        % extract p
    if rm.mauchly.pValue > 0.05
        statement = ['p = ranova_' TEP_peak{k} '_latency.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' TEP_peak{k} '_latency.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' TEP_peak{k} '_latency  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
end

% append the output variables to the general MATLAB file
save(output_file, 'AGSICI_TEP', '-append');

clear TEP_peak c k s peak_x peak_amplitude peak_latency names modify...
    fig figure_name idx data_table statement wd rm p

%% 10) PLOT ABSOLUTE PEAK VALUES - boxplot
% ----- section input -----
peaks2plot = 1:5;
comparison = 'paired';
% -------------------------
% load input if necessary
if exist('AGSICI_TEP') ~= 1
    load(output_file, 'AGSICI_TEP', 'AGSICI_TEP_default')
end

% choose datasets to display
switch comparison
    case 'single'
        comp_conditions = [2, 3, 4, 1];
        col = colours(comp_conditions, :);
    case 'paired'
        comp_conditions = [1, 5, 6, 7];
        col = colours;
end

% loop through peaks
for k = peaks2plot
    % load data
    data_amplitude = []; data_latency = [];
    for c = 1:length(comp_conditions)     
        data_amplitude(c, :) = squeeze(AGSICI_TEP.amplitude_peak(comp_conditions(c), :, k));
        data_latency(c, :) = squeeze(AGSICI_TEP.latency(comp_conditions(c), :, k));
    end

    % plot amplitude
    fig = plot_box(data_amplitude', 'amplitude', condition(comp_conditions), col, figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('AGSICI_TEP_amplitude_%s_%s', comparison, AGSICI_TEP_default.peak{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear k c fig figure_name data_amplitude data_latency comparison peaks2plot comp_conditions col

%% 11) PLOT SICI AMPLITUDE - boxplot
% ----- section input -----
peaks2plot = 1:5;
% -------------------------
% load input if necessary
if exist('AGSICI_TEP') ~= 1
    load(output_file, 'AGSICI_TEP', 'AGSICI_TEP_default')
end

% calculate amplitude change due to SICI
for c = 1:3
    for s = 1:size(AGSICI_TEP.amplitude_peak, 2)
        for k = 1:length(AGSICI_TEP_default.peak)
            AGSICI_TEP.SICI(c, s, k) = AGSICI_TEP.amplitude_peak(4 + c, s, k) - AGSICI_TEP.amplitude_peak(1, s, k);
        end
    end
end
save(output_file, 'AGSICI_TEP', '-append');

% loop through peaks
for k = peaks2plot
    % subset data
    data_amplitude = squeeze(AGSICI_TEP.SICI(:, :, k));

    % plot amplitude
    fig = plot_box(data_amplitude', 'SICI', condition([5:7]), colours(2:4, :), figure_counter)
    hold off

    % name and save figure
    figure_name = sprintf('AGSICI_TEP_SICI_%s', AGSICI_TEP_default.peak{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear k c fig figure_name data_amplitude peaks2plot

%% 12) PLOT SICI TIMECOURSE 
% ----- section input -----
channel = {'Pz'}; 
y_lim = [-1.9 1.9];
% -------------------------
% compute SICI timecourse
for c = 1:3
    % calculate mean data
    data = [];
    for s = 1:size(data_individual, 2)
        for e = 1:size(data_individual, 3)
            for i = 1:size(data_individual, 4)
                data_SICI(c, s, e, i) = data_individual(4 + c, s, e, i) - data_individual(1, s, e, i);
                data(s, e, 1, 1, 1, i) = data_SICI(c, s, e, i);
            end
        end
    end
    
%     % save for lwtswave
%     name = sprintf('merged SICI %s', condition{4 + c});
%     savelw(data, header, name, [])
end
clear c s e i name data

% plot SICI as butterfly plot
for e = 1:length(channel)
    % identify the channel
    channel_n = find(strcmp(electrodes, channel{e}));
    
    % plot separately for each CS
    for c = 1:3
        %  butterfly plot
        data_visual = squeeze(mean(data_SICI(c, :, :, :), 2));
        fig = plot_TEP(figure_counter, x*1000, data_visual, y_lim, 'channel', channel_n, 'colour', colours(1 + c, :));

        % name figure & save
        figure_name = sprintf('AGSICI_SICI_butterfly_%s', intensity{c});
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
    end
end
clear e c channel_n data_visual fig y_lim figure_name

% compute GFP of SICI & plot
for c = 1:3
    % calculate GFP
    gfp(c, :) = std(squeeze(mean(data_SICI(c, :, :, :), 2)), 1);
    
    % plot 
    fig = plot_GFP(figure_counter, x*1000, gfp(c, :), [0, 1.5], 'colour', colours(1 + c, :));

    % name & save
    figure_name = sprintf('AGSICI_SICI_GFP_%s', intensity{c});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
    figure_counter = figure_counter + 1;
end

% plot SICI for the chosen channel
for e = 1:length(channel)
    % identify the channel
    e_n = find(strcmp(electrodes, channel{e}));
        
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % cycle through conditions
    for c = 1:3
        % prepare the data
        data_visual = squeeze(mean(data_SICI(c, :, e_n, :), 2))';
        SEM_visual = squeeze(std(data_SICI(c, :, e_n, :), 0, 2)/sqrt(size(data_individual, 2)))';

        % plot        
        P(c) = plot(x*1000, data_visual, 'Color', colours(1 + c, :), 'LineWidth', 3);
        F(c) = fill([x*1000 fliplr(x*1000)],[data_visual + SEM_visual fliplr(data_visual - SEM_visual)], ...
            colours(1 + c, :), 'FaceAlpha', shade, 'linestyle', 'none');
    end
   
    % shade interpolated interval 
    yl = get(gca, 'ylim'); 
    rectangle('Position', [-5, yl(1), 15, yl(2) - yl(1)], 'FaceColor', [0.75 0.73 0.73], 'EdgeColor', 'none')

    % mark TMS stimulus and zero line
    line([0, 0], yl, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    line(x*1000, zeros(1, length(x)), 'Color', [0 0 0], 'LineWidth', 1.5, 'LineStyle', ':')

    % add other parameters
    xlabel('time (ms)')
    ylabel('\Delta amplitude (\muV \pm SEM)')
    set(gca, 'FontSize', 14)
    set(gca, 'Layer', 'Top')
    xlim([x(1)*1000 x(end)*1000])

    % save figure
    figure_name = sprintf('AGSICI_SICI_%s', channel{e});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'])
    figure_counter = figure_counter + 1 ;
end

% append the output variables to the general MATLAB file
save(output_file, 'data_SICI', '-append');
clear channel c e e_n data_visual SEM_visual fig yl P F figure_name 

%% 13) EXTRACT SICI GFP AMPLITUDE
% ----- section input -----
GFP_peak = {'15ms' '50ms' '100ms' '180ms'};
TOI = [10 24; 35 55; 80 120; 150 200];
percent = 20;                               
map_lims = [-1.5 1.5];                   
%-------------------------
% calculate extraction parameters
for k = 1:length(GFP_peak)
    GFP_span(k) = (TOI(k, 2) - TOI(k, 1) + 1)/1000;
    GFP_center(k) = TOI(k, 1)/1000 + GFP_span(k)/2;
end

% visualization style
col_fig = [0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.3647    0.2078    0.9882];
cmap = blue2red(120);

% loop through subjects
for s = 1:size(data_individual, 2)
    % setup names   
    idx = find(~ismember(subject, outliers));
    figure_title = sprintf('Subject n. %d', idx(s));    
       
    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;
        
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
    % loop through peaks
    for k = 1:length(GFP_peak)               
        % choose data
        for c = 1:3
            % data for timeseries visualization
            data_visual(c, :) = std(squeeze(data_SICI(c, s, :, :)), 1);

            % data for topoplot
            for e = 1:32
                data_topoplot(c, e, 1, 1, 1, :) = squeeze(data_SICI(c, s, e, :));
            end
        end
                       
        % define default TOI 
        center = GFP_center(k);
        span = GFP_span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;   
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 6, 1:18);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, GFP_peak{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 6, 1:18)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % extract amplitude and latency
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, x_step, time_window(1), 1);

            % update the figure            
            subplot(4, 6, 1:18)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
            end
            for c = 1:3
                plot(x, data_visual(c, :), 'Color', [0.55 0.55 0.55] - (c - 1)*0.2, 'LineWidth', 2.5)
            end
            plot(lat_peak, y_max, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')
            
            % add the topoplot  
            for a = 1:3
                subplot(4, 6, 18 + a);
                header.chanlocs = header.chanlocs(1:32);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims, cmap) 
            end

            % ask for approval
            answer = questdlg('Do you want to proceed?', GFP_peak{k},...
                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');

            % switch action
            switch answer
                case 'Yes, extract outcome values.'
                    % close the figure
                    close(fig_1)

                    % exit the while loop
                    finish = 1;

                case 'No, I will adjust TOI manually.'
                    % assign previous center and span
                    choose_center = center;  
                    choose_span = 2 * span;  

                    % identify the limits for visualisation of current peak
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / x_step);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / x_step);
                    choose_x = (choose_center - choose_span/2) : x_step : (choose_center + choose_span/2);

                    % prepare data and header for visualization
                    choose_data = data_visual(:, choose_x1 : choose_x2);
                    choose_header = header;
                    choose_header.datasize(6) = length(choose_data);  
                    choose_header.xstart = choose_center - choose_span/2;

                    % check if vector size matches
                    if size(choose_data, 2) ~= length(choose_x)
                        diff = size(choose_data, 2) - length(choose_x);
                        if diff > 0
                            choose_data = choose_data(:, 1:end - diff);
                        elseif diff < 0
                            choose_x = choose_x(1:end + diff);
                        end
                    end

                    % launch the choosing figure                 
                    choose_figure_name = ['Choose manually peak ' GFP_peak{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2, 'Color', [0.0549    0.5216    0.8118])
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    for a = 1:3
                        choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims, cmap);
                    end           

                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    for a = 1
                        subplot(4, 4, a) 
                        cla(choose_axesHandles(2))
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims, cmap);
                    end
                    hold off

                    % update the central latency
                    center = pos_x;

                    % close the choosing figure
                    pause(2)
                    close(choose_fig)

                    % close the the main figure
                    close(fig_1)
                end
        end
        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
            choose_data choose_center choose_axesHandles answer choose_figure_name

        % record outcome variables
        AGSICI_SICI_GFP.latency(:, s, k) = lat_peak; 
        AGSICI_SICI_GFP.amplitude_peak(:, s, k) = y_max; 
        AGSICI_SICI_GFP.amplitude_mean(:, s, k) = y_mean; 

        % set up the main figure
        figure(fig)
        subplot(length(GFP_peak), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', GFP_peak{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('GFP (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = intensity;
        for a = 1:3
            subplot(length(GFP_peak), 7, 3 + a + 7*(axis_counter-1))
            topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims, cmap);   
            if axis_counter == 1
                 set(get(gca, 'title'), 'string', topoplot_titles{a}, 'FontSize', 14);
            end
        end
            
        % update axis counter
        axis_counter = axis_counter + 1;
    end                                       

    % finalize the summary figure
    figure(fig)  
    subplot(length(GFP_peak), 7, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if subject(idx(s)) < 10
        subj = ['0' num2str(subject(idx(s)))];
    else
        subj = num2str(subject(idx(s)));
    end
    figure_name = ['AGSICI_' study '_amplitude_' subj];
    savefig([folder_figures '\SICI-GFP amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\SICI-GFP amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
end   
    
% append progressively the output variables to the general MATLAB file
save(output_file, 'AGSICI_SICI_GFP', '-append');  

clear k s c a e TOI GFP_peak GFP_span GFP_center percent map_lims idx figure_title fig axis_counter center span ...
    data_visual data_topoplot fig_1 yl  y_mean y_max lat_peak col_fig1 col_fig pos finish topoplot_titles ...
    yl_sgtitle figure_name cmap

%% 14) PLOT SICI GFP & PERFORM RANOVA
% ----- section input -----
GFP_peak = {'15ms' '50ms' '100ms' '180ms'};
modify = 'off';
y_lim_topo = [-1 1]; 
y_lim_box = [0 3.5];
%-------------------------
% modify output structure if needed
if strcmp(modify, 'on')
    struct_temp = AGSICI_SICI_GFP;
    clear AGSICI_SICI_GFP
    AGSICI_SICI_GFP.latency.data = struct_temp.latency;
    AGSICI_SICI_GFP.amplitude_peak.data = struct_temp.amplitude_peak;
    AGSICI_SICI_GFP.amplitude_mean.data = struct_temp.amplitude_mean;
    modify = 'off';
    clear struct_temp 
end

% calculate mean data + plot peak topoplots
cmap = blue2red(48);
for c = 1:3
    for k = 1:length(GFP_peak) 
        % extract individual values
        for s = 1:size(data_individual, 2)
            peak_x(c, k, s) = (AGSICI_SICI_GFP.latency.data(c, s, k) - time_window(1))/x_step;
            peak_amplitude(c, k, s) = AGSICI_SICI_GFP.amplitude_peak.data(c, s, k);
            peak_latency(c, k, s) = AGSICI_SICI_GFP.latency.data(c, s, k) * 1000;
        end
        
        % compute mean values
        AGSICI_SICI_GFP.amplitude_peak.mean(c, k) = mean(peak_amplitude(c, k, :));
        AGSICI_SICI_GFP.amplitude_peak.SD(c, k) = std(peak_amplitude(c, k, :));
        AGSICI_SICI_GFP.amplitude_peak.SEM(c, k) = std(peak_amplitude(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_SICI_GFP.amplitude_peak.extremes(c, k, :) = [min(peak_amplitude(c, k, :)) max(peak_amplitude(c, k, :))];
        AGSICI_SICI_GFP.latency.mean(c, k) = mean(peak_latency(c, k, :));
        AGSICI_SICI_GFP.latency.SD(c, k) = std(peak_latency(c, k, :));
        AGSICI_SICI_GFP.latency.SEM(c, k) = std(peak_latency(c, k, :))/sqrt(size(data_individual, 2));
        AGSICI_SICI_GFP.latency.extremes(c, k, :) = [min(peak_latency(c, k, :)) max(peak_latency(c, k, :))];
        
        % extract peak topographies normalized by peak GFP
        data_topoplot = [];
        for s = 1:size(data_individual, 2)
            data_topoplot(s, :) = data_SICI(c, s, 1:32, ceil(peak_x(c, k, s)))/peak_amplitude(c, k, s);
        end
        
        % prepare data 
        data_topoplot = double(mean(data_topoplot, 1));
        chanlocs = header.chanlocs;

        % plot peak topoplot
        fig = figure(figure_counter);
        topoplot(data_topoplot, chanlocs, 'maplimits', y_lim_topo, 'shading', 'interp', 'whitebk', 'on',...
            'colormap', cmap, 'style', 'map', 'electrodes', 'off')
        set(gcf,'color',[1 1 1]);

        % name & save
        figure_name = sprintf('AGSICI_SICI_topo_%s_%s', GFP_peak{k}, intensity{c});
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        figure_counter = figure_counter + 1;
    end
end

% plot boxplots
for k = 1:length(GFP_peak)    
    % plot GFP amplitude
    fig = plot_box(squeeze(AGSICI_SICI_GFP.amplitude_mean.data(:, :, k))', 'GFP', intensity, colours([2, 3, 4], :), figure_counter)
    ylim(y_lim_box)
    hold off

    % name and save figure
    figure_name = sprintf('AGSICI_SICI_GFP_amplitude_%s', GFP_peak{k});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')   

    % update figure counter
    figure_counter = figure_counter + 1;
end

% run RANOVA for amplitude and latency
idx = find(~ismember(subject, outliers));
for k = 1:length(GFP_peak)
    % ----- amplitude ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = 1:3
        statement = ['data_table.' intensity{c} ' = AGSICI_SICI_GFP.amplitude_peak.data(c, :, k)'';'];
        eval(statement)
    end
    
    % within-subjects design
    wd = table(intensity','VariableNames',{'stimulus'});
    wd.stimulus = categorical(wd.stimulus); 
    
    % rm-ANOVA
    rm = fitrm(data_table, 'CS1-CS3 ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' GFP_peak{k} '_amplitude = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)
    
    % extract p
    if rm.mauchly.pValue > 0.07
        statement = ['p = ranova_' GFP_peak{k} '_amplitude.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' GFP_peak{k} '_amplitude.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' GFP_peak{k} '_amplitude  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
    
    % ----- latency ----- 
    % prepare data table
    data_table = table; 
    data_table.subject = idx'; 
    for c = 1:3
        statement = ['data_table.' intensity{c} ' = AGSICI_SICI_GFP.latency.data(c, :, k)'';'];
        eval(statement)
    end
    
    % rm-ANOVA
    rm = fitrm(data_table, 'CS1-CS3 ~ 1', 'WithinDesign', wd);
    
    % check for approx. normal distribution
    figure(figure_counter)
    histogram(reshape(data_table{:,2:end},[],1), 10)
    figure_counter = figure_counter + 1;
    
    % check for sphericity
    rm.mauchly
    
    % results table
    statement = ['ranova_' GFP_peak{k} '_latency = ranova(rm, ''WithinModel'', ''stimulus'');'];
    eval(statement)    
    
        % extract p
    if rm.mauchly.pValue > 0.05
        statement = ['p = ranova_' GFP_peak{k} '_latency.pValue(3);'];
        eval(statement)
    else
        statement = ['p = ranova_' GFP_peak{k} '_latency.pValueGG(3);'];
        eval(statement)
    end
    
    % if ssignificant, compute post-hoc paired t-tests
    if p < 0.05
        statement = ['posthoc_' GFP_peak{k} '_latency  = multcompare(rm,''stimulus'');'];
        eval(statement)
    end
end

% append the output variables to the general MATLAB file
save(output_file, 'AGSICI_SICI_GFP', '-append');
clear GFP_peak c k s peak_x peak_amplitude peak_latency data_topoplot cmap y_lim_topo y_lim_box chanlocs ...
    fig figure_name idx data_table statement wd rm p

%% functions
function savelw(data, header, name, target)
    % modify header
    header.name = name;
    header.datasize = size(data);
    header.xstart = -0.05;
    for e = 1:size(data, 1)
        header.events(e).code = 'AG_TEP';
        header.events(e).latency = 0;
        header.events(e).epoch = e;
    end
    
    % add target if necessary
    if strcmp(target, 'target')
        header.chanlocs(length(header.chanlocs)+1).labels = 'target';
    end
    
    % save 
    save([name '.mat'], 'data')
    save([name '.lw6'], 'header')
end
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
    
    % shade interpolated interval 
    rectangle('Position', [int(1), y_lim(1), int(2) - int(1), y_lim(2) - y_lim(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')

    % loop through channels to plot
    for a = 1:size(data_visual, 1)     
        P(a) = plot(x, data_visual(a, :), 'Color', colour, 'LineWidth', 1);
    end

    % highlight channel
    if ~isempty(channel_n)
        P(end+1) =  plot(x, data_visual(channel_n, :), 'Color', [0 0 0], 'LineWidth', 3);
    end

    % TMS stimulus
    line([0, 0], y_lim, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % other parameters
    xlim([x(1) - length(x)*(x(2) - x(1))*0.05, x(end) + length(x)*(x(2) - x(1))*0.05])
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
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
    [pks, locs] = findpeaks(y, 'MinPeakDistance', 5, 'MinPeakProminence', 0.01);
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
function topo_plot(header, data, x_pos, x_start, map_lims, cmap)
    % define varargins    
    varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on' 'colormap' cmap 'style' 'map' 'electrodes' 'off'};

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
    
    % shade interpolated interval 
    rectangle('Position', [-5, y_lim(1), 15, y_lim(2) - y_lim(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')
    
    % shade GFP
    F = fill([x fliplr(x)],[data_visual zeros(1, length(x))], colour, 'linestyle', 'none');

    % plot GFP line
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);

    % TMS stimulus
    line([0, 0], y_lim, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
    
    % other parameters
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
end
function [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, xstart, polarity)
% identify TOI data
start = ceil(((center-span/2) - xstart)/xstep) + 1;
stop = ceil(((center+span/2) - xstart)/xstep);
x_TOI = start : stop;
y_TOI = data_visual(:, x_TOI);

% calculate number of points to average
x_include_n = ceil((percent/100) * size(y_TOI, 2));

% identify peak values for all datasets  
if polarity > 0
    [y_max, x_max] = max(y_TOI, [], 2);
else
    [y_max, x_max] = min(y_TOI, [], 2);
end
x_peak = start + x_max - 1; 

% identify timepoints that will be included in the mean amplitude
x_TOI_include = [];
for a = 1:size(data_visual, 1)
    x_TOI_include_i = x_max(a);
    t_left = ceil((x_include_n - 1)/2); tx_left = x_max(a) - t_left : x_max(a) - 1;
    t_right = x_include_n - 1 - t_left; tx_right = x_max(a) + 1 : x_max(a) + t_right;
    
    % first look left, check for limit
    n = length(find(tx_left <= 0));
    if n > 0
        tx_left = tx_left(tx_left > 0);
        for b = 1:n
            tx_right(end + 1) = tx_right(end) + 1;
        end
    end
    
    % look right, check for limit
    n = length(find(tx_right > size(y_TOI, 2)));
    if n > 0
        tx_right = tx_right(tx_right <= size(y_TOI, 2));
        for b = 1:n
            tx_left(end + 1) = tx_left(1) - b;
        end
    end
    
    % append approved datapoints
    x_TOI_include_i = sort([x_TOI_include_i tx_left tx_right]);
    x_TOI_include = [x_TOI_include; x_TOI_include_i];    
end
x_include = x_TOI_include + start - 1;

% calculate the mean value
for a = 1:size(data_visual, 1)
    y_mean(a, 1) = mean(data_visual(a, x_include(a, :)));
end

% calculate peak latency
lat_peak = x_peak * xstep + xstart;
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
function fig = plot_box(data, datatype, condition, col, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % prepare data for a swarm
    data_group = [];
    data_y = [];
    for d = 1:size(data, 2)
        data_group(end + 1:end + size(data, 1), 1) = d;
        data_y = [data_y; data(:, d)];
    end
    
    % plot the markers
    beeswarm(data_group, data_y, 'colormap', col)
    
    % determine x limits
    xl = [0.25, size(data, 1)-0.25];
    xlim(xl)

    % add zero line
    line(xl, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.9) 

    % boxplot
    boxplot(data, condition, 'colors', col)

    % other parameters
    xlabel('TMS stimulus')
    set(gca, 'FontSize', 14) 
    set(gca, 'layer', 'top');    
    
    % y label
    if strcmp(datatype, 'amplitude')
        ylabel('amplitude (\muV)')
    elseif strcmp(datatype, 'GFP')
        ylabel('GFP (\muV)')
    elseif strcmp(datatype, 'SICI')
        ylabel('change in amplitude (\muV)')
    elseif strcmp(datatype, 'latency')
        ylabel('latency (ms)')
        view([90 -90])
    end
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
function cmap = blue2red(m)
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    cmap = [r g b]; 
end