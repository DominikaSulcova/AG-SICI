%% AG SICI: TMS-EVOKED POTENTIALS
% Written by Dominika for AG-SICI project (2021)
% 
% preliminary version:
% 1) load the data
% 2) GFP
% 3) extract TEP amplitude

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
study = 'P2';
subject = [1:7];
protocol = {'spTMS', 'ppTMS'};
intensity = {'CS1', 'CS2', 'CS3'};
prefix = 'avg bl avgchan icfilt ica visual crop but fft-notchfilt prefilt prea'; 
filename = 'AG-SICI_P2';

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

% define conditions
conditions{1} = [protocol{1} ' TS'];
for p = 1:length(protocol)
    for i = 1:length(intensity)
        conditions{(p-1)*length(intensity) + i + 1} = [protocol{p} ' ' intensity{i}];
    end
end
clear p i

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        for a = 1:length(conditions)
           colours(a, :) = uisetcolor; 
        end
    case 'NO'
        if exist('P2_colours.mat') > 0
            load('P2_colours.mat')
        else
            disp('No colour scheme found in this directory!')    
        end
end
save('colours_P2.mat', 'colours'); 
clear a answer

% load a random header
load([prefix ' spTMS TS P2 01.lw6'], '-mat')
labels = {header.chanlocs.labels};
header.chanlocs(33) = [];

% visualization 
figure_counter = 1;
xstep = header.xstep;
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

% create folder for figures
folder_results = uigetdir(pwd, 'Choose the Results folder');
output_file = [folder_results '\' filename '.mat'];
folder_figures = [folder_results '\' filename '_figures'];
if ~exist(folder_figures) 
    mkdir(folder_figures);
end

%% 1) load the data - SAUTER --> LOAD DATA FROM MATLAB FILE
for s = 1:length(subject)
    for c = 1:length(conditions)
        if subject(s) < 10
            load([prefix ' ' conditions{c} ' ' study ' 0' num2str(subject(s)) '.mat'])
        else
            load([prefix ' ' conditions{c} ' ' study ' ' num2str(subject(s)) '.mat'])
        end
        A(s, c, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));
    end
end
disp(['Datasize: ' num2str(size(A))])
clear s c data 

% save dataset to the global MATLAB file
if exist(output_file) == 0
    AGSICI_P2_data = A;
    save(output_file, 'AGSICI_P2_data');
else
    load(output_file)
    if exist('AGSICI_P2_data') == 0
        AGSICI_P2_data = A;
        save(output_file, 'AGSICI_P2_data', '-append');
    else
        AGSICI_data = cat(1, AGSICI_P2_data, A);    
        save(output_file, 'AGSICI_P2_data', '-append');
    end
end
clear A

%% 2) GFP - SAUTER
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------
% compute individual GFP (exclude target channel)
AGSICI_P2_GFP = struct;
for s = 1:length(subject)
    for c = 1:length(conditions)
        AGSICI_P2_GFP.individual(s, c, :) = std(squeeze(AGSICI_P2_data(s, c, 1:32, :)), 1);  
    end
end
clear s c 

%% 3) TEP amplitude extraction
% ----- decide output parameters -----
AGSICI_P2_default.peak = {'N15' 'P25' 'N40' 'N75' 'N100' 'P180'};                                           % peaks of interest
AGSICI_P2_default.center = [0.015 0.025 0.04 0.075 0.1 0.2];                                                % default starting latencies
AGSICI_P2_default.span = [0.01 0.015 0.015 0.03 0.06 0.06];                                                            % span of the extraction window
AGSICI_P2_default.eoi = {{'target'} {'C1' 'C3' 'CP5'} {'target'} {'Cz'} {'Fz' 'FC2' 'F4'} {'Cz'}};          % electrodes of interest
percent = 20;                                                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                                                                                          % y limits for topoplots 
% ------------------------------------
% set colours for visualisation - figure 1
col_fig1 = [0    0.4471    0.7412;
    0.3922    0.8314    0.0745;
    0.3804    0.5882    0.1059;
    0.0275    0.4392    0.1647;
    0.9882    0.6353    0.6353;
    0.9882    0.2392    0.2392;
    0.6353    0.0784    0.1843];

% set colours for visualisation - figure 2
col_fig2 = [1.0000    0.4118    0.1608;
    0.8902    0.0549    0.0549;
    0.6353    0.0784    0.1843;
    0.4941    0.1843    0.5569;
    0.2941    0.2510    0.8706;
    0.0745    0.6235    1.0000];

% loop through subjects
for s = 2:length(subject)
    % setup name  
    figure_title = sprintf('Subject n. %d', subject(s)); 
    
    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    axis_counter = 1;
    fig = figure(figure_counter);
    
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    % choose data for topoplot
    for c = 1:length(conditions)                    
        for e = 1:32
            data_topoplot(c, e, 1, 1, 1, :) = squeeze(AGSICI_P2_data(s, c, e, :));
        end
    end
    
    % loop through peaks
    for k = 1:length(AGSICI_P2_default.peak)  
        % identify peak polarity
        if AGSICI_P2_default.peak{k}(1) == 'P'
            polarity = 1;
        else
            polarity = -1;
        end

        % identify EOIs
        eoi = [];
        for e = 1:length(AGSICI_P2_default.eoi{k})
            eoi(e) = find(strcmp(labels, AGSICI_P2_default.eoi{k}{e}));
        end

        % choose data
        data_visual = [];
        for c = 1:length(conditions)
            data_visual(c, :, :) = squeeze(AGSICI_P2_data(s, c, eoi, :)); 
        end

        % average across eois
        if length(eoi) > 1
            data_visual = squeeze(mean(data_visual, 2));
        end

        % define default TOI 
        center = AGSICI_P2_default.center(k);
        span = AGSICI_P2_default.span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;                
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 7, 1:21);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, AGSICI_P2_default.peak{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 7, 1:21)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % extract amplitude and latency
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, time_window(1), polarity);

             % update the figure
            subplot(4, 7, 1:21)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
                plot(x, data_visual(a, :), 'Color', col_fig1(a, :), 'LineWidth', 2.5)
                plot(lat_peak(a), y_max(a), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
                hold on
            end

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')

            % add the topoplot   
            for a = 1:size(data_visual, 1)
                subplot(4, 7, 21+a);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims) 
            end

            % ask for approval
            answer = questdlg('Do you want to proceed?', AGSICI_P2_default.peak{k},...
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
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / xstep);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / xstep);
                    choose_x = (choose_center - choose_span/2) : xstep : (choose_center + choose_span/2);

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
                    choose_figure_name = ['Choose manually peak ' AGSICI_P2_default.peak{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2)
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, 1)];
                    topo_plot(header, data_topoplot(1, :, :, :, :, :), choose_center, time_window(1), map_lims);
                        
                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    subplot(4, 4, 1) 
                    cla(choose_axesHandles(2))
                    topo_plot(header, data_topoplot(1, :, :, :, :, :), pos_x, time_window(1), map_lims);
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
        
        % record outcome variables
        for c = 1:length(conditions)            
            AGSICI_P2_TEP.amplitude_peak(s, c, k) = y_max(c); 
            AGSICI_P2_TEP.amplitude_mean(s, c, k) = y_mean(c); 
            AGSICI_P2_TEP.latency(s, c, k) = lat_peak(c); 
        end

        % set up the main figure
        figure(fig)
        n_row = length(AGSICI_P2_default.peak);
        n_col = 3 + length(AGSICI_P2_default.peak);
        subplot(n_row, n_col, [1 2 3] + n_col*(axis_counter - 1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', AGSICI_P2_default.peak{k}), 'Color', col_fig2(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig2(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = conditions;
        for a = 1:length(AGSICI_P2_default.peak)
            subplot(n_row, n_col, 3 + a + n_col*(axis_counter - 1))
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
    subplot(n_row, n_col, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if subject(s) < 10
        subj = ['0' num2str(subject(s))];
    else
        subj = num2str(subject(s));
    end
    figure_name = ['AGSICI_P2_TEP_amp_' subj];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])
    close(fig)
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'AGSICI_P2_TEP', '-append');

    % update the figure counter
    figure_counter = figure_counter + 1;  
end
clear finish fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
    choose data choose_center choose_axesHandles answer col_fig1 col_fig2 s c figure_title fig axis_counter k choose_data choose_figure_name...
    polarity e i data_visual data_topoplot a subj figure_name n_col n_row eoi lat_peak span center yl yl_sgtitle y_mean y_max p 

% append new variables to the general MATLAB file
save(output_file, 'AGSICI_P2_default', '-append');
clear percent map_lims

%% 4) TEP amplitude - visualization
% ----- decide output parameters -----
index = logical([0 1 0 0 1 1 1;
    1 1 1 1 1 1 0;
    1 1 1 1 1 1 1;
    0 1 1 1 1 1 1;
    1 1 1 1 1 1 1;
    1 1 1 1 1 1 1]);                                                                                        % y limits for topoplots 
% ------------------------------------
% plot mean spTMS amplitudes
for p = 1:size(AGSICI_P2_TEP.amplitude_peak, 3)     
    %  get mean data
    data_visual = []; sem_visual = [];
    for c = 1:4
        data_visual(c) = mean(AGSICI_P2_TEP.amplitude_peak(index(p, :), c, p));
        sem_visual(c) = std(AGSICI_P2_TEP.amplitude_peak(index(p, :), c, p))/sqrt(length(subject(index(p, :))));
    end

    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % plot the amplitudes
    for c = 1:4
        bar(c, data_visual(c), 'EdgeColor', 'none', 'FaceColor', colours(c, :))
    end
    
    % errorbars
    errorbar(1:4, data_visual, sem_visual, sem_visual, 'k', 'linestyle', 'none');  
    
    % add parameters
    title(sprintf('peak amplitude: %s', AGSICI_P2_default.peak{p}))
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    
    % name and save figure
    figure_name = ['AGSICI_TEP_'  AGSICI_P2_default.peak{p}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])
    
    % upgrade counter
    figure_counter = figure_counter + 1;
end
clear p c data_visual sem_visual figure_name

%% 5) SICI: TEPs
% ----- decide output parameters -----
electrode = {'target' 'Cz'};
linestyle = {'-' ':'};
% ------------------------------------
% calculate individual SICI
AGSICI_P2_SICI = struct;
for s = 1:length(subject)
    for c = 5:7
        for e = 1:size(AGSICI_P2_data, 3)
            % calculate GFP (exclude target channel)
            AGSICI_P2_SICI.individual(s, c-4, e, :) = squeeze(AGSICI_P2_data(s, c, e, :)) - squeeze(AGSICI_P2_data(s, 1, e, :));  
        end
    end
end
clear s c e

% calculate mean SICI
for c = 1:3
    for e = 1:size(AGSICI_P2_data, 3)
        for i = 1:size(AGSICI_P2_data, 4)
            AGSICI_P2_SICI.mean(c, e, i) = mean(squeeze(AGSICI_P2_SICI.individual(:, c, e, i)));
            AGSICI_P2_SICI.CI(c, e, i) = (std(squeeze(AGSICI_P2_SICI.individual(:, c, e, i)))/sqrt(length(subject))) * z;
        end
    end  
end 
clear c e i 

% plot mean SICI 
for c = 1:3
    % prepare data
    data_visual = squeeze(AGSICI_P2_SICI.mean(c, :, :));
    
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    yl = [-2.5 2.5];

    % shade interpolated interval
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [-0.005, yl(1), 0.015, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
    
    % plot all channels
    for e = 1:size(data_visual, 1)
        plot(x, data_visual(e, :), 'Color', colours(c, :))
    end
    
    % plot electrode(s) of interest
    for e = 1:length(electrode)
        e_n = find(strcmp(labels, electrode{e}));
        e_p(e) = plot(x, data_visual(e_n, :), 'Color', [0 0 0], 'LineWidth', 2.5, 'LineStyle', linestyle{e})
    end

    % mark TMS stimulus and zerol line
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    line(time_window, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1.5)

    % add other parameters
    title(sprintf('SICI in AG: %s', conditions{4+c}))
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(yl)
    
    % add legend
    lgd = legend(e_p, electrode, 'Location', 'southeast');
    lgd.FontSize = 14;
    hold off
 
    % change figure size
    fig.Position = [250 250 600 350];

    % name and save figure
    figure_name = ['AGSICI_TEP_SICI_'  conditions{4+c}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear c data_visual fig yl lgd figure_name 

% save SICI datasets for letswave
for c = 1:3
    % create data
    for s = 1:length(subject)
        data(s, :, 1, 1, 1, :) = squeeze(AGSICI_P2_SICI.individual(s, c, :, :));
    end

    % create header
    header.name = ['merged SICI ' conditions{4+c}];
    header.datasize = size(data);
    header.xstart = time_window(1);

    % save
    save([header.name '.mat'], 'data')
    save([header.name '.lw6'], 'header')
end
clear c s

% calculate GFP of SICI
for s = 1:length(subject)
    for c = 1:3
        % calculate GFP (exclude target channel)
        AGSICI_P2_SICI.GFP(s, c, :) = std(squeeze(AGSICI_P2_SICI.individual(s, c, 1:32, :)), 1);  
    end
end
clear s c 
clear electrode linestyle

% append to the general MATLAB file
save(output_file, 'AGSICI_P2_SICI', '-append');

%% 6) SICI: change in the amplitude
% calculate individual change for peak amplitude of each TEP peak 
for p = 1:length(AGSICI_P2_default.peak) 
    for s = subject(index(p, :))
        for c = 1:3        
            AGSICI_P2_SICI.peak_change(s, c, p) = AGSICI_P2_TEP.amplitude_peak(s, 4+c, p) - AGSICI_P2_TEP.amplitude_peak(s, 1, p);
        end
    end
end
clear s c p

% append to the general MATLAB file
save(output_file, 'AGSICI_P2_SICI', '-append');

% get data for visualization
data_visual = []; sem_visual = [];
for p = 1:length(AGSICI_P2_default.peak)
    for c = 1:3
        % identify peak polarity, get the data
        if AGSICI_P2_default.peak{p}(1) == 'P'
            polarity = 1;
        else
            polarity = -1;
        end
        
        %
        data_visual(p, c) = mean(AGSICI_P2_SICI.peak_change(subject(index(p, :)), c, p)) * polarity;
        sem_visual(p, c) = std(AGSICI_P2_SICI.peak_change(subject(index(p, :)), c, p))/sqrt(length(subject(index(p, :))));
    end
end
clear p c polarity     

% launch the figure
fig = figure(figure_counter);
hold on

% plot errorbars
ngroups = size(data_visual, 1);
nbars = size(data_visual, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x_bar = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x_bar, data_visual(:, i), sem_visual(:, i), sem_visual(:, i), 'k', 'linestyle', 'none');
end 

% plot the data
barplot = bar(data_visual, 'EdgeColor', 'none');
col = colours([2:4], :);
for b = 1:size(data_visual, 2)
    barplot(b).FaceColor = col(b, :);
end

% set other parameters
title('AG SICI: change in TEP peak amplitude')
ylabel('amplitude (\muV \pmSEM)');
xlabel('TEP component')
set(gca, 'xtick', 1:6, 'xticklabel', AGSICI_P2_default.peak)
set(gca, 'Fontsize', 14)
legend(barplot, intensity, 'Location', 'southeast', 'fontsize', 14)

% name and save figure
figure_name = 'SICI_amplitude';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear data_visual sem_visual fig i nbars ngroups x_bar b barplot col figure_name

%% functions
function peak_x = gfp_plot(x, y, time_window, xstep, labeled, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% launch the figure  
plot(x, y)
yl = get(gca, 'ylim');
cla

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data, mark TMS stimulus
plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.015);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.015
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.220
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
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
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

