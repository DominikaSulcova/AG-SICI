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
%       - plot GFP alone --> for pub figure
% 
% 4) IDENTIFY EOIs


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
% --------------------------------

% define conditions
condition{1} = [protocol{1} ' TS'];
for p = 1:length(protocol)
    for i = 1:length(intensity)
        condition{(p-1)*length(intensity) + i + 1} = [protocol{p} ' ' intensity{i}];
    end
end
clear p i

% load a header
load([prefix ' ' study ' 01 ' condition{1} '.lw6'], '-mat')

% visualization parameters
figure_counter = 1;
shade = 0.2;
electrodes = {header.chanlocs.labels};
x_step = header.xstep;
x_datastart = header.xstart;
x = [time_window(1):x_step:time_window(2)];
x_start = (time_window(1) - x_datastart)/x_step;
x_end = (time_window(2) - x_datastart)/x_step;

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
outliers = [8, 11, 15];
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
    savelw(data, header, name)        
end
clear c s e t data name

% compute average and SEM
data_mean = squeeze(mean(data_individual, 2));
data_SEM =  squeeze(std(data_individual, 0, 2)/sqrt(size(data_individual, 2)));

% append to the outcome MATLAB file
save(output_file, 'data_MS', 'data_individual', '-append');
clear target 

%% 2) TEP BUTTERFLY - individual conditions
% ----- section input -----
channel = 'Cz';
y_lim = {[-4 4] [0 2.2]};
% -------------------------
% identify channel to highlight
for e = 1:length(electrodes)
    if strcmp(electrodes{e}, channel)
        channel_n = e;
    end
end
clear e

% colour scheme
colours_all = [0.65, 0.65, 0.65; colours(2:4, :); colours(2:4, :)];

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

%% 3) GFP - all conditions
% ----- section input -----
labeled = 'off';
max_peaks = 6;
y_lim = [0 2.5];
% -------------------------
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

clear labeled max_peaks data gfp h_axis fig figure_name pos

%% 4) IDENTIFY EOIs
% launch the default structure
eoi_n = [3,3,3,3,1];
AGSICI_TEP_default = struct;
AGSICI_TEP_default.peak = {'P25' 'P30' 'N75' 'N100' 'P180'};
AGSICI_TEP_default.center = [22.5 31 72.5 116 187.5];
AGSICI_TEP_default.span = [14 30 30 40 60];

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

%% functions
function savelw(data, header, name)
    % modify header
    header.name = name;
    header.datasize = size(data);
    header.xstart = -0.05;
    header.chanlocs(length(header.chanlocs)+1).labels = 'target';
    for e = 1:size(data, 1)
        header.events(e).code = 'AG_TEP';
        header.events(e).latency = 0;
        header.events(e).epoch = e;
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


