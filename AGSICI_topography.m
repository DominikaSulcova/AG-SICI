%% AG-SICI: ANALYSIS OF RESPONSE TOPOGRAPHY
% written by Dominika for AG-SICI project (2021)
% 
% 1) normalizes data 

%% parameters
clear all, clc

% dataset
subject = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};

% load data
data_path = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
% data_path = 'C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_var = 'AGSICI_data';
load(data_path, data_var)

% load a header
load('E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')
% load('C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')

% visualization 
figure_counter = 1;
colours = [0.902 0.235 0.235; 0.075 0.624 1; 0.922 0.471 0.831; 0.467 0.675 0.188];
lines = {':' '--' '-'};
time_window = [-0.05, 0.3]; 
x = [time_window(1):header.xstep:time_window(2)];

% create output folders
folder_fig = [uigetdir '\AG-SICI_plus_figs'];
if ~exist(folder_fig) 
    mkdir(folder_fig)
end    
folder_exp = [pwd '\AG-SICI_export'];
if ~exist(folder_exp) 
    mkdir(folder_exp)
end

%% 1) export data for Ragu
% ----- adjustable parameters -----
folder_name = {'AG-SICI_all' 'AG-SICI_orientation'};
% ----- adjustable parameters -----

% create output folder
for a = 1:length(folder_name)
    mkdir([folder_exp '\' folder_name{a}])
end

% write text files for Ragu --> all conditions
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % choose data to write, remove 'target' channel
                data = squeeze(AGSICI_data(p, c, i, s, 1:32, :))';
                
                % define subject
                if subject(s) < 10
                    subj = ['S0' num2str(subject(s))];
                else
                    subj = ['S' num2str(subject(s))];
                end
                
                % save as .csv               
                filename = ['AGSICI_' subj '_' position{p} '_' current{c} '_' intensity{i}(end-2:end) '.csv']; 
                writematrix(data, [folder_exp '\' folder_name{1} '\' filename])
            end
        end
    end
end
clear p c i s data subj filename     

% write text files for Ragu --> average over intensities
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)
            % choose data to write, remove 'target' channel
            data = squeeze(mean(AGSICI_data(p, c, :, s, 1:32, :), 3))';

            % define subject
            if subject(s) < 10
                subj = ['S0' num2str(subject(s))];
            else
                subj = ['S' num2str(subject(s))];
            end

            % save as .csv               
            filename = ['AGSICI_' subj '_' position{p} '_' current{c} '.csv']; 
            writematrix(data, [folder_exp '\' folder_name{2} '\' filename])
        end
    end
end
clear p c s data subj filename     

% create the montage file
filename = [folder_exp '\AGSICI_montage.xyz'];
fileID = fopen(filename, 'a');
fprintf(fileID, '32\r\n');
for a = 1:32
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear filename fileID a folder_name

%% 3) DISS - between levels of ORIENTATION
% ----- adjustable parameters -----
comp = {[1 1; 1 2] [2 1; 2 2] [1 1; 2 1] [1 2; 2 2]};
% ----- adjustable parameters -----

% calculate DISS and spatial correlation C  
AGSICI_DISS = struct;
for c = 1:length(comp)
    AGSICI_DISS(c).comparison = [position{comp{c}(1, 1)} '-' current{comp{c}(1, 2)} ' X ' position{comp{c}(2, 1)} '-' current{comp{c}(2, 2)}];
    for s = 1:length(subject)
        for i = 1:size(AGSICI_data_norm, 5)
            % calculate difference in normalized data
            diff = squeeze(AGSICI_data_norm(comp{c}(1, 1), comp{c}(1, 2), s, :, i) - AGSICI_data_norm(comp{c}(2, 1), comp{c}(2, 2), s, :, i));
            
            % calculate outcome values
            AGSICI_DISS(c).data.DISS(s, i) = sqrt(mean(diff.^2));
            AGSICI_DISS(c).data.corr(s, i) = 1 - (AGSICI_DISS(c).data.DISS(s, i)^2)/2;
        end
    end
end
clear c s i diff

% plot the figures
for c = 1:length(comp)
    % choose data 
    data_visual(1, :) = mean(AGSICI_DISS(c).data.DISS, 1);
    data_visual(2, :) = mean(AGSICI_DISS(c).data.corr, 1);

    % plot DISS and correlation
    fig = figure(figure_counter);
    plot_DISS(x, data_visual, 'legend', {'DISS' 'spatial correlation'})
    title(['DISS: ' AGSICI_DISS(c).comparison], 'FontSize', 16, 'FontWeight', 'bold')

    % save the figure
    figure_name = ['AGSICI_DISS_' AGSICI_DISS(c).comparison];
    saveas(fig, [folder_fig '\' figure_name '.png'])
    figure_counter = figure_counter + 1;
end
clear c data_visual fig figure_name

% append new variables to the general MATLAB file
save(data_path, 'AGSICI_DISS', '-append');

%% 2) data normalization - all levels of ORIENTATION
% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)
            for i = 1:size(AGSICI_data, 6)
                data_temp_i(p, c, s, :, i) = mean(squeeze(AGSICI_data(p, c, :, s, :, i)), 1);
            end
        end
    end
end
clear p c s i

% calculate GFP for each subject/conil position
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)    
            AGSICI_GFP_subject(p, c, s, :) = std(squeeze(data_temp_i(p, c, s, :, :)), 1);
        end
    end
end
clear p s c

% normalize data by GFP
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)  
            for i = 1:size(AGSICI_data, 6)
                % devide data at each time point by GMFP
                AGSICI_data_norm(p, c, s, :, i) = squeeze(data_temp_i(p, c, s, :, i)) / AGSICI_GFP_subject(p, c, s, i);
            end
        end
    end
end
clear p c s i 

% append new variables to the general MATLAB file
save(data_path, 'AGSICI_GFP_subject', 'AGSICI_data_norm', '-append');

%% 4) plot ICA - average timecourse
% ----- adjustable parameters -----
folder_ica = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\ICA2';
% ----- average across all categories -----

% prepare data
load('avg merged avg ica-um ica2.mat')
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;
data_visual = squeeze(data(1, :, 1, 1, 1, x_start:x_end));
clear data x_start x_end

% plot overall average ICA timecourse
for a = 1:size(data_visual, 1)
    % plot tthe figure
    fig = figure(figure_counter)
    plot_TEP(x, data_visual(a, :), 'peak_latency', true)
    
    figure_name = ['IC' num2str(a) '_course'];
    saveas(fig, [folder_ica '\' figure_name '.png'])
    figure_counter = figure_counter + 1;
end
clear a fig figure_name data_visual

% ----- differences in orientation -----
% prepare data
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;
for p = 1:length(position)
    for c = 1:length(current)
        load(['avg merged avg ica-um ica2 ' position{p} ' ' current{c} '.mat'])
        data_visual((p-1)*2 + c, :, :) = squeeze(data(1, :, 1, 1, 1, x_start:x_end));
    end
end
clear p c data x_start x_end

% plot overall average ICA timecourse
for a = 1:size(data_visual, 2)
    % plot the figure
    fig = figure(figure_counter)
    plot_TEP(x, squeeze(data_visual(:, a, :)), ...
        'colours', colours, 'legend', {'along normal' 'along reversed' 'across normal' 'across reversed'})
    
    figure_name = ['IC' num2str(a) '_course_orient'];
    saveas(fig, [folder_ica '\' figure_name '.png'])
    figure_counter = figure_counter + 1;
end
clear a fig figure_name data_visual

% ----- differences in intensity -----
% prepare data
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;
for i = 1:length(intensity)
    load(['avg merged avg ica-um ica2 ' intensity{i} '.mat'])
    data_visual(i, :, :) = squeeze(data(1, :, 1, 1, 1, x_start:x_end));
end
clear i data x_start x_end

% plot overall average ICA timecourse
for a = 1:size(data_visual, 2)
    % plot the figure
    fig = figure(figure_counter)
    plot_TEP(x, squeeze(data_visual(:, a, :)), 'legend', {'100 %rMT' '120 %rMT' '140 %rMT'})
    
    figure_name = ['IC' num2str(a) '_course_int'];
    saveas(fig, [folder_ica '\' figure_name '.png'])
    figure_counter = figure_counter + 1;
end
clear a fig figure_name data_visual

%% 5) hilbert
% ----- adjustable parameters -----
prefix = 'avg hilbert highpass_35Hz bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1';
% ----- average across all categories -----

% extract the data
AGSICI_hilbert = struct;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % define subject
                if subject(s) < 10
                    subj = ['0' num2str(subject(s))];
                else
                    subj = [num2str(subject(s))];
                end
                
                % load the data
                name = [prefix ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
                load(name)
                
                % append to the outcome structure
                AGSICI_hilbert.data(p, c, i, s, :, :) = squeeze(data(:, 1:32, :, :, :, :));
            end
        end
    end
end
clear p c i s subj name data

% calculate RMS
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % for each electrode separately
                for e = 1:size(AGSICI_hilbert.data, 5)
                    AGSICI_hilbert.RMS(p, c, i, s, e) =  sqrt(mean(squeeze(AGSICI_hilbert.data(p, c, i, s, e, :)).^2));
                end
                
                % average across electrodes
                AGSICI_hilbert.RMS_avg(p, c, i, s) = mean(AGSICI_hilbert.RMS(p, c, i, s, :));
            end
        end
    end
end
clear p c i s e   

% calculate mean values to plot
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_visual(i, (p-1)*2 + c) = mean(AGSICI_hilbert.RMS_avg(p, c, i, :));
            data_sem(i, (p-1)*2 + c) = std(AGSICI_hilbert.RMS_avg(p, c, i, :))/sqrt(length(subject));
        end
    end
end
clear p c i 

% plot average RMS by condition
fig = figure(figure_counter)
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for a = 1:size(data_visual, 1)
    barplot(a).FaceColor = colours(a, :)
end
legend(barplot, {'along - normal', 'along - reversed', 'across - normal', 'across - reversed'}, ...
    'Location', 'bestoutside', 'fontsize', 14)
xlabel('intensity of stimulation (%rMT)')
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
ylim([0, 1.6])
hold off

%% hilbert - correlation with N45 amplitude
% prepare data
load(data_path, 'AGSICI_outcome')
data_cor = [];
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_cor(end+1 : end+length(subject), 1) = AGSICI_hilbert.RMS_avg(p, c, i, :);
        end
    end
end
for a = 1:height(AGSICI_outcome)
    if strcmp(AGSICI_outcome.peak{a}, 'N45') 
        data_temp(a) = AGSICI_outcome.amplitude(a);
    end
end

% prepare the model
data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

% choose only correlations that show TEP-MEP interactions
if col(a) == 6
% plot data + regression line
fig = figure(figure_counter);
hold on
plot_cor = plotAdded(data_model);

% adjust parameters    
title([medication{m} ' : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
set(gca, 'FontSize', 14)
plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
legend off
if data_model.Coefficients.Estimate(2) > 0
    text_pos = [0.95 0.85 0.75];
else
    text_pos = [0.25 0.15 0.05];
end
T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.3f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic');   
set(T(2), 'fontsize', 14); 
set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 
hold off

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
function plot_DISS(x, data_visual, varargin)
% check for colours
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'colours'));
    if ~isempty(a)
        col = varargin{a + 1};
    else
        col = [0.24 0.49 0.99; 0.87 0.16 0.40];
    end
else
    col = [0.24 0.49 0.99; 0.87 0.16 0.40];
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

% plot spatial correlation
hold on
P(1) = fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    col(1, :) , 'facealpha', 0.5, 'linestyle', 'none');

% plot DISS
P(2) = plot(x, data_visual(1, :), 'Color', col(2, :), 'LineWidth', 3)

% add lines
line([x(1) x(end)], [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line([x(1) x(end)], [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 4) 

% plot legend if required
if length(legend_on) > 0
    lgd = legend(P([2 1]), legend_on, 'Location', 'southeast');
    lgd.FontSize = 14;
end

% set up the parameters
set(gca, 'fontsize', 14); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
ylim([-0.75 1.75]) 
xlim([x(1) x(end)])
end
