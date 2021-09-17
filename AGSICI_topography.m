%% AG-SICI: ANALYSIS OF THE TOPOGRAPHY
% written by Dominika for AG-SICI project (2021)
% 
% 1) normalizes data 

%% parameters
clear all, clc

% datafile location
% data_path = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_path = 'C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_var = 'AGSICI_data';

% dataset
subject = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};

% load a header
% load('E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')
load('C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')

% visualization 
figure_counter = 8;
colours = [0 0.48 0.74; 0.99 0.18 0.18; 0.96 0.68 0.68];
lines = {':' '-'};
time_window = [-0.05, 0.3]; 
x = [time_window(1):header.xstep:time_window(2)];

% create output folders
foldername = 'AG-SICI_plus_figs';
folderpath = [pwd '\' foldername];
if ~exist(folderpath) 
    mkdir(folderpath)
end     

%% 1) data normalization 
% p = 1; s = 1; i = 1; c = 1;
% load the data
load(data_path, data_var)

% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)
            for i = 1:size(AGSICI_data, 6)
                data_temp(p, c, s, :, i) = mean(squeeze(AGSICI_data(p, c, :, s, :, i)), 1);
            end
        end
    end
end
clear p c s i

% average across current
for p = 1:length(position)
    for s = 1:length(subject)
        for i = 1:size(AGSICI_data, 6)
            data_temp_c(p, s, :, i) = mean(squeeze(data_temp(p, :, s, :, i)), 1);
        end
    end
end
data_temp_c = data_temp_c(:, :, 1:32, :);
clear p s i data_temp

% calculate GMFP for each subject/conil position
for p = 1:length(position)
    for s = 1:length(subject)    
        AGSICI_GMFP_subject(p, s, :) = std(squeeze(data_temp_c(p, s, :, :)), 1);
    end
end
clear p s 

% normalize data by GMFP
for p = 1:length(position)
    for s = 1:length(subject)  
        for i = 1:size(AGSICI_data, 6)
            % devide data at each time point by GMFP
            AGSICI_data_norm(p, s, :, i) = squeeze(data_temp_c(p, s, :, i)) / AGSICI_GMFP_subject(p, s, i);
        end
    end
end
clear p s i data_temp_c

% append new variables to the general MATLAB file
save(data_path, 'AGSICI_GMFP_subject', 'AGSICI_data_norm', '-append');

%% difference in coil position
% calculate DISS and spatial correlation C
AGSICI_DISS = struct;
for s = 1:length(subject)
    for i = 1:size(AGSICI_data_norm, 4)
        % between TS and CS stimuli
        diff = squeeze(AGSICI_data_norm(2, s, :, i) - AGSICI_data_norm(1, s, :, i));
        AGSICI_DISS.position(1, s, i) = sqrt(mean(diff.^2));
        AGSICI_DISS.position(2, s, i) = 1 - (AGSICI_DISS.position(1, s, i)^2)/2;
    end
end
clear s i diff

% plot DISS and correlation
data_visual = squeeze(mean(AGSICI_DISS.position(:, :, :), 2));
fig = figure(figure_counter);
hold on
title('Global dissimilarity - coil position', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(1, :) , 'facealpha', 0.5, 'linestyle', 'none');
plot(x, data_visual(1, :), 'Color', colours(2, :), 'LineWidth', 2.5)
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'spatial correlation' 'DISS index'}, 'Location', 'southeast');
lgd.FontSize = 12;
hold off
figure_name = 'AGSICI_DISS_position';
savefig([folderpath '\' figure_name '.fig'])
saveas(fig, [folderpath '\' figure_name '.png'])
figure_counter = figure_counter + 1;
clear data_visual fig ax lgd figure_name
