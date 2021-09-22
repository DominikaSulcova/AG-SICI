%% AG-SICI: ANALYSIS OF THE TOPOGRAPHY
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
% data_path = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_path = 'C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_var = 'AGSICI_data';
load(data_path, data_var)

% load a header
% load('E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')
load('C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')

% visualization 
figure_counter = 1;
colours = [0 0.48 0.74; 0.99 0.18 0.18; 0.96 0.68 0.68];
lines = {':' '-'};
time_window = [-0.05, 0.3]; 
x = [time_window(1):header.xstep:time_window(2)];

% create output folders
folder_fig = [pwd '\AG-SICI_plus_figs'];
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

%% 2) group-level GFP
% ----- adjustable parameters -----
TOI = [0.015, 0.100];
% ----- adjustable parameters -----

% average data across subjects and intensities
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            % average across subjects
            for e = 1:size(AGSICI_data, 5)
                for t = 1:size(AGSICI_data, 6)
                    data_mean(p, c, i, e, t) = mean(squeeze(AGSICI_data(p, c, i, :, e, t)));
                end
            end 
        end        
        % average across intensities
        for e = 1:size(AGSICI_data, 5)
            for t = 1:size(AGSICI_data, 6)
                data_gmfp(p, c, e, t) = mean(squeeze(data_mean(p, c, :, e, t)));
            end
        end 
    end
end
clear p c i e t 

%calculate GMFP and plot separate figures
for p = 1:length(position)
    for c = 1:length(current)
        % calculate GMFP (exclude target channel)
        GFP(p, c, :) = std(squeeze(data_gmfp(p, c, 1:size(data_gmfp, 3) - 1, :)), 1); 
        
        % launch the figure
        fig = figure(figure_counter);
        title(sprintf('%s STS - %s current', position{p}, current{c}), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'fontsize', 12)
        xlabel('time (s)')
        ylabel('GFP (\muV)')
        hold on

        % set limits of the figure
        plot(x, squeeze(GFP(p, c, :)))
        yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
        xl = get(gca, 'xlim');
        cla
        ylim(yl); xlim(xl);
        
        % define colour
        col = colours_cond((p-1)*2 + c, :);
        
        % mark the area of interest
        x_start = (TOI(1) - time_window(1))/header.xstep;
        x_end = (TOI(2) - time_window(1))/header.xstep;
        x_fill = [TOI(1) : header.xstep: TOI(2)];
        fill([x_fill fliplr(x_fill)],[squeeze(GFP(p, c, x_start:x_end))' zeros(1, length(squeeze(GFP(p, c, x_start:x_end))))], ...
            col , 'facealpha', 0.3, 'linestyle', 'none');
        
        % plot interpolated part
        rectangle('Position', [0, yl(1)+0.01, 0.01, yl(2) - yl(1)], ...
            'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none') 
        
        % plot signal 
        P = plot(x, squeeze(GFP(p, c, :)), 'Color', col, 'LineWidth', 2.5);
        
        % plot lines
        line(xl, [0 0], 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.5)
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
        hold off
        
        % update position
        pos = get(gcf, 'Position');
        fig.Position = [pos(1) pos(2)-200 600 300];
        
        % save figure
        figure_name = sprintf('AGSICI_GFP_%s_%s', position{p}, current{c});
        savefig([folder_fig '\' figure_name '.fig'])
        saveas(fig, [folder_fig '\' figure_name '.png'])
        
        % update counter
        figure_counter = figure_counter + 1;
    end
end
clear p c fig zl xl x_start x_end x_fill col P figure_name pos

%% ) data normalization - all levels of ORIENTATION
% p = 1; s = 1; i = 1; c = 1;
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

% calculate GMFP for each subject/conil position
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)    
            AGSICI_GMFP_subject(p, c, s, :) = std(squeeze(data_temp_i(p, c, s, :, :)), 1);
        end
    end
end
clear p s c

% normalize data by GMFP
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)  
            for i = 1:size(AGSICI_data, 6)
                % devide data at each time point by GMFP
                AGSICI_data_norm(p, c, s, :, i) = squeeze(data_temp_i(p, c, s, :, i)) / AGSICI_GMFP_subject(p, c, s, i);
            end
        end
    end
end
clear p c s i 

% append new variables to the general MATLAB file
save(data_path, 'AGSICI_GMFP_subject', 'AGSICI_data_norm', '-append');

%% 3) data normalization - pool current conditions

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

% average across current
for p = 1:length(position)
    for s = 1:length(subject)
        for i = 1:size(AGSICI_data, 6)
            data_temp_c(p, s, :, i) = mean(squeeze(data_temp_i(p, :, s, :, i)), 1);
        end
    end
end
data_temp_c = data_temp_c(:, :, 1:32, :);
clear p s i data_temp_i

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
clear p s i d

% append new variables to the general MATLAB file
save(data_path, 'AGSICI_GMFP_subject', 'AGSICI_data_norm', '-append');

%% 4) plot 

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
