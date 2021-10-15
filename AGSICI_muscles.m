%% AG-SICI: MUSCULAR CONTRIBUTION TO TEP COMPONENT N45
% Written by Dominika for GABA-AD project (2021)
% 
% ----- muscle contraction X N45 -----
% 1) prepares the data
%    - loads preprocessed data from letswave
%       - raw datasets to identify early muscular contraction  
%       - GFP of processed datasets to identify N45 
%    - crops the raw data in predefined time window
%       - default [-0.05, 0.3]s
%   - saves to the output structure --> 'AGSICI_muscle_activity'
% 
% 2) calculates GFP
%   - GFP of individual datasets
%   - mean GFP per condition
%   - exports in letswave format for browsing
%   - saves to the output structure and appends to global matlab file
% 
% 3) plots mean GFP 
%   - GFP of muscle contraction and N45 in one figure
%   - shades TOIs of interest
%   - saves as .fig and .png to the output folder
% 
% 4) 
% 
% 5)
% 
% 6) exports for R 
%   - saves extracted values in a long format .csv table 'AGSICI_cont_X_N45'
% 
% 7) plots mean values
%   - separate line plots for muscle contraction and N45
%   - saves as .fig and .png to the output folder
% 
% 8) overall correlation
%   -

%% parameters
clear all; clc

% adjustable parameters
prompt = {'Include subjects:' 'Crop the data:'};
dlgtitle = 'Adjustable parameters';
dims = [1 50];
definput = {'[1, 3:18, 20, 21]' '[-0.05, 0.3]'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
eval(['subject = ' answer{1} ';']); 
eval(['time_window = ' answer{2} ';']);  
clear prompt dlgtitle dims definput answer

% dataset
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

% create output folder
output_folder = [uigetdir(pwd, 'Choose the directory with the global output file') '\AG-SICI_muscles'];
if ~exist(output_folder) 
    mkdir(pwd, 'AG-SICI_muscles')
end    

% identify matlab file with results
results_file = [uigetdir(pwd, 'Choose the directory with the global output file') '\AG-SICI_plus.mat'];

% load template header
load(results_file, 'header')

% visualization 
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

%% 1) muscle contraction X N45: prepare individual GFP
% ----- adjustable parameters -----
prefix_m = 'crop avg art-sup raw P1'; 
electrodes = 1:2:31;
% ----- adjustable parameters -----

% load GFP of cleaned datasets
load(results_file, 'AGSICI_GFP_subject')

% create structure with individual GFP
AGSICI_muscle_activity = struct;
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
                
                % load and crop raw data --> muscular contraction 
                name = [prefix_m ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
                load(name)
                AGSICI_muscle_activity.contraction.data(p, c, i, s, :, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));  
                
                % calculate GFP of raw data - only selected electrodes
                AGSICI_muscle_activity.contraction.GFP(p, c, i, s, :) = std(squeeze(AGSICI_muscle_activity.contraction.data(p, c, i, s, electrodes, :)), 1);
                
                % fill in GFP of processed data --> N45
                AGSICI_muscle_activity.N45.GFP(p, c, i, s, :) = AGSICI_GFP_subject.data(p, c, i, s, :);  
            end 
        end
    end
end
clear p c i s name prefix_m subj AGSICI_GFP_subject electrodes

% append new variables to the general MATLAB file
save(results_file, 'AGSICI_muscle_activity', '-append');

%% 3) muscle contraction X N45: calculate AUC  
% ----- adjustable parameters -----
TOI_m = [0.003, 0.01];
TOI_N45 = [0.035, 0.055];
% ----- adjustable parameters -----

% set cropping limits for TOIs
x_start_m = (TOI_m(1) - time_window(1))/xstep;
x_end_m = (TOI_m(2) - time_window(1))/xstep;
x_start_N45 = (TOI_N45(1) - time_window(1))/xstep;
x_end_N45 = (TOI_N45(2) - time_window(1))/xstep;

% calculate AUC over TOI for each subject/condition 
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)  
                AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, s) = trapz(squeeze(AGSICI_muscle_activity.contraction.GFP(p, c, i, s, x_start_m:x_end_m)));
                AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, s) = trapz(squeeze(AGSICI_muscle_activity.N45.GFP(p, c, i, s, x_start_N45:x_end_N45)));
            end
        end
    end
end
clear p c i s x_start_m x_end_m x_start_N45 x_end_N45

% append new variables to the general MATLAB file
save(results_file, 'AGSICI_muscle_activity', '-append');

%% 4) muscle contraction X N45: export AUC for R 
% write in a long-format table
AGSICI_cont_X_N45 = table;
AGSICI_cont_X_N45.subject = zeros(0); AGSICI_cont_X_N45.orientation = {}; AGSICI_cont_X_N45.intensity = zeros(0); 
AGSICI_cont_X_N45.contraction = zeros(0); AGSICI_cont_X_N45.N45 = zeros(0); 
row_cnt = 1;
for s = 1:length(subject)
    for p = 1:length(position)
        for c = 1:length(current)
            for i = 1:length(intensity)
                % fill in the row
                AGSICI_cont_X_N45.subject(row_cnt) = subject(s);             
                AGSICI_cont_X_N45.orientation(row_cnt) = {[position{p} '-' current{c}]};
                AGSICI_cont_X_N45.intensity(row_cnt) = str2double(intensity{i}(end-2:end));
                AGSICI_cont_X_N45.contraction(row_cnt) = AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, s);
                AGSICI_cont_X_N45.N45(row_cnt) = AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, s);
                
                % update row counter
                row_cnt = row_cnt + 1;
            end
        end
    end
end

% save as .csv
writetable(AGSICI_cont_X_N45, [output_folder '\AGSICI_cont_X_N45.csv'])
clear p c s i row_cnt

%% 5) muscle contraction X N45: plot mean values - all subjects 
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
            y(i) = mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, :));
            SEM(i) = std(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, :)) / sqrt(length(subject));
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
title(sprintf('Muscular activity (GFP - AUC): %d - %dms', TOI_m(1) * 1000, TOI_m(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('AUC (\muV*s \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);

% save the figure       
savefig([output_folder '\AGSICI_cont_all.fig'])
saveas(fig, [pwd '\AGSICI_cont_all.png'])

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
            y(i) = mean(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, :));
            SEM(i) = std(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, :)) / sqrt(length(subject));
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
title(sprintf('Component N45 (GFP - AUC): %d - %dms', TOI_N45(1) * 1000, TOI_N45(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('AUC (\muV*s \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);

% save the figure       
savefig([output_folder '\AGSICI_N45_all.fig'])
saveas(fig, [pwd '\AGSICI_N45_all.png'])

% update the counter
figure_counter = figure_counter + 1;  
clear fig leg perr y SEM

%% 6) muscle contraction: distribution across conditions - all subjects 
% calculate mean values to plot
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_visual(i, (p-1)*2 + c) = mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, :));
        end
    end
end
clear p c i 

% plot average GFP-AUC by condition
fig = figure(figure_counter)
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along - normal', 'along - reversed', 'across - normal', 'across - reversed'}, ...
    'Location', 'bestoutside', 'fontsize', 14)
xlabel('intensity of stimulation (%rMT)')
ylabel('muscular contraction')
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(0.75, yl(2) - 250, 'all subjects', 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_size_all.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size_all.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig yl barplot a data_visual 

% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        data_visual(:, (p-1)*2 + c) = squeeze(mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, :, :), 3));
    end
end
clear p c

% plot average RMS - boxplot
fig = figure(figure_counter); ax = gca;       
hold on
boxplot(data_visual, 'color', colours)
ax.XTickLabel = '';
label_array = {'along' 'along' 'across' 'across'; 'normal' 'reversed' 'normal' 'reversed'}; 
for i = 1:length(label_array)
    text(i, ax.YLim(1), sprintf('%s\n%s', label_array{:, i}), 'FontSize', 14, ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');    
end
ylabel('muscular contraction - GFP-AUC')
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(2, yl(2) - 1500, 'all subjects', 'FontSize', 14)
clear i label_array yl

% plot the markers
for b = 1:4
    scat(b) = scatter(repelem(b, size(data_visual, 1)), data_visual(:, b),...
        75, colours(b, :), 'filled');
end

% mark outliers
h_out = flipud(findobj(gcf,'tag','Outliers'));
for h = 1:length(h_out)
    x_out =  get(h_out(h), 'XData');
    y_out =  get(h_out(h), 'YData');
    for i = 1:length(x_out)
        if ~(isnan(x_out(i)))
            index_out(h, i) = find(data_visual(:, h) == y_out(i));
            text(x_out(i) + 0.1, double(y_out(i)), sprintf('%d', subject(index_out(h, i))))
        end
    end
end
clear h_out h x_out y_out i index_out

% save the figure
savefig([output_folder '\AGSICI_cont_size.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig ax scat b data_visual 

%% 7) muscle contraction X N45: calculate mean GFP, save for LW
% ----- adjustable parameters -----
outliers = [5 8 12];
% ----- adjustable parameters -----

% calculate index to exclude outliers
index = subject ~= outliers(1) & subject ~= outliers(2) & subject ~= outliers(3);

% calculate mean GFP per condition
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for t = 1:size(AGSICI_muscle_activity.contraction.GFP, 5)    
                % muscular contraction 
                AGSICI_muscle_activity.contraction.GFP_mean(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.contraction.GFP(p, c, i, :, t)));
                AGSICI_muscle_activity.contraction.GFP_mean_WO(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.contraction.GFP(p, c, i, index, t)));
                
                % N45
                AGSICI_muscle_activity.N45.GFP_mean(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.N45.GFP(p, c, i, :, t)));
                AGSICI_muscle_activity.N45.GFP_mean_WO(p, c, i, t) = mean(squeeze(AGSICI_muscle_activity.N45.GFP(p, c, i, index, t)));
            end
        end
    end
end
clear p c i t

% append new variables to the general MATLAB file
save(results_file, 'AGSICI_muscle_activity', '-append');

%% 8) muscle contraction X N45: export for letswave for browsing
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            % put individual data in letswave format
            for s = 1:length(subject) 
                % muscle activity
                data_m(s, 1, 1, 1, 1, :) = AGSICI_muscle_activity.contraction.GFP(p, c, i, s, :);
                
                % N45
                data_c(s, 1, 1, 1, 1, :) = AGSICI_muscle_activity.N45.GFP(p, c, i, s, :);
            end
            
            % modify header
            header.datasize = size(data_c);
            header.xstart = time_window(1);
            header.chanlocs(2:end) = [];
            header.chanlocs.labels = 'GFP';
            header.chanlocs.topo_enabled = 0;
            
            % save individual datasets
            data = data_m; 
            header.name = ['merged muscle contraction ' position{p} ' ' current{c} ' ' intensity{i}];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header')
            data = data_c; 
            header.name = ['merged muscle N45 ' position{p} ' ' current{c} ' ' intensity{i}];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header')
            
            % save averaged data
            clear data
            data(1, 1, 1, 1, 1, :) = AGSICI_muscle_activity.contraction.GFP_mean(p, c, i, :);
            header.size = size(data);
            header.name = ['avg merged muscle contraction ' position{p} ' ' current{c} ' ' intensity{i}];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header')            
            data(1, 1, 1, 1, 1, :) = AGSICI_muscle_activity.N45.GFP_mean(p, c, i, :);
            header.name = ['avg merged muscle N45 ' position{p} ' ' current{c} ' ' intensity{i}];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header') 
        end
    end
end
clear p c i s data_m data_c data 

%% 9) muscle contraction X N45: plot mean GFP - without outliers
% ----- adjustable parameters -----
cont_window = [-0.01, 0.05]; 
% ----- adjustable parameters -----

% recalculate x and limits for muscular contraction
x_cont = [cont_window(1):xstep:cont_window(2)];
x_start_cont = (cont_window(1) - time_window(1))/xstep;
x_end_cont = (cont_window(2) - time_window(1))/xstep;

% plot mean GFP per condition
for p = 1:length(position)
    for c = 1:length(current)
        % launch the figure
        fig = figure(figure_counter);        
        sgtitle(sprintf('%s STS - %s current', position{p}, current{c}), 'FontSize', 16, 'FontWeight', 'bold')
        
        % plot GFP of muscular activity, mark TOI
        subplot(2, 1, 1)
        hold on
        rectangle('Position', [TOI_m(1), 2, TOI_m(2) - TOI_m(1), 200-2], 'FaceColor', [0.98, 0.83, 0.83], 'EdgeColor', 'none')
        plot_TEP(x_cont, squeeze(AGSICI_muscle_activity.contraction.GFP_mean_WO(p, c, :, x_start_cont:x_end_cont)), [-0.002 0.002], ...
            'limit', [0 200], 'legend', {'100 %rMT' '120 %rMT' '140 %rMT'})
        title('GFP of raw data, left hemisphere')
        
        % plot GFP of N45 component, mark TOI
        subplot(2, 1, 2)
        hold on
        rectangle('Position', [TOI_N45(1), 0.025, TOI_N45(2) - TOI_N45(1), 3.5-0.025], 'FaceColor', [0.98, 0.83, 0.83], 'EdgeColor', 'none')
        plot_TEP(x, squeeze(AGSICI_muscle_activity.N45.GFP_mean_WO(p, c, :, :)), [-0.005 0.01], 'limit', [0 3.5])
        title('GFP of processed data')    
        
        % save figure
        figure_name = sprintf('AGSICI_GFP_%s_%s_WO', position{p}, current{c});
        savefig([output_folder '\' figure_name '.fig'])
        saveas(fig, [output_folder '\' figure_name '.png'])
        
        % update counter
        figure_counter = figure_counter + 1;
    end
end
clear p c fig cont_window x_cont x_start_cont x_end_cont

%% 10) muscle contraction X N45: plot mean values - without outliers
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
            y(i) = mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index));
            SEM(i) = std(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index)) / sqrt(length(subject(index)));
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
title(sprintf('Muscular activity (GFP - AUC): %d - %dms', TOI_m(1) * 1000, TOI_m(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('AUC (\muV*s \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);
yl = get(gca, 'ylim');
text(2.1, yl(2) - 250, sprintf('without outliers:\nsubj. %d, %d, %d', outliers(1), outliers(2), outliers(3)), 'FontSize', 14)

% save the figure       
savefig([output_folder '\AGSICI_cont_all_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_all_WO.png'])

% update the counter
figure_counter = figure_counter + 1;  
clear fig leg perr y SEM yl

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
            y(i) = mean(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index));
            SEM(i) = std(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index)) / sqrt(length(subject(index)));
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
title(sprintf('Component N45 (GFP - AUC): %d - %dms', TOI_N45(1) * 1000, TOI_N45(2) * 1000), ...
    'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (%rMT)'); ylabel('AUC (\muV*s \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
set(leg, 'Location','northwest', 'FontSize', 14);
yl = get(gca, 'ylim');
text(2.1, yl(2) - 8, sprintf('without outliers:\nsubj. %d, %d, %d', outliers(1), outliers(2), outliers(3)), 'FontSize', 14)

% save the figure       
savefig([output_folder '\AGSICI_N45_all.fig'])
saveas(fig, [output_folder '\AGSICI_N45_all.png'])

% update the counter
figure_counter = figure_counter + 1;  
clear fig leg perr y SEM yl

%% 11) muscle contraction: distribution across conditions - without outliers
% calculate mean values to plot
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_visual(i, (p-1)*2 + c) = mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index));
        end
    end
end
clear p c i 

% plot average GFP-AUC by condition
fig = figure(figure_counter)
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along - normal', 'along - reversed', 'across - normal', 'across - reversed'}, ...
    'Location', 'bestoutside', 'fontsize', 14)
xlabel('intensity of stimulation (%rMT)')
ylabel('muscular contraction')
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(0.75, yl(2) - 250, sprintf('without outliers:\nsubj. %d, %d, %d', outliers(1), outliers(2), outliers(3)), 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_size_all_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size_all_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig yl barplot a data_visual 

% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        data_visual(:, (p-1)*2 + c) = squeeze(mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, :, index), 3));
    end
end
clear p c

% plot average RMS - boxplot
fig = figure(figure_counter); ax = gca;       
hold on
boxplot(data_visual, 'color', colours)
ax.XTickLabel = '';
label_array = {'along' 'along' 'across' 'across'; 'normal' 'reversed' 'normal' 'reversed'}; 
for i = 1:length(label_array)
    text(i, ax.YLim(1), sprintf('%s\n%s', label_array{:, i}), 'FontSize', 14, ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');    
end
ylabel('muscular contraction - GFP-AUC')
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(1, yl(2) - 500, sprintf('without outliers:\nsubj. %d, %d, %d', outliers(1), outliers(2), outliers(3)), 'FontSize', 14)
clear i label_array yl

% plot the markers
for b = 1:4
    scat(b) = scatter(repelem(b, size(data_visual, 1)), data_visual(:, b),...
        75, colours(b, :), 'filled');
end

% save the figure
savefig([output_folder '\AGSICI_cont_size_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig ax scat b data_visual 

%% 12) muscle contraction X N45: overall correlation - without outliers
% extract data
data_corr_m = [];
data_corr_N45 = [];
for p = 1:length(position)
    for c = 1:length(current)        
        for i = 1:length(intensity)
            data_corr_m = [data_corr_m; squeeze(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index))];
        end
    end
end
data_corr = [data_corr_m data_corr_N45];
clear p c i data_corr_m data_corr_N45

% choose colours
marker_col = [];
for a = 1:length(position) * length(current)
    for b = 1:length(intensity) * length(subject(index))
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
savefig([output_folder '\AGSICI_cont_corr_all_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_corr_all_WO.png'])

% update the counter
figure_counter = figure_counter + 1;  

% ----- ranked correlation -----
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
savefig([output_folder '\AGSICI_cont_corr_all_WO_ranked.fig'])
saveas(fig, [output_folder '\AGSICI_cont_corr_all_WO_ranked.png'])

% update the counter
figure_counter = figure_counter + 1;  

clear data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 13) muscle contraction X N45: correlation per orientation - without outliers
for p = 1:length(position)
    for c = 1:length(current)       
        % extract data
        data_corr_m = [];
        data_corr_N45 = []; 
        for i = 1:length(intensity)
            data_corr_m = [data_corr_m; squeeze(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index))];
        end
        data_corr = [data_corr_m data_corr_N45];
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
        savefig([output_folder '\AGSICI_cont_corr_' position{p} '-' current{c} '_WO.fig'])
        saveas(fig, [output_folder '\AGSICI_cont_corr_' position{p} '-' current{c} '_WO.png'])

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
        savefig([output_folder '\AGSICI_cont_corr_ranked_' position{p} '-' current{c} '_WO.fig'])
        saveas(fig, [output_folder '\AGSICI_cont_corr_ranked_' position{p} '-' current{c} '_WO.png'])

        % update the counter
        figure_counter = figure_counter + 1;  
    end
end
clear p c data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 14) muscle contraction X N45: correlation per intensity - without outliers 
% choose colours
marker_col = [];
for a = 1:length(position) * length(current)
    for b = 1:length(subject(index))
        marker_col = [marker_col; colours(a, :)];
    end
end
clear a b

% plot
for i = 1:length(intensity)         
    % extract data
    data_corr_m = [];
    data_corr_N45 = []; 
    for p = 1:length(position)
        for c = 1:length(current)   
            data_corr_m = [data_corr_m; squeeze(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index))];
        end
    end
    data_corr = [data_corr_m data_corr_N45];
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
    savefig([output_folder '\AGSICI_cont_corr_int_' intensity{i}(end-2:end) '_WO.fig'])
    saveas(fig, [output_folder '\AGSICI_cont_corr_int_' intensity{i}(end-2:end) '_WO.png'])

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
    savefig([output_folder '\AGSICI_cont_corr_int_' intensity{i}(end-2:end) '_WO_ranked.fig'])
    saveas(fig, [output_folder '\AGSICI_cont_corr_int_' intensity{i}(end-2:end) '_WO_ranked.png'])

    % update the counter
    figure_counter = figure_counter + 1;  
end
clear data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% ) tonic muscular activity - extract/export data
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

% export as a long-format table
AGSICI_tonic = table;
AGSICI_tonic.subject = zeros(0); AGSICI_tonic.orientation = {}; AGSICI_tonic.intensity = zeros(0); 
AGSICI_tonic.baseline_RMS = zeros(0);  
row_cnt = 1;
for s = 1:length(subject)
    for p = 1:length(position)
        for c = 1:length(current)
            for i = 1:length(intensity)
                % fill in the row
                AGSICI_tonic.subject(row_cnt) = subject(s);             
                AGSICI_tonic.orientation(row_cnt) = {[position{p} '-' current{c}]};
                AGSICI_tonic.intensity(row_cnt) = str2double(intensity{i}(end-2:end));
                AGSICI_tonic.baseline_RMS(row_cnt) = AGSICI_muscle_activity.tonic(p, c, i, s);
                
                % update row counter
                row_cnt = row_cnt + 1;
            end
        end
    end
end
writetable(AGSICI_tonic, [output_folder '\AGSICI_tonic.csv'])
clear p c s i row_cnt

% append new variables to the general MATLAB file
save(results_file, 'AGSICI_muscle_activity', 'AGSICI_tonic', '-append');

%% ) tonic muscular activity - plot
% calculate mean values to plot
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_visual(i, (p-1)*2 + c) = mean(AGSICI_muscle_activity.tonic(p, c, i, :));
        end
    end
end
clear p c i 

% plot average RMS by condition
fig = figure(figure_counter)
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along - normal', 'along - reversed', 'across - normal', 'across - reversed'}, ...
    'Location', 'bestoutside', 'fontsize', 14)
xlabel('intensity of stimulation (%rMT)')
ylabel('baseline RMS')
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_tonic_bl_all.fig'])
saveas(fig, [output_folder '\AGSICI_tonic_bl_all.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig barplot a data_visual 

% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        data_visual(:, (p-1)*2 + c) = squeeze(mean(AGSICI_muscle_activity.tonic(p, c, :, :), 3));
    end
end
clear p c

% plot average RMS - boxplot
fig = figure(figure_counter); ax = gca;       
hold on
boxplot(data_visual, 'color', colours)
ax.XTickLabel = '';
label_array = {'along' 'along' 'across' 'across'; 'normal' 'reversed' 'normal' 'reversed'}; 
for i = 1:length(label_array)
    text(i, ax.YLim(1), sprintf('%s\n%s', label_array{:, i}), 'FontSize', 14, ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');    
end
ylabel('baseline RMS')
set(gca, 'Fontsize', 14)
clear i label_array 

% plot the markers
for b = 1:4
    scat(b) = scatter(repelem(b, size(data_visual, 1)), data_visual(:, b),...
        75, colours(b, :), 'filled');
end

% save the figure
savefig([output_folder '\AGSICI_tonic_bl.fig'])
saveas(fig, [output_folder '\AGSICI_tonic_bl.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig ax scat b data_visual 

%% 13) muscle contraction: stratification, all subjects included
% ----- adjustable parameters -----
prefix_m_cont = 'avg avgchan bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
% ----- adjustable parameters -----
% split data to groups by size of contraction
data_table = AGSICIcontXN45;
data_1 = []; data_2 = []; data_3 = []; data_4 = []; 
case_counter = ones(4, 2, 3);
for a = 1:height(data_table)
    % identify dataset
    orientation = char(data_table.orientation(a));
    if strcmp(orientation(1:3), 'alo')
        p = 1; 
    elseif strcmp(orientation(1:3), 'acr')
        p = 2;
    end
    if strcmp(orientation(end-2:end), 'mal')
        c = 1;
    elseif strcmp(orientation(end-2:end), 'sed')
        c = 2;
    end
    for b = 1:length(intensity)
        if str2double(intensity{b}(end-2:end)) == data_table.intensity(a)
            i = b;
        end
    end    
    if data_table.subject(a) < 10
        subj = ['0' num2str(data_table.subject(a))];
    else
        subj = num2str(data_table.subject(a));
    end
    
    % fill in counter
    case_counter(data_table.contractioncategory(a), p, i) =  case_counter(data_table.contractioncategory(a), p, i) + 1;   
    
%     % load the processed dataset
%     dataset_name = [prefix_m_cont ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
%     load(dataset_name)
%     
%     % append data to the group variable  
%     statement = ['data_' num2str(data_table.contractioncategory(a)) '(end+1, :, 1, 1, 1, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));'];
%     eval(statement)
end
clear a b c p i orientation 

% modify header
header.datasize = size(data_1);
header.xstart = time_window(1);
header.chanlocs = header.chanlocs(1:32);

% save for letswave
for d = 1:4
    statement = ['data = data_' num2str(d) ';'];
    eval(statement)
    header.name = ['merged muscle contraction group_' num2str(d)];
    save([header.name '.mat'], 'data')
    save([header.name '.lw6'], 'header')
end
clear d statement data_1 data_2 data_3 data_4

% ----- counts per position -----
data_visual = squeeze(sum(case_counter(:, :, :), 3));
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'stacked', 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along' 'across'}, 'Location', 'bestoutside', 'fontsize', 14)
xlabel('muscular contraction')
ylabel('count')
set(gca, 'xtick', 1:4, 'xticklabel', {'smallest' '' '' 'largest'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(0.75, yl(2) - 250, 'all subjects', 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_counts_position.fig'])
saveas(fig, [output_folder '\AGSICI_cont_counts_position.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear a fig barplot data_visual 

% ----- counts per intensity -----
col = [0.24 0.49 0.99; 0.72 0.27 1; 0.87 0.16 0.40];
data_visual = squeeze(sum(case_counter(:, :, :), 2));
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'stacked', 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = col(a, :);
end
legend(barplot, {'100 %rMT' '120 %rMT' '140 %rMT'}, 'Location', 'bestoutside', 'fontsize', 14)
xlabel('muscular contraction')
ylabel('count')
set(gca, 'xtick', 1:4, 'xticklabel', {'smallest' '' '' 'largest'})
set(gca, 'Fontsize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_counts_intensity.fig'])
saveas(fig, [output_folder '\AGSICI_cont_counts_intensity.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear col a fig barplot data_visual 

%% 14) muscle contraction: stratification per intensity, all subjects included

%% 15) muscle contraction: distribution across conditions, without outliers
% ----- adjustable parameters -----
outliers = [5 12];
% ----- adjustable parameters -----
% calculate subject index
index = subject ~= outliers(1) & subject ~= outliers(2);

% calculate mean values to plot
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            data_visual(i, (p-1)*2 + c) = mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index));
        end
    end
end
clear p c i 

% plot average GFP-AUC by condition
fig = figure(figure_counter)
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along - normal', 'along - reversed', 'across - normal', 'across - reversed'}, ...
    'Location', 'bestoutside', 'fontsize', 14)
xlabel('intensity of stimulation (%rMT)')
ylabel('muscular contraction')
set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(0.75, yl(2) - 250, sprintf('without outliers:\nsubj. %d and %d', outliers(1), outliers(2)), 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_size_all_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size_all_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear yl fig barplot a data_visual 

% average across intensities
for p = 1:length(position)
    for c = 1:length(current)
        data_visual(:, (p-1)*2 + c) = squeeze(mean(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, :, index), 3));
    end
end
clear p c

% plot average RMS - boxplot
fig = figure(figure_counter); ax = gca;       
hold on
boxplot(data_visual, 'color', colours)
ax.XTickLabel = '';
label_array = {'along' 'along' 'across' 'across'; 'normal' 'reversed' 'normal' 'reversed'}; 
for i = 1:length(label_array)
    text(i, ax.YLim(1), sprintf('%s\n%s', label_array{:, i}), 'FontSize', 14, ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');    
end
ylabel('muscular contraction')
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
text(2, yl(2) - 1000, sprintf('without outliers:\nsubj. %d and %d', outliers(1), outliers(2)), 'FontSize', 14)
clear i label_array 

% plot the markers
for b = 1:4
    scat(b) = scatter(repelem(b, size(data_visual, 1)), data_visual(:, b),...
        75, colours(b, :), 'filled');
end

% save the figure
savefig([output_folder '\AGSICI_cont_size_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_size_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear yl fig ax scat b data_visual 

%% 16) muscle contraction: stratification, without outliers
% ----- adjustable parameters -----
prefix_m_cont = 'avg avgchan bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
% ----- adjustable parameters -----

% split data to groups by size of contraction
data_table = AGSICIcontXN45;
data_1 = []; data_2 = []; data_3 = []; data_4 = []; 
case_counter = ones(4, 2, 3);
for a = 1:height(data_table)
    % identify dataset
    orientation = char(data_table.orientation(a));
    if strcmp(orientation(1:3), 'alo')
        p = 1; 
    elseif strcmp(orientation(1:3), 'acr')
        p = 2;
    end
    if strcmp(orientation(end-2:end), 'mal')
        c = 1;
    elseif strcmp(orientation(end-2:end), 'sed')
        c = 2;
    end
    for b = 1:length(intensity)
        if str2double(intensity{b}(end-2:end)) == data_table.intensity(a)
            i = b;
        end
    end    
    if data_table.subject(a) < 10
        subj = ['0' num2str(data_table.subject(a))];
    else
        subj = num2str(data_table.subject(a));
    end
    
    % fill in counter
    case_counter(data_table.contractioncategory(a), p, i) =  case_counter(data_table.contractioncategory(a), p, i) + 1;   
    
    % load the processed dataset
    dataset_name = [prefix_m_cont ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
    load(dataset_name)
    
    % append data to the group variable  
    statement = ['data_' num2str(data_table.contractioncategory(a)) '(end+1, :, 1, 1, 1, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));'];
    eval(statement)
end
clear a b c p i orientation 

% modify header
header.datasize = size(data_1);
header.xstart = time_window(1);
header.chanlocs = header.chanlocs(1:32);

% save for letswave
for d = 1:4
    statement = ['data = data_' num2str(d) ';'];
    eval(statement)
    header.name = ['merged muscle contraction group_' num2str(d) ' WO'];
    save([header.name '.mat'], 'data')
    save([header.name '.lw6'], 'header')
end
clear d statement data_1 data_2 data_3 data_4

% ----- counts per position -----
data_visual = squeeze(sum(case_counter(:, :, :), 3));
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'stacked', 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = colours(a, :);
end
legend(barplot, {'along' 'across'}, 'Location', 'bestoutside', 'fontsize', 14)
xlabel('muscular contraction')
ylabel('count')
set(gca, 'xtick', 1:4, 'xticklabel', {'smallest' '' '' 'largest'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
ylim([yl(1) yl(2) + 10])
text(1.5, yl(2)+5, 'without outliers', 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_counts_position_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_counts_position_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear a fig barplot yl ata_visual 

% ----- counts per intensity -----
col = [0.24 0.49 0.99; 0.72 0.27 1; 0.87 0.16 0.40];
data_visual = squeeze(sum(case_counter(:, :, :), 2));
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'stacked', 'EdgeColor', 'none');
for a = 1:size(data_visual, 2)
    barplot(a).FaceColor = col(a, :);
end
legend(barplot, {'100 %rMT' '120 %rMT' '140 %rMT'}, 'Location', 'bestoutside', 'fontsize', 14)
xlabel('muscular contraction')
ylabel('count')
set(gca, 'xtick', 1:4, 'xticklabel', {'smallest' '' '' 'largest'})
set(gca, 'Fontsize', 14)
yl = get(gca, 'ylim');
ylim([yl(1) yl(2) + 10])
text(1.5, yl(2)+5, 'without outliers', 'FontSize', 14)
hold off

% save the figure
savefig([output_folder '\AGSICI_cont_counts_intensity_WO.fig'])
saveas(fig, [output_folder '\AGSICI_cont_counts_intensity_WO.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear col a fig yl barplot data_visual 

%% 17) muscle contraction: stratification per factors, without outliers
% ----- adjustable parameters -----
prefix_m_cont = 'avg avgchan bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
% ----- adjustable parameters -----

% split by intensity
for i = 1:length(intensity)
    data_1 = []; data_2 = []; 
    
    % subset the table
    rows = data_table.intensity == str2double(intensity{i}(end-2:end));
    data_table_int = data_table(rows, :);
    
    % re-assign to groups
    data_table_int.contractioncategory(1:height(data_table_int)/2) = 1;
    data_table_int.contractioncategory(height(data_table_int)/2+1:end) = 2;
    
    % extract the data
    for a = 1:height(data_table_int)
        % identify dataset
        orientation = char(data_table_int.orientation(a));
        if strcmp(orientation(1:3), 'alo')
            p = 1; 
        elseif strcmp(orientation(1:3), 'acr')
            p = 2;
        end
        if strcmp(orientation(end-2:end), 'mal')
            c = 1;
        elseif strcmp(orientation(end-2:end), 'sed')
            c = 2;
        end 
        if data_table_int.subject(a) < 10
            subj = ['0' num2str(data_table_int.subject(a))];
        else
            subj = num2str(data_table_int.subject(a));
        end 

        % load the processed dataset
        dataset_name = [prefix_m_cont ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
        load(dataset_name)
        
        % append data to the group variable 
        statement = ['data_' num2str(data_table_int.contractioncategory(a)) '(end+1, :, 1, 1, 1, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));'];
        eval(statement)

    end
    
    % save for letswave
    for d = 1:2
        statement = ['data = data_' num2str(d) ';'];
        eval(statement)
        header.name = ['merged muscle contraction ' intensity{i} ' group_' num2str(d) ' WO'];
        save([header.name '.mat'], 'data')
        save([header.name '.lw6'], 'header')
    end
    clear d statement data_1 data_2 
end
clear a c p i orientation 

% split by orientation
orientation = {'along-normal' 'along-reversed' 'across-normal' 'across-reversed'};
for p = 1:length(position)
    for c = 1:length(current)
        data_1 = []; data_2 = []; 

        % subset the table
        rows = data_table.orientation == orientation{(p-1)*2+c};
        data_table_int = data_table(rows, :);

        % re-assign to groups
        data_table_int.contractioncategory(1:round(height(data_table_int)/2)) = 1;
        data_table_int.contractioncategory(round(height(data_table_int)/2)+1:end) = 2;

        % extract the data
        for a = 1:height(data_table_int)
            % identify intensity
            for b = 1:length(intensity)
                if str2double(intensity{b}(end-2:end)) == data_table_int.intensity(a)
                    i = b;
                end
            end 
            
            % identify subject
            if data_table_int.subject(a) < 10
                subj = ['0' num2str(data_table_int.subject(a))];
            else
                subj = num2str(data_table_int.subject(a));
            end 

            % load the processed dataset
            dataset_name = [prefix_m_cont ' ' subj ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'];
            load(dataset_name)

            % append data to the group variable 
            statement = ['data_' num2str(data_table_int.contractioncategory(a)) '(end+1, :, 1, 1, 1, :) = squeeze(data(:, 1:32, :, :, :, x_start:x_end));'];
            eval(statement)

        end

        % save for letswave
        for d = 1:2
            statement = ['data = data_' num2str(d) ';'];
            eval(statement)
            header.name = ['merged muscle contraction ' position{p} ' ' current{c} ' group_' num2str(d) ' WO'];
            save([header.name '.mat'], 'data')
            save([header.name '.lw6'], 'header')
        end
        clear d statement data_1 data_2 
    end
end
clear a b c p i orientation 

%% functions
function plot_TEP(x, data_visual, interpol_lim, varargin)
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
rectangle('Position', [interpol_lim(1), lim(1), interpol_lim(2) - interpol_lim(1), lim(2) - lim(1)], ...
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
xlabel('muscular activity (GFP - AUC)'); ylabel('N45 (GFP - AUC)');
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