%% AG-SICI: TMS-EVOKED POTENTIALS
% Written by Dominika for AG-SICI project (2021)
% 
% Colection of scripts to analyse preprocessed TMS-EEG data:
% 1) load the data
%       - loads individual averaged data and trims it in a predefined time
%       window
%       - saves the data from all subjects and conditions in a 6D matrix:
%       2 (position) X 2 (current) X 3 (intensity) X number of subjects X
%       number of channels X number of timepoints
%       - the sequence of all processed subjects is saved as 'subject_order'
%           --> in case there are some data added to already existing
%           dataset, subject = newly added datasets 
% 
% 2) preliminary visualization of TEPs 
%       - averages data from all participants for all conditions, calculates CI 
%       - plots separately 4 coil placements to show differences in
%       intensities
%       - averages across intensities and plots mean + CI of all placement
%       conditions in one figure
% 
% 3) calculate global field power from grand average
%       - calculates GFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in the outcome structure
% 
% 4) define peak widths using grand average GFP - just informative 
%       - for each condition separately, plots and saves the figure 
%       - uses function findpeaks awith the 'halfheight' option
%       - saves peak widths in the outcome structure
% 
% 5) peak tracking: topographies + electrodes of interest
%       - lets user decide default latency and span for each TEP peak
%       --> starting peak point to track peak progress
%       - for each peak/condition plots the window of interest (defined by
%       default peak ceter and span) and lets user adjust manually the
%       exact peak latency --> predefined seed electrodes for visualization
%       - puts 0 timepoint to the peak latency --> 'centered' dataset
%       - subtracts the value at t0 from the rest of the data to track the
%       progress of the peak --> 'subtracted' dataset
%       - averages the datasets across subjects and saves for letswave
%       - extracts individual peak latency and averages across subjects
% 
% 6) GFP amplitude extraction
%
% 7) TEP amplitude extraction
 

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
subject = [1, 3:18, 20, 21];
subject_order = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};
prefix = 'avg avgchan bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
filename = 'AG-SICI_plus';

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

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
clear a answer

% load a random header
load([prefix ' 01 across normal stim_100.lw6'], '-mat')
labels = {header.chanlocs.labels};
header.chanlocs(33) = [];

% visualization 
figure_counter = 1;
xstep = header.xstep;
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

% create folder for figures
output_file = [pwd '\' filename '.mat'];
folder_figures = [pwd '\' filename '_figs'];
if ~exist(folder_figures) 
    mkdir(folder_figures);
end

clear c p 
%% 1) load the data
% load data based on the prefix + conditions
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject)
                % load individual data from letswave
                if subject(s) < 10
                    load([prefix ' 0' num2str(subject(s)) ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'])
                else
                    load([prefix ' ' num2str(subject(s)) ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'])
                end
                
                % append the data in the data matrix
                A(p, c, i, s, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));
            end
        end
    end 
end
disp(['Datasize: ' num2str(size(A))])
clear p c i s data 

% save dataset to the global MATLAB file
if exist([filename '.mat']) == 0
    AGSICI_data = A;
    save([filename '.mat'], 'AGSICI_data');
else
    load([filename '.mat'])
    if exist('AGSICI_data') == 0
        AGSICI_data = A;
        save([filename '.mat'], 'AGSICI_data', '-append');
    else
        AGSICI_data = cat(4, AGSICI_data, A);    
        save([filename '.mat'], 'AGSICI_data', '-append');
    end
end
clear A

%% 2) preliminary TEP visualization 
% ----- decide output parameters -----
electrode = {'target'};
y_limits = [-4, 5.5];
line_type = {':' '-.' '-'};
% ------------------------------------
% read saved data
load(filename)

% average data, calculate CI
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            % average across subjects
            for e = 1:size(AGSICI_data, 5)
                for t = 1:size(AGSICI_data, 6)
                    AGSICI_data_mean(p, c, i, e, t) = mean(squeeze(AGSICI_data(p, c, i, :, e, t)));
                    AGSICI_data_CI(p, c, i, e, t) = (std(squeeze(AGSICI_data(p, c, i, :, e, t)))/sqrt(length(subject))) * z;
                end
            end 
        end        
        % average across intensities
        for e = 1:size(AGSICI_data, 5)
            for t = 1:size(AGSICI_data, 6)
                gfp_data(p, c, e, t) = mean(squeeze(AGSICI_data_mean(p, c, :, e, t)));
                gfp_CI(p, c, e, t) = mean(squeeze(AGSICI_data_CI(p, c, :, e, t)));
            end
        end 
    end
end
clear p c i e t 

% plot TEPs of all conditions
for e = 1:length(electrode)     
    % launch the figure
    fig = figure(figure_counter);
    hold on
        
    % loop through datasets
    for p = 1:length(position)
        for c = 1:length(current)
            % plot each tested position in a separate subplot
            subplot(2, 2, (p - 1)*2 + c)
            
            for i = 1:length(intensity)  
                % prepare data
                data_visual = squeeze(AGSICI_data_mean(p, c, i, find(contains(labels, electrode{e})), :));  
                
                if i == 1
                    % shade interpolated interval 
                    plot(x, data_visual, 'b:', 'LineWidth', 0.5);
                    hold on
                    rectangle('Position', [0, y_limits(1), 0.01, y_limits(2) - y_limits(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
                end

                % plot data     
                P(i) = plot(x, data_visual, 'Color', colours((p - 1)*2 + c, :), 'LineWidth', 2, 'LineStyle', line_type{i});
                hold on
            end 
            
            % add other parameters
            title([position{p} ' STS, ' current{c} ' current'])
            xlabel('time (s)')
            ylabel('amplitude (\muV)')
            set(gca, 'FontSize', 14)
            xlim(time_window)
            ylim(y_limits)
            line([0, 0], y_limits, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
%             lgd = legend(P, {'100 %rMT' '120 %rMT' '140 %rMT'}, 'Location', 'southwest');
%             lgd.FontSize = 14; title(lgd, 'Stimulation intensity')
        end
    end    
    hold off
    
    % name and save figure
    figure_name = ['AGSICI_TEP_' electrode{e}];
    savefig([pwd '\' filename '_figs\' figure_name '.fig'])
    saveas(fig, [pwd '\' filename '_figs\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear e p c i fig data_visual P lgd figure_name

% plot averaged TEPs in one figure
for e = 1:length(electrode)     
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % shade interpolated interval 
    plot(x, squeeze(gfp_data(1, 1, find(contains(labels, electrode{e})), :)), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, y_limits(1), 0.01, y_limits(2) - y_limits(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
        
    % loop through datasets
    for p = 1:length(position)
        for c = 1:length(current) 
            % prepare data
            data_visual = squeeze(gfp_data(p, c, find(contains(labels, electrode{e})), :))'; 
            CI_visual = squeeze(gfp_CI(p, c, find(contains(labels, electrode{e})), :))'; 

            % plot data     
            P((p - 1)*2 + c) = plot(x, data_visual, 'Color', colours((p - 1)*2 + c, :), 'LineWidth', 2.5);
            F((p - 1)*2 + c) = fill([x fliplr(x)],[data_visual + CI_visual fliplr(data_visual - CI_visual)], ...
                colours((p - 1)*2 + c, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        end
    end
    
    % add other parameters
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(y_limits)
    line([0, 0], y_limits, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
%     lgd = legend(P, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'}, 'Location', 'south', 'numcolumns', 2);
%     lgd.FontSize = 14; title(lgd, 'Tested conditions')
    
    % name and save figure
    figure_name = ['AGSICI_TEP_all_' electrode{e}];
    savefig([pwd '\' filename '_figs\' figure_name '.fig'])
    saveas(fig, [pwd '\' filename '_figs\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear e p c fig data_visual P F lgd figure_name CI_visual

% save dataset to the global MATLAB file
save([filename '.mat'], 'AGSICI_data_mean', 'AGSICI_data_CI', 'gfp_data', 'gfp_CI', '-append');
clear electrode AGSICI_data_CI 

%% 3) GFP
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------

% compute individual GFP
AGSICI_GFP_subject = struct;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for s = 1:length(subject_order)
                % calculate GFP (exclude target channel)
                AGSICI_GFP_subject.data(p, c, i, s, :) = std(squeeze(AGSICI_data(p, c, i, s, 1:32, :)), 1);  
            end
        end
    end
end
clear p c i s

% compute grand average GFP, plot
row_count = 1;
for p = 1:length(position)
    for c = 1:length(current)
        % calculate GFP (exclude target channel)
        AGSICI_GFP(p, c, :) = std(squeeze(gfp_data(p, c, 1:size(gfp_data, 3) - 1, :)), 1);  

        % set dataset name + figure title
        fig_name = ['AG-SICI_GFP_' position{p} '_' current{c}];
        fig_title = ['GFP: ' position{p} ' STS, ' current{c} ' current'];

        % plot GFP and extract peak latencies
        fig = figure(figure_counter);
        hold on

        if ~isempty(max_peaks)
        % plot GFP
        h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
        AGSICI_TEP(row_count).latencies = gfp_plot(x, squeeze(AGSICI_GFP(p, c, :)), time_window, xstep, labeled, 'max_peaks', max_peaks);
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')

        % add topoplots
        for t = 1:length(AGSICI_TEP(row_count).latencies)
            % choose data for topoplot 
            for e = 1:size(AGSICI_data, 5)
                data_topoplot(1, e, 1, 1, 1, :) = squeeze(gfp_data(p, c, e, :));
            end
            
            % plot the topoplot
            h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);            
            topo_plot(header, data_topoplot, AGSICI_TEP(row_count).latencies(t), time_window(1), [-2, 2])

            % shift down
            pos = get(h_axis(1 + t), 'Position');
            pos(2) = pos(2) - 0.05;
            set(h_axis(1 + t), 'Position', pos);

            % add timing
            text(-0.3, -0.8, sprintf('%1.0f ms', AGSICI_TEP(row_count).latencies(t)*1000), 'Color', [1 0 0], 'FontSize', 14)
        end

        else
            AGSICI_TEP(row_count).latencies(:) = gfp_plot(x, squeeze(AGSICI_GFP(p, c, :)), time_window, xstep, labeled);
        end  
        hold off

        % save figure, update    
        savefig([pwd '\' filename '_figs\' fig_name '.fig'])
        saveas(fig, [pwd '\' filename '_figs\' fig_name '.png'])
        figure_counter = figure_counter + 1;
        
        % update row counter
        row_count = row_count + 1;
    end
end
clear p c t h_axis pos fig_name fig_title row_count fig data_topoplot

% append new variables to the general MATLAB file
save([filename '.mat'], 'AGSICI_GFP', 'AGSICI_TEP', 'AGSICI_GFP_subject', '-append')
clear labeled max_peaks gfp_CI

%% 4) peak widths
% ----- decide output parameters -----
x_limit = [0.015, 0.220];
% ------------------------------------

% calculate narrow window parameters
x_start_narrow = (x_limit(1) - time_window(1))/xstep;
x_end_narrow = (x_limit(2) - time_window(1))/xstep;
x_narrow = [x_limit(1) : xstep: x_limit(2)];

% loop through conditions
row_count = 1;
for p = 1:length(position)
    for c = 1:length(current)
        % choose data and x
        data = squeeze(AGSICI_GFP(p, c, x_start_narrow:x_end_narrow));

        % identify peak widths
        [P, L, AGSICI_TEP(row_count).widths, R] = findpeaks(data, 'Annotate','extents', ...
            'WidthReference', 'halfheight','MinPeakDistance', 10, 'MinPeakProminence', 0.015);
        AGSICI_TEP(row_count).widths = ceil(AGSICI_TEP(row_count).widths) * header.xstep;
        if length(AGSICI_TEP(row_count).widths) ~= length(AGSICI_TEP(row_count).latencies)
            disp('ATTENTION: number of identified peaks does not match with peaks extracted in the previous step!')
        end

        % plot figure
        fig = figure(figure_counter);
        hold on
        findpeaks(data, x_narrow, 'Annotate', 'extents', 'WidthReference', 'halfheight', 'MinPeakProminence', 0.015);
        set(gca, 'fontsize', 14)
        xlabel('time(s)')
        ylabel('power (\muV^2)')
        title([position{p} ' STS, ' current{c} ' current'])
        grid off
        legend off

        % add width denotation
        for k = 1:numel(AGSICI_TEP(row_count).latencies)
            if k == numel(AGSICI_TEP(row_count).latencies)                
                text(AGSICI_TEP(row_count).latencies(k) - 0.005, -0.25, ...
                    sprintf('%1.0f ms', AGSICI_TEP(row_count).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
            else
                text(AGSICI_TEP(row_count).latencies(k) - 0.005, -0.25, ...
                    sprintf('%1.0f', AGSICI_TEP(row_count).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
            end
        end
        ylim([-0.5, 2])
        hold off
        
        % save figure, update
        savefig([pwd '\' filename '_figs\AG-SICI_widths_' position{p} '_' current{c}])
        saveas(fig, [pwd '\' filename '_figs\AG-SICI_widths' position{p} '_' current{c} '.png'])
        figure_counter = figure_counter + 1;
        
        % update row counter
        row_count = row_count + 1;
    end
end

% append new variables to the general MATLAB file
save([filename '.mat'], 'AGSICI_TEP', '-append');
clear x_limit x_start_narrow x_end_narrow x_narrow c p k data fig P L R row_count

%% 5) peak tracking
% ----- decide output parameters -----
seed_electrode = {'target' 'Cz'};                                   % electrode that will be used to set 0 timepoint
seed_peaks = {1:2; 3:6};                                            % which peeks use which seed electrode
AGSICI_TEP_avg.peak = {'P25' 'N40' 'N45' 'P75' 'N100' 'P180'};      % choose peak names
AGSICI_TEP_avg.center = [0.024, 0.038, 0.047, 0.075, 0.118, 0.200]; % choose default peak centers
AGSICI_TEP_avg.width = [0.02, 0.02, 0.015, 0.04, 0.06, 0.07];       % choose default peak widths
buffer = 0.5;                                                       % a margin of the window for peak visualisation 
% ------------------------------------

check = questdlg('Do you want to track peaks ?', 'Peak tracking', 'YES', 'NO', 'NO'); 
if strcmp(check, 'YES')
    % check if numbers match
    if numel(seed_electrode) ~= numel(seed_peaks)
        disp('Number of seed electrodes does not correspond to the seed peak distribution!')
    end

    % select data to process
    index = ismember(subject_order, subject);
    AGSICI_data_subset = AGSICI_data(:, :, :, index, :, :);
    for p = 1:length(position)
        for c = 1:length(current)     
            for s = 1:size(AGSICI_data_subset, 4)
                % average across intensities
                for e = 1:size(AGSICI_data, 5)
                    for t = 1:size(AGSICI_data, 6)
                        eoi_data(p, c, s, e, t) = mean(squeeze(AGSICI_data_subset(p, c, :, s, e, t)));
                    end
                end 
            end        
        end
    end
    clear index p c s e t 

    % track selected peaks
    % a = 1; k = 1; p = 1; c = 1;    
    for k = 1:length(AGSICI_TEP_avg.peak)  
        % index seed electrode
        for a = 1:numel(seed_electrode)
            if any(seed_peaks{a} == k)
                seed = find(contains(labels, seed_electrode{a}));
            end
        end

        % check if this peek should be tracked
        answer = questdlg(['Do you want to track ' AGSICI_TEP_avg.peak{k} ' ?'], ...
            [seed_electrode{a} ' electrode , peak ' AGSICI_TEP_avg.peak{k}], 'YES', 'NO', 'NO'); 
        if strcmp(answer, 'NO')
            continue
        else
            % loop through the datasets
            for p = 1:length(position)
                for c = 1:length(current)
                    for s = 1:length(subject)              
                        % choose data 
                        for e = 1:size(eoi_data, 4)
                            data(1, e, 1, 1, 1, :) = squeeze(eoi_data(p, c, s, e, :));
                        end

                        % identify peak latency for current dataset
                        [peak_center, data_c, data_c_sub] = track_peak(data, header, time_window, ...
                            k, subject(s), AGSICI_TEP_avg, buffer, seed);

                        % fill in outcome structure
                        s_index = length(subject_order) - length(subject) + s;
                        AGSICI_TEP_subject(s_index).latency(p, c, k) = peak_center;

                        % append centered data
                        statement = ['AGSICI_' AGSICI_TEP_avg.peak{k} '_' position{p} '_' current{c} ...
                            '_centered(s, :, :) = data_c;'];
                        eval(statement)

                        % append centered subtracted data
                        statement = ['AGSICI_' AGSICI_TEP_avg.peak{k} '_' position{p} '_' current{c} ...
                            '_subtracted(s, :, :) = data_c_sub;'];
                        eval(statement)
                    end
                end
            end
        end
    end
    clear seed answer s_index a k p c e data statement peak_center data_c data_c_sub 

    % average across subjects, save for letswave
    for k = 1:length(AGSICI_TEP_avg.peak)    
        % fill in peak name 
        AGSICI_eoi(k).peak = AGSICI_TEP_avg.peak{k};

        % loop through datasets
        for p = 1:length(position)
            for c = 1:length(current)
                % skip round if the dataset doesn't exists
                data_name = ['AGSICI_' AGSICI_TEP_avg.peak{k} '_' position{p} '_' current{c} '_centered'];
                if ~exist(data_name)
                    disp(['Peak ' AGSICI_TEP_avg.peak{k} ' - ' position{p} ', ' current{c} ' : dataset not found.'])
                    continue
                end

                % fill in mean latency
                for s = 1:length(subject_order)
                    latency(s) = AGSICI_TEP_subject(s).latency(p, c, k);
                end
                AGSICI_eoi(k).latency(p, c) = mean(latency);

                % choose centered data
                statement = ['data_i = ' data_name ';'];
                eval(statement)

                % average across trials (subject x medication)
                for e = 1:size(data_i, 2)
                    for i = 1:size(data_i, 3)
                        data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
                        AGSICI_eoi(k).data_centred(p, c, e, i) = mean(squeeze(data_i(:, e, i)));        
                    end
                end

                % save for LW
                fname = ['AGSICI TEP tracked_2 ' position{p} ' ' current{c} ' '  AGSICI_TEP_avg.peak{k} ' centered'];
                header.name = fname; 
                header.datasize(6) = size(data, 6);
                span = (1 + buffer) * AGSICI_TEP_avg.width(k);
                header.xstart = - span/2;
                save([fname '.mat'], 'data');
                save([fname '.lw6'], 'header');    
                clear data

                % choose centered subtracted data
                data_name = ['AGSICI_' AGSICI_TEP_avg.peak{k} '_' position{p} '_' current{c} '_subtracted;'];
                statement = ['data_i = ' data_name ';'];
                eval(statement)

                % average across trials (subject x medication)
                for e = 1:size(data_i, 2)
                    for i = 1:size(data_i, 3)
                        data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
                        AGSICI_eoi(k).data_sub(p, c, e, i) = mean(squeeze(data_i(:, e, i)));        
                    end
                end

                % save for LW
                fname = ['AGSICI TEP tracked ' position{p} ' ' current{c} ' '  AGSICI_TEP_avg.peak{k} ' subtracted'];
                header.name = fname; 
                save([fname '.mat'], 'data');
                save([fname '.lw6'], 'header');    
                clear data  
            end
        end
    end
    clear k p c s e i data_name  statement fname span

    % append new variables to the general MATLAB file
    save([filename '.mat'], 'AGSICI_TEP_avg',  'AGSICI_TEP_subject', 'AGSICI_eoi', '-append');
    clear seed_electrode seed_peaks buffer eoi_n eoi_data 
else
    if exist('AG-SICI.mat') > 0
        load('AG-SICI.mat', 'AGSICI_TEP_subject')
        AGSICI_TEP_subject = rmfield(AGSICI_TEP_subject, {'latency_real', 'amplitude'});
        save([filename '.mat'], 'AGSICI_TEP_avg',  'AGSICI_TEP_subject', '-append');
    end
end

%% 6) GFP amplitude extraction
% ----- decide output parameters -----
POI = {'P25'};                                      % peaks of interest
percent = 20;                                       % % of timepoints included in the mean amplitude calculation
% ------------------------------------

% make sure that the dataset contains all tested participants
if size(AGSICI_data, 4) ~= length(subject_order) 
    error('%d subjects were found in the dataset, but %d were expected!', ...
        size(AGSICI_data, 4), length(subject_order))
end

% loop through subjects and conditions
% p = 1; c = 1; s = 1; i = 1; k = 1;
for s = 1:length(subject_order) 
    for p = 1:length(position)  
        for c = 1:length(current)
            % setup names   
            figure_title = ['Subject n. ' num2str(subject_order(s)) ' : ' position{p} ' STS, ' current{c} ' current'];                                       
            if subject_order(s) < 10
                file_name = ['AG-SICI_S0' num2str(subject_order(s)) '_' position{p} '_' current{c}];
            else
                file_name = ['AG-SICI_S' num2str(subject_order(s)) '_' position{p} '_' current{c}];
            end
            
            % launch summary figure 
            if figure_counter < 3
                figure_counter = 3;
            end
            fig = figure(figure_counter); 
                        
            % loop through intensities
            for i = 1:length(intensity)
                % loop through peaks
                for k = 1:length(POI)                      
                    % choose data 
                    data_visual = double(squeeze(AGSICI_GFP_subject.data(p, c, i, s, :))); 
                    
                    % choose data for topoplot 
                    for e = 1:size(AGSICI_data, 5)
                        data_topoplot(1, e, 1, 1, 1, :) = squeeze(AGSICI_data(p, c, i, s, e, :));
                    end

                    % define default TOI 
                    peak_n = find(contains(AGSICI_TEP_avg.peak, POI{k}));
                    center = AGSICI_TEP_avg.center(peak_n);
                    span = AGSICI_TEP_avg.width(peak_n);
                    
                    % set the final TOI
                    finish = 0;
                    while finish == 0;
                        % launch the figure
                        fig_1 = figure(1); 
                        hold on                      
                        
                        % plot the background 
                        subplot(4, 6, [7:24]); 
                        plot(x, data_visual, 'b:', 'LineWidth', 0.5)
                        yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
                        xlim([time_window(1), time_window(2)])
                        rectangle('Position', [0, yl(1), 0.01, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
                        title([figure_title ' : ' POI{k}], 'FontSize', 16)
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])

                        % set limits for topoplot colourmap
                        map_lims = yl;
                        
                        % visualize default peak TOI
                        subplot(4, 6, [7:24])
                        rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')
                        
                        % calculate mean amplitude
                        [amplitude, averaged_x, averaged_data] = GFP_amplitude(data_visual', center, span, percent, xstep, time_window(1)); 

                        % calculate  real peak latency (round up)
                        central_latency = averaged_x(ceil(length(averaged_x)/2));

                        % update the figure
                        subplot(4, 6, [7:24])
                        for a = 1:length(averaged_x)
                            line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)], 'Color', 'red', 'LineWidth', 1)
                            hold on
                        end
                        
                        % add the topoplot    
                        subplot(4, 6, 1);
                        topo_plot(header, data_topoplot, central_latency, time_window(1), map_lims) 
                        hold on 
                        
                        % replot the data to make it visible
                        subplot(4, 6, [7:24])
                        plot(x, data_visual, 'b', 'LineWidth', 2.5)
                        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
                        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
                        hold on
                        
                        % ask for approval
                        answer = questdlg('Do you want to proceed?', POI{k},...
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
                                choose_data = double(squeeze(AGSICI_GFP_subject.data(p, c, i, s, choose_x1 : choose_x2)));
                                choose_header = header;
                                choose_header.datasize(6) = length(choose_data);  
                                choose_header.xstart = choose_center - choose_span/2;
                                
                                % check if vector size matches
                                if length(choose_data) ~= length(choose_x)
                                    diff = length(choose_data) - length(choose_x);
                                    if diff > 0
                                        choose_data = choose_data(1:end - diff);
                                    elseif diff < 0
                                        choose_x = choose_x(1:end + diff);
                                    end
                                end

                                % launch the choosing figure                 
                                choose_figure_name = ['Choose manually peak ' POI{k}];
                                choose_axesHandles = [];
                                choose_fig = figure(2);   
                                choose_axesHandles = [choose_axesHandles subplot(3, 3, [4:9])];  
                                plot(choose_x, choose_data, 'LineWidth', 2)
                                xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                                title(choose_figure_name, 'FontSize', 16)
                                hold on                

                                % plot the line at the center
                                l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                                hold on    

                                % plot the central topography 
                                choose_map_lims = get(choose_axesHandles(1), 'ylim');
                                choose_axesHandles = [choose_axesHandles subplot(3, 3, 2)];
                                topo_plot(header, data_topoplot, choose_center, time_window(1), choose_map_lims);
                                hold on            

                                % choose the peak position
                                pos_x = get_position(choose_axesHandles);  

                                % update the figure
                                set (choose_fig, 'WindowButtonMotionFcn', '');
                                subplot(3, 3, [4:9])
                                set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                                subplot(3, 3, 2) 
                                cla(choose_axesHandles(2))
                                topo_plot(header, data_topoplot, pos_x, time_window(1), choose_map_lims);
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
                    clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l
                    
                    % record outcome variables
                    AGSICI_GFP_subject.latency(p, c, i, s, peak_n) = central_latency;  
                    AGSICI_GFP_subject.amplitude(p, c, i, s, peak_n) = amplitude; 
                                        
                    % set up figure
                    r = (i - 1)*6;
                    f = (1 + 3 * (k - 1) + r);
                    h = (2 + 3 * (k - 1) + r);  
                    g = (3 + 3 * (k - 1) + r); 
                                        
                    % plot the data
                    figure(fig)
                    hold on
                    
                    %plot the data
                    subplot(length(intensity), 6, [f, h])
                    plot(x, data_visual, 'b:', 'LineWidth', 0.5)
                    yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
                    xlim([time_window(1), time_window(2)])
                    rectangle('Position', [0, yl(1), 0.01, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
                    title([intensity{i}(end-2 : end) ' % rMT - ' POI{k}], 'FontSize', 16)
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                    
                    % set limits for topoplot colourmap
                    map_lims = yl;

                    % visualize final peak TOI
                    rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')
                    
                    % update the figure
                    subplot(length(intensity), 6, [f, h])
                    for a = 1:length(averaged_x)
                        line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)], 'Color', 'red', 'LineWidth', 1)
                        hold on
                    end

                    % add the topoplot    
                    subplot(length(intensity), 6, g);
                    topo_plot(header, data_topoplot, central_latency, time_window(1), map_lims);
                    hold on 

                    % replot the data to make it visible
                    subplot(length(intensity), 6, [f, h])
                    plot(x, data_visual, 'b', 'LineWidth', 2)
                    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2)
                    line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)                                     
                end                                       
            end
            
            % save the summary figure, close
            figure(fig)
            sgtitle(figure_title)  
            hold off
            savefig(fig, [pwd '\' filename '_figs\' file_name '.fig'])
            saveas(fig, [pwd '\' filename '_figs\' file_name '.png'])
            close(fig)

            % update the figure counter
            figure_counter = figure_counter + 1;   
        end
    end
    
    % play a celebratory sound at the end of each participant
    tune = load('handel.mat');
    sound(tune.y, tune.Fs)
end
clear p c i s k e a f g h r electrode_n peak_n figure_title file_name data_visual data_topoplot ...
    center span fig yl map_lims polarity amplitude averaged_x averaged_data central_latency...
    finish answer choose_fig choose_axesHandles choose_center choose_data choose_figure_name pos_x tune    

% save data in a R-compatible table 
if ~exist('AGSICI_GFP_amp')
    AGSICI_GFP_amp = table;
end
row_counter = height(AGSICI_GFP_amp) + 1;
for s = 1:length(subject_order) 
    for p = 1:length(position)  
        for c = 1:length(current)
            for i = 1:length(intensity)
                for k = 1:length(POI) 
                    %fill in the table
                    AGSICI_GFP_amp.subject(row_counter) = subject_order(s);
                    AGSICI_GFP_amp.position(row_counter) = position(p);
                    AGSICI_GFP_amp.current(row_counter) = current(c);
                    AGSICI_GFP_amp.intensity(row_counter) = str2double(intensity{i}(end-2 : end));
                    AGSICI_GFP_amp.peak(row_counter) = POI(k);
                    AGSICI_GFP_amp.amplitude(row_counter) = AGSICI_GFP_subject.amplitude(p, c, i, s, k);
                    AGSICI_GFP_amp.latency(row_counter) = AGSICI_GFP_subject.latency(p, c, i, s, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear s p c i k row_counter
writetable(AGSICI_GFP_amp, 'AGSICI_GFP_amp.csv')

% append new variables to the general MATLAB file
save([filename '.mat'], 'AGSICI_GFP_subject', 'AGSICI_GFP_amp', '-append');
clear POI percent

%% 7) TEP amplitude extraction
% ----- decide output parameters -----
AGSICI_TEP_default.peak = {'P25' 'N45'};                                    % peaks of interest
AGSICI_TEP_default.center = [0.025 0.045];                                  % default starting latencies
AGSICI_TEP_default.span = [0.015 0.015];                                    % span of the extraction window
AGSICI_TEP_default.eoi = {{'C1' 'C3' 'CP5'} {'F3' 'FC1' 'FC5'}};            % electrodes of interest
percent = 20;                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                                                          % y limits for topoplots 
% ------------------------------------
% set colours for visualisation
col_fig1 = [0.0902   0.3725    0.5608; 0.2549    0.8000    0.8000; 0.2549    0.8000    0.5451];

% loop through subjects and conditions
for s = 1
    % set axis counter
    axis_counter = 1;
    
    % loop through conditions            
    for p = 1:length(position)
        for c = 1:length(current)
            % setup names   
            figure_title = sprintf('Subject n. %d: %s STS, %s current', subject_order(s), position{p}, current{c});    
            
            % choose data for topoplot
            for i = 1:length(intensity)                    
                for e = 1:32
                    data_topoplot(i, e, 1, 1, 1, :) = squeeze(AGSICI_data(p, c, i, s, e, :));
                end
            end

            % launch summary figure 
            if figure_counter < 3
                figure_counter = 3;
            end
            fig = figure(figure_counter);

            % adjust figure size
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % loop through peaks
            for k = 1:length(AGSICI_TEP_default.peak)  
                % identify peak polarity
                if AGSICI_TEP_default.peak{k}(1) == 'P'
                    polarity = 1;
                else
                    polarity = -1;
                end

                % identify EOIs
                eoi = [];
                for e = 1:length(AGSICI_TEP_default.eoi{k})
                    eoi(e) = find(strcmp(labels, AGSICI_TEP_default.eoi{k}{e}));
                end

                % choose data
                data_visual = [];
                for i = 1:length(intensity)
                    % data for timeseries visualization
                    data_visual(i, :, :) = squeeze(AGSICI_data(p, c, i, s, eoi, :)); 
                end

                % average across eois
                if length(eoi) > 1
                    data_visual = squeeze(mean(data_visual, 2));
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
                    [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, time_window(1), polarity);

                    % update the figure
                    subplot(4, 6, 1:18)
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
                        subplot(4, 6, 18 + a);
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
                            choose_figure_name = ['Choose manually peak ' AGSICI_TEP_default.peak{k}];
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
                            for a = 1:size(data_visual, 1)
                                choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                                topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims);
                            end           

                            % choose the peak position
                            pos_x = get_position(choose_axesHandles);  

                            % update the figure
                            set (choose_fig, 'WindowButtonMotionFcn', '');
                            subplot(4, 4, [5:16])
                            set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                            for a = 1:size(data_visual, 1)
                                subplot(4, 4, a) 
                                cla(choose_axesHandles(2))
                                topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims);
                            end
                            hold off

                            % update the central latency
                            center = pos_x;

                            % close the choosing figure
                            pause(1)
                            close(choose_fig)

                            % close the the main figure
                            close(fig_1)
                    end
                end
                clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
                    choose data choose_center choose_axesHandles answer

                % record outcome variables
                for i = 1:length(intensity)
                    AGSICI_TEP_peaks.latency(p, c, i, s, k) = lat_peak(i); 
                    AGSICI_TEP_peaks.amplitude_peak(p, c, i, s, k) = y_max(i); 
                    AGSICI_TEP_peaks.amplitude_mean(p, c, i, s, k) = y_mean(i); 
                end

                % set up the main figure
                figure(fig)
                subplot(8, 6, [1 2 3] + 6*(axis_counter-1))
                hold on
                plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
                yl = get(gca, 'ylim'); 
                xlim(time_window);
                rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
                line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
                line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
                text(-0.14, 0, sprintf('%s', AGSICI_TEP_default.peak{k}), 'Color', colours((p-1)*2 + c, :), 'FontSize', 16, 'FontWeight', 'bold')
                set(gca, 'Fontsize', 10)
                ylabel('amplitude (\muV)')
                xlabel('time (s)') 

                % mark peaks 
                for a = 1:size(data_visual, 1)
                    line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', colours((p-1)*2 + c, :), 'LineWidth', 2)
                    plot(x, data_visual(a, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.75)
                    hold on
                end 
                
                % add title
                if mod(axis_counter, 2) == 1
                    set(get(gca, 'title'), 'string', sprintf('%s STS, %s current', position{p}, current{c}), 'FontSize', 14);
                end

                % add topoplots
                for a = 1:size(data_visual, 1)
                    subplot(8, 6, 3 + a + 6*(axis_counter-1))
                    topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims);   
                    if axis_counter == 1
                         set(get(gca, 'title'), 'string', sprintf('%s%% rMT', intensity{a}(end-2:end)), 'FontSize', 14);
                    end
                end

                % update axis counter
                axis_counter = axis_counter + 1;
            end                                       
        end
    end 
    % finalize the summary figure
    figure(fig)
    sgtitle(sprintf('subject n. %d', subject_order(s)))           
    hold off

    % name and save figure
    if subject_order(s) < 10
        subj = ['0' num2str(subject_order(s))];
    else
        subj = num2str(subject_order(s));
    end
    figure_name = ['AGSICI_TEP_' subj '_amplitude'];
    savefig([folder_figures '\TEP amplitude\' figure_name '.fig'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'AGSICI_TEP_peaks', '-append');
end
clear col_fig1 s p c figure_title fig axis_counter k polarity e i data_visual data_topoplot a subj figure_name

% save data in a R-compatible table 
if ~exist('AGSICI_TEP_peaks_table')
    AGSICI_TEP_peaks_table = table;
end
row_counter = height(AGSICI_TEP_peaks_table) + 1;
for s = 1:length(subject_order) 
    for p = 1:length(position)  
        for c = 1:length(current)
            for i = 1:length(intensity)
                for k = 1:length(AGSICI_TEP_default.peak) 
                    %fill in the table
                    AGSICI_TEP_peaks_table.subject(row_counter) = subject_order(s);
                    AGSICI_TEP_peaks_table.position(row_counter) = position(p);
                    AGSICI_TEP_peaks_table.current(row_counter) = current(c);
                    AGSICI_TEP_peaks_table.intensity(row_counter) = {intensity{i}(end - 2:end)};
                    AGSICI_TEP_peaks_table.peak(row_counter) = AGSICI_TEP_default.peak(k);
                    AGSICI_TEP_peaks_table.amplitude_peak(row_counter) = AGSICI_TEP_peaks.amplitude_peak(p, c, i, s, k);
                    AGSICI_TEP_peaks_table.amplitude_mean(row_counter) = AGSICI_TEP_peaks.amplitude_mean(p, c, i, s, k);
                    AGSICI_TEP_peaks_table.latency(row_counter) = AGSICI_TEP_peaks.latency(p, c, i, s, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m t s p k row_counter
writetable(AGSICI_TEP_peaks_table, 'AGSICI_TEP_peaks_table.csv')

% append new variables to the general MATLAB file
save(output_file, 'AGSICI_TEP_default', 'AGSICI_TEP_peaks', '-append');
clear percent map_lims

%% 8) plot mean TEP amplitude
% ----- peak amplitude -----
for k = 1:length(AGSICI_TEP_default.peak) 
    % launch the figure
    fig = figure(figure_counter); 
    hold on

    % plot the data
    counter = 1;
    for p = 1:length(position)
        for c = 1:length(current)
            % calculate mean and SEM
            for i = 1:length(intensity)            
                y(i) = mean(AGSICI_TEP_peaks.amplitude_peak(p, c, i, :, k));
                SEM(i) = std(AGSICI_TEP_peaks.amplitude_peak(p, c, i, :, k)) / sqrt(length(subject_order));
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
    
    % add features, adjust parameters
    set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
    set(gca, 'Fontsize', 14)
    title(sprintf('Peak amplitude: %s', AGSICI_TEP_default.peak{k}), ...
        'FontWeight', 'bold', 'FontSize', 16)
    xlabel('stimulation intensity (%rMT)'); ylabel('amplitude (\muV \pm SEM)');
    xlim([0.75, length(intensity) + 0.25])
    leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
    set(leg, 'Location','northwest', 'FontSize', 14);

    % save the figure       
    savefig([folder_figures '\AGSICI_amplitude_' AGSICI_TEP_default.peak{k} '.fig'])
    saveas(fig, [folder_figures '\AGSICI_amplitude_' AGSICI_TEP_default.peak{k} '.png'])

    % update the counter
    figure_counter = figure_counter + 1;  
end
clear p c i k fig counter y SEM perr leg

% ----- peak latency -----
for k = 1:length(AGSICI_TEP_default.peak) 
    % launch the figure
    fig = figure(figure_counter); 
    hold on

    % plot the data
    counter = 1;
    for p = 1:length(position)
        for c = 1:length(current)
            % calculate mean and SEM
            for i = 1:length(intensity)            
                y(i) = mean(AGSICI_TEP_peaks.latency(p, c, i, :, k));
                SEM(i) = std(AGSICI_TEP_peaks.latency(p, c, i, :, k)) / sqrt(length(subject_order));
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
    
    % add features, adjust parameters
    set(gca, 'xtick', 1:length(intensity), 'xticklabel', {'100' '120' '140'})
    set(gca, 'Fontsize', 14)
    title(sprintf('Peak latency: %s', AGSICI_TEP_default.peak{k}), ...
        'FontWeight', 'bold', 'FontSize', 16)
    xlabel('stimulation intensity (%rMT)'); ylabel('latency (ms \pm SEM)');
    xlim([0.75, length(intensity) + 0.25])
    leg = legend(perr, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'});
    set(leg, 'Location','northwest', 'FontSize', 14);

    % save the figure       
    savefig([folder_figures '\AGSICI_latency_' AGSICI_TEP_default.peak{k} '.fig'])
    saveas(fig, [folder_figures '\AGSICI_latency_' AGSICI_TEP_default.peak{k} '.png'])

    % update the counter
    figure_counter = figure_counter + 1;  
end
clear p c i k fig counter y SEM perr leg

%% 9) muscle contraction X N45: correlation per intensity - without outliers 
% ----- adjustable parameters -----
outliers = [5 8 12];
% ----- adjustable parameters -----

% calculate index to exclude outliers
index = subject ~= outliers(1) & subject ~= outliers(2) & subject ~= outliers(3);

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
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_TEP_peaks.amplitude_peak(p, c, i, index, 2))];
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
    savefig([folder_figures '\AGSICI_corr_TEP-N45Xcont_' intensity{i}(end-2:end) '_WO.fig'])
    saveas(fig, [folder_figures '\AGSICI_corr_TEP-N45Xcont_' intensity{i}(end-2:end) '_WO.png'])

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
    title(sprintf('Ranked data: %s %%rMT', intensity{i}(end-2:end)), 'FontWeight', 'bold', 'FontSize', 16)

    % save the figure       
    savefig([folder_figures '\AGSICI_corr_TEP-N45Xcont_' intensity{i}(end-2:end) '_WO_ranked.fig'])
    saveas(fig, [folder_figures '\AGSICI_corr_TEP-N45Xcont_' intensity{i}(end-2:end) '_WO_ranked.png'])

    % update the counter
    figure_counter = figure_counter + 1;  
end
clear data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp

%% 10) muscle contraction X N45: correlation per position/intensity - without outliers
% ----- adjustable parameters -----
outliers = [5 8 12];
% ----- adjustable parameters -----

% calculate index to exclude outliers
index = subject ~= outliers(1) & subject ~= outliers(2) & subject ~= outliers(3);

% plot the correlations
fig = figure(figure_counter);
for p = 1:length(position)  
    for i = 1:length(intensity)
        % extract data
        data_corr_m = [];
        data_corr_N45 = []; 
        for c = 1:length(current)
            data_corr_m = [data_corr_m; squeeze(AGSICI_muscle_activity.contraction.GFP_AUC(p, c, i, index))];
            data_corr_N45 = [data_corr_N45; squeeze(AGSICI_muscle_activity.N45.GFP_AUC(p, c, i, index))];
        end
        data_corr = [data_corr_m data_corr_N45];
        clear data_corr_muscle data_corr_N45

        % choose colours
        marker_col = [];
        for a = 1:size(data_corr, 1)
            marker_col = [marker_col; colours((p-1)*2 + 1, :)];
        end
        clear a 

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
        subplot(length(intensity), length(position), (i-1)*2 + p)
        hold on
        plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')
        title(sprintf('%s - %s', position{p}, intensity{i}), 'FontWeight', 'bold', 'FontSize', 14)
    end
end
clear p c data_corr data_corr_ranked fig data_model data_model_ranked marker_col temp
% 
% % save the figure       
% savefig([output_folder '\AGSICI_cont_corr_ranked_' position{p} '-' current{c} '_WO.fig'])
% saveas(fig, [output_folder '\AGSICI_cont_corr_ranked_' position{p} '-' current{c} '_WO.png'])
% 
% % update the counter
% figure_counter = figure_counter + 1;  

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
function [pos_x, data, sub_data] = track_peak(data, header, time_window, k, s, TEP, buffer, seed)
% figure params 
figure_name = ['Subject n. ' num2str(s) ' - peak ' TEP.peak{k}] ;
figure_center = TEP.center(k);
span = ((1 + buffer) * TEP.width(k));

% set the peak timepoint manually
finish = 0;
while finish == 0
    % identify the TOI of current peak                    
    x1 = ceil((figure_center - span/2 - time_window(1)) / header.xstep);
    x2 = ceil((figure_center + span/2 - time_window(1)) / header.xstep);
    x = (figure_center - span/2) : header.xstep : (figure_center + span/2);

    % prepare data and header for visualization
    data_visual = data(1, :, :, :, :, [x1 : x2]);
    header_visual = header;
    header_visual.datasize(6) = length(data_visual);  
    header_visual.xstart = figure_center - span/2;

    % launch the figure                    
    axesHandles = [];
    fig = figure;   
    hold on          
    
    % check the length
    if length(squeeze(data_visual(1, seed, :, :, :, :))) > numel(x)
        delta = length(squeeze(data_visual(1, seed, :, :, :, :))) -  numel(x);
        data_visual = data_visual(1, :, :, :, :, 1:end-delta);
    elseif length(squeeze(data_visual(1, seed, :, :, :, :))) < numel(x)
        delta = numel(x) - length(squeeze(data_visual(1, seed, :, :, :, :)));
        data_visual(1, :, :, :, :, end + 1:end + delta) = zeros(1, delta);
    end
    
    % plot windowed data
    axesHandles = [axesHandles subplot(3, 3, [4:9])];  
    plot(x, squeeze(data_visual(1, seed, :, :, :, :)), 'LineWidth', 2)
    xlim([(figure_center - span/2), (figure_center + span/2)])
    title(figure_name, 'FontSize', 16)

    % plot the line at the center
    h = line([figure_center, figure_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');   

    % plot the central topography 
    map_lims = get(axesHandles(1), 'ylim');
    axesHandles = [axesHandles subplot(3, 3, 2)];
    TEP_topoplot(header_visual, data_visual, figure_center, map_lims);

    % make the topography change with mouse movement 
    set(fig, 'WindowButtonMotionFcn', {@mouse_move, axesHandles, header_visual, data_visual});              

    % choose the peak position
    pos_x = get_position(axesHandles);             

    % update the figure
    set (fig, 'WindowButtonMotionFcn', '');
    subplot(3, 3, [4:9])
    set(h, 'XData', [pos_x, pos_x], 'YData', get(gca,'ylim'), 'LineStyle', '-');
    subplot(3, 3, 2) 
    cla(axesHandles(2))
    TEP_topoplot(header_visual, data_visual, pos_x, map_lims);

    pause(2)

    % ask for approval
    answer = questdlg('Do you want to proceed?', figure_name,...
        'Yes', 'No, I want to keep fiddling', 'Yes');

    % switch action
    switch answer
        case 'Yes'
            % close the figure
            close(fig)

            % exit the while loop
            finish = 1;

        case 'No, I want to keep fiddling'
            % assign previous center
            figure_center = pos_x;  

            % close the figure
            close(fig)                                
    end
end   

% index the peak timepoint  
n_peak = ceil((pos_x - time_window(1)) / header.xstep); 

% define final length
window_length = numel((pos_x - span/2) : header.xstep : (pos_x + span/2));

% loop through the data and subtract the value at the peak latency from all
% timepoints
sub_data = data;
topo_vector = squeeze(data(1, :, 1, 1, 1, n_peak));                                             
for i = 1:header.datasize(2)    % all channels
    for n = 1:length(data)      % all datapoints 
        sub_data(1, i, :, :, :, n) = data(1, i, :, :, :, n) - topo_vector(i);
    end
end

% crop the data
x_out1 = ceil((pos_x - span/2 - time_window(1)) / header.xstep);
x_out2 = ceil((pos_x + span/2 - time_window(1)) / header.xstep);
sub_data = squeeze(sub_data(1, :, :, :, :, [x_out1 : x_out2]));
data = squeeze(data(1, :, :, :, :, [x_out1 : x_out2]));

function TEP_topoplot(header, data, x_pos, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - header.xstart)/header.xstep);
vector = squeeze(data(1,:,1,1,1,x_visual));

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i = 1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i = i + 1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
function mouse_move(hObject,eventdata, axesHandles, header_visual, data_visual)

% get the position of the mouse
CP = get(hObject, 'CurrentPoint');
position = get(hObject, 'Position');
xCursor = CP(1,1)/position(1,3); % normalize
yCursor = CP(1,2)/position(1,4); % normalize

% get the position of the axes within the GUI
axesPos = get(axesHandles(1),'Position');
minx    = axesPos(1);
miny    = axesPos(2);
maxx    = minx + axesPos(3);
maxy    = miny + axesPos(4);

% check if the mouse is within the axes 
if xCursor >= minx && xCursor <= maxx && yCursor >= miny && yCursor <= maxy
    % get the cursor position within the lower axes
    currentPoint = get(axesHandles(1),'CurrentPoint');
    x = currentPoint(2,1);
    y = currentPoint(2,2);
    % update the topoplot in the uper axes
    map_lims = get(axesHandles(1), 'ylim');
    cla(axesHandles(2))
    subplot(3, 3, 2)
    TEP_topoplot(header_visual, data_visual, x, map_lims); 
    hold on
end
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
end
function [amplitude, averaged_x, averaged_data] = GFP_amplitude(data, center, span, percent, step, xstart)
% ------------------------------------------------------------------------
% Fnc: Calculates mean amplitude of the most prominent <percent> of the TOI 
%       --> TOI defined by the center value and the span of the time window
% 
% Input: 
%   data - vector of values (double), GFP 
%   center, span - num values that define the TOI
%   percent - how many top percent of datapoints will be included in the average
%   step, xstart - defines properties of time axes 
% 
% Author: Dominika (2020) 
% ------------------------------------------------------------------------

% prepare the interval vector 
interval = false(1, length(data));

% define the boundaries of the TOI - both limits are included in the interval!
start = round(((center-span/2) - xstart)/step);
stop = round(((center+span/2) - xstart)/step);
data_crop = data(start : stop);

% calculate number of points to average
points_number = ceil((percent/100) * length(data_crop));

% sort the data 
data_sorted = sort(data_crop, 'descend');        

% calculate the mean value
points_included = data_sorted(1:points_number);
amplitude = mean(points_included);

% index averaged points
for i = 1:points_number
    point_pos = find(data_crop == points_included(i)); 
    % in case there are two same values, take the earlier one
    index(i) = point_pos(1);
end
index = sort(index);

% fill ones in the interval vector at indexed positions
for i = 1:points_number
    interval(start + index(i) - 1) = true;
end

% calculate the outcome vectors for visualisation
averaged_data = data .* interval;                                  % keep only values ov averaged datapoints, the rest is set to 0
averaged_index = find(averaged_data);                               % index the averaged datapoints
averaged_x = averaged_index * step + xstart;                        % set the time interval
averaged_data = averaged_data(find(averaged_data));                 % get rid of the zeros

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
function plot_corr(data_model, data_corr, marker_col, corr_type)
% calculate correlation coefficient and p
[cor_coef, cor_p] = corr(data_corr, 'Type', corr_type);

% plot correlation
plot_cor = plotAdded(data_model);

% adjust parameters    
set(gca, 'FontSize', 14)
xlabel('muscular activity (GFP - AUC)'); ylabel('amplitude N45 (µV)');
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