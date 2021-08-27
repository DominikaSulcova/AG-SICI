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
% 
% 2) preliminary visualization of TEPs 
%       - averages data from all participants for all conditions, calculates CI 
%       - plots separately 4 coil placements to show differences in
%       intensities
%       - averages across intensities and plots mean + CI of all placement
%       conditions in one figure
% 
% 3) calculate global mean field power from grand average
%       - calculates GMFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in the outcome structure
% 
% 4) define peak widths using grand average GMFP - just informative 
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
% 6) extract peak amplitudes
%       - identifies amplitude of each timepoint from the window of
%       interest (can be manually adjusted around the peak)
%       - calculates mean amplitude of top 25% of these timepoints
%       - outputs as viusal representation (matlab figure and .png)
%       - saves amplitude and central latency of each peak to matlab var
%       - creates and saves a long-format table (.csv) with outcome data 
 

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
subject = [9, 11, 13:17, 20, 21];
subject_order = [1, 3:8, 10, 12, 9, 11, 13:17, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};
prefix = 'avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1'; 
filename = 'AG-SICI.mat';

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
save('colours.mat', 'colours'); clear a answer

% load a random header
load([prefix ' 01 across normal stim_100.lw6'], '-mat')
labels = {header.chanlocs.labels};

% visualization 
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

clear c p 
%% 1) load the data
% load data based on the prefix + conditions
% for p = 1:length(position)
%     for c = 1:length(current)
%         for i = 1:length(intensity)
%             for s = 1:length(subject)
%                 % load individual data from letswave
%                 if subject(s) < 10
%                     load([prefix ' 0' num2str(subject(s)) ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'])
%                 else
%                     load([prefix ' ' num2str(subject(s)) ' ' position{p} ' ' current{c} ' ' intensity{i} '.mat'])
%                 end
%                 
%                 % append the data in the data matrix
%                 A(p, c, i, s, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));
%             end
%         end
%     end 
% end
% disp(['Datasize: ' num2str(size(A))])
% clear p c i s data 
% 
% % save dataset to the global MATLAB file
% if exist(filename) == 0
%     save(filename, 'AGSICI_data');
% else
%     load(filename)
%     AGSICI_data = cat(4, AGSICI_data, A);    
%     save(filename, 'AGSICI_data', '-append');
% end

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
                gmfp_data(p, c, e, t) = mean(squeeze(AGSICI_data_mean(p, c, :, e, t)));
                gmfp_CI(p, c, e, t) = mean(squeeze(AGSICI_data_CI(p, c, :, e, t)));
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
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])

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
    plot(x, squeeze(gmfp_data(1, 1, find(contains(labels, electrode{e})), :)), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, y_limits(1), 0.01, y_limits(2) - y_limits(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
        
    % loop through datasets
    for p = 1:length(position)
        for c = 1:length(current) 
            % prepare data
            data_visual = squeeze(gmfp_data(p, c, find(contains(labels, electrode{e})), :))'; 
            CI_visual = squeeze(gmfp_CI(p, c, find(contains(labels, electrode{e})), :))'; 

            % plot data     
            P((p - 1)*2 + c) = plot(x, data_visual, 'Color', colours((p - 1)*2 + c, :), 'LineWidth', 2.5);
%             F((p - 1)*2 + c) = fill([x fliplr(x)],[data_visual + CI_visual fliplr(data_visual - CI_visual)], ...
%                 colours((p - 1)*2 + c, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        end
    end
    
    % add other parameters
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(y_limits)
    line([0, 0], y_limits, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
%     lgd = legend(P, {'along - normal' 'along - reversed' 'across - normal' 'across - reversed'}, 'Location', 'southwest');
%     lgd.FontSize = 14; title(lgd, 'Tested conditions')
    
    % name and save figure
    figure_name = ['AGSICI_TEP_all_' electrode{e}];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear e p c fig data_visual P F lgd figure_name CI_visual

% save dataset to the global MATLAB file
save(filename, 'AGSICI_data_mean', 'AGSICI_data_CI', 'gmfp_data', 'gmfp_CI', '-append');
clear electrode AGSICI_data_CI 

%% 3) GMFP
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------

% loop through conditions, compute GMFP, plot
row_count = 1;
for p = 1:length(position)
    for c = 1:length(current)
        % calculate GMFP (exclude target channel)
        AGSICI_GMFP(p, c, :) = std(squeeze(gmfp_data(p, c, 1:size(gmfp_data, 3) - 1, :)), 1);  

        % set dataset name + figure title
        fig_name = ['AG-SICI_GMFP_' position{p} '_' current{c}];
        fig_title = ['GMFP: ' position{p} ' STS, ' current{c} ' current'];

        % plot GMFP and extract peak latencies
        fig = figure(figure_counter);
        hold on

        if ~isempty(max_peaks)
        % plot GMFP
        h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
        AGSICI_TEP(row_count).latencies = gmfp_plot(x, squeeze(AGSICI_GMFP(p, c, :)), time_window, xstep, labeled, 'max_peaks', max_peaks);
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')

        % add topoplots
        for t = 1:length(AGSICI_TEP(row_count).latencies)
            % choose data for topoplot 
            for e = 1:size(AGSICI_data, 5)
                data_topoplot(1, e, 1, 1, 1, :) = squeeze(gmfp_data(p, c, e, :));
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
            AGSICI_TEP(row_count).latencies(:) = gmfp_plot(x, squeeze(AGSICI_GMFP(p, c, :)), time_window, xstep, labeled);
        end  
        hold off

        % save figure, update    
        savefig(fig_name)
        saveas(fig, [fig_name '.png'])
        figure_counter = figure_counter + 1;
        
        % update row counter
        row_count = row_count + 1;
    end
end
clear p c t h_axis pos fig_name fig_title row_count fig data_topoplot

% append new variables to the general MATLAB file
save(filename, 'AGSICI_GMFP', 'AGSICI_TEP', '-append')
clear labeled max_peaks gmfp_data gmfp_CI

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
        data = squeeze(AGSICI_GMFP(p, c, x_start_narrow:x_end_narrow));

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
        savefig(['AG-SICI_widths_' position{p} '_' current{c}])
        saveas(fig, ['AG-SICI_widths' position{p} '_' current{c} '.png'])
        figure_counter = figure_counter + 1;
        
        % update row counter
        row_count = row_count + 1;
    end
end

% append new variables to the general MATLAB file
save(filename, 'AGSICI_TEP', '-append');
clear x_limit x_start_narrow x_end_narrow x_narrow c p k data fig P L R row_count

%% 5) peak tracking
% ----- decide output parameters -----
seed_electrode = {'target' 'Cz'};                                   % electrode that will be used to set 0 timepoint
seed_peaks = {1:2; 3:6};                                            % which peeks use which seed electrode
AGSICI_TEP_avg.peak = {'P25' 'N40' 'P50' 'P75' 'N100' 'P180'};      % choose peak names
AGSICI_TEP_avg.center = [0.024, 0.038, 0.047, 0.075, 0.118, 0.200]; % choose default peak centers
AGSICI_TEP_avg.width = [0.02, 0.02, 0.015, 0.04, 0.06, 0.07];        % choose default peak widths
buffer = 0.5;                                                       % a margin of the window for peak visualisation 
% ------------------------------------
% loads previously saved variables
load(filename)

% check if numbers match
if numel(seed_electrode) ~= numel(seed_peaks)
    disp('Number of seed electrodes does not correspond to the seed peak distribution!')
end

% select data to process
index = ismember(subject_order, subject);
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:size(AGSICI_data(:, :, :, index, :, :), 4)
            % average across intensities
            for e = 1:size(AGSICI_data, 5)
                for t = 1:size(AGSICI_data, 6)
                    eoi_data(p, c, s, e, t) = mean(squeeze(AGSICI_data(p, c, :, s, e, t)));
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
                    
                    % in case there is already an average
                    if exist('AGSICI_eoi')
                    end
                end
            end
        end
    end
end
clear seed answer s_index a k p c e data statement peak_center data_c data_c_sub 

save()

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
                    AGSICI_eoi(k).data_2(p, c, e, i) = mean(squeeze(data_i(:, e, i)));        
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
                    AGSICI_eoi(k).data_2(e, i) = mean(squeeze((:, e, i)));        
                end
            end

            % save for LW
            fname = ['AGSICI TEP tracked_2 ' position{p} ' ' current{c} ' '  AGSICI_TEP_avg.peak{k} ' subtracted'];
            header.name = fname; 
            save([fname '.mat'], 'data');
            save([fname '.lw6'], 'header');    
            clear data  
        end
    end
end
clear k p c s e i data_name  statement fname span

% append new variables to the general MATLAB file
save(filename, 'AGSICI_TEP_avg',  'AGSICI_TEP_subject', 'AGSICI_eoi', '-append');
clear seed_electrode seed_peaks buffer eoi_n eoi_data 

%% 6) amplitude
% ----- decide output parameters -----
POI = {'P25' 'N40'};                                % 2 peaks of interest
EOI = {'target' 'target'};                          % 2 electrodes of interest
percent = 25;                                       % % of timepoints included in the mean amplitude calculation
% ------------------------------------

% loop through subjects and conditions
% p = 1; c = 1; s = 1; i = 1; k = 1;
for s = 1:length(subject) 
    for p = 1:length(position)  
        for c = 1:length(current)
            % setup names   
            figure_title = ['Subject n. ' num2str(subject(s)) ' : ' position{p} ' STS, ' current{c} ' current'];                                       
            if subject(s) < 10
                file_name = ['AG-SICI_amplitude_S0' num2str(subject(s)) '_' position{p} '_' current{c}];
            else
                file_name = ['AG-SICI_amplitude_S' num2str(subject(s)) '_' position{p} '_' current{c}];
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
                    % choose TEP data 
                    peak_n = find(contains(AGSICI_TEP_avg.peak, POI{k}));
                    electrode_n = find(contains(labels, EOI{peak_n}));
                    data_visual = double(squeeze(AGSICI_data(p, c, i, s, electrode_n, :))); 
                    
                    % choose data for topoplot 
                    for e = 1:size(AGSICI_data, 5)
                        data_topoplot(1, e, 1, 1, 1, :) = squeeze(AGSICI_data(p, c, i, s, e, :));
                    end

                    % define default TOI 
                    center = AGSICI_TEP_subject(s).latency(p, c, peak_n);
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

                        % define the polarity of current peak 
                        switch POI{k}(1)
                            case 'P'
                                polarity = 'positive';
                            case 'N'
                                polarity = 'negative';
                        end  
                        
                        % calculate mean amplitude
                        [amplitude, averaged_x, averaged_data] = TEP_amplitude(data_visual', polarity, center, span, percent, xstep, time_window(1)); 

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
                                choose_data = double(squeeze(AGSICI_data(p, c, i, s, electrode_n, choose_x1 : choose_x2)));
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
                    AGSICI_TEP_subject(s).latency_real(p, c, i, k) = central_latency;  
                    AGSICI_TEP_subject(s).amplitude(p, c, i, k) = amplitude; 
                                        
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
            savefig(fig, [file_name '.fig'])
            saveas(fig, [file_name '.png'])
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
AGSICI_outcome = table;
row_counter = 1;
for s = 1:length(subject) 
    for p = 1:length(position)  
        for c = 1:length(current)
            for i = 1:length(intensity)
                for k = 1:length(POI) 
                    %fill in the table
                    AGSICI_outcome.subject(row_counter) = subject(s);
                    AGSICI_outcome.position(row_counter) = position(p);
                    AGSICI_outcome.current(row_counter) = current(c);
                    AGSICI_outcome.intensity(row_counter) = str2double(intensity{i}(end-2 : end));
                    AGSICI_outcome.peak(row_counter) = POI(k);
                    AGSICI_outcome.amplitude(row_counter) = AGSICI_TEP_subject(s).amplitude(p, c, i, k);
                    AGSICI_outcome.latency(row_counter) = AGSICI_TEP_subject(s).latency_real(p, c, i, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear s p c i k row_counter
writetable(AGSICI_outcome, 'AGSICI_outcome.csv')

% append new variables to the general MATLAB file
save(filename, 'AGSICI_TEP_subject', 'AGSICI_outcome', '-append');
clear POI EOI percent

%% functions
function peak_x = gmfp_plot(x, y, time_window, xstep, labeled, varargin)
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
ylabel('power (\muV^2)')
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

topo_plot(vector2,chanlocs2,varargin{:});
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

function [amplitude, averaged_x, averaged_data] = TEP_amplitude(data, polarity, center, span, percent, step, xstart)
% ------------------------------------------------------------------------
% Fnc: Calculates mean amplitude of the most prominent <percent> of the TOI 
%       --> TOI defined by the center value and the span of the time window
% 
% Input: 
%   data - vector of values (double), signal from target electrode
%   polarity -  string, values: 'positive', 'negative' - according to TEP peak  
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

% sort the data according to peak polarity
switch polarity
    case 'positive'
        data_sorted = sort(data_crop, 'descend');        
    case 'negative'
        data_sorted = sort(data_crop); 
end

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

topo_plot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end