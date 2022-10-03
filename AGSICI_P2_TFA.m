%% AG-SICI: TMS-EVOKED POTENTIALS - time-frequency analysis
% Written by Dominika for AG-SICI project (2022)
% 
% Colection of scripts to perform TF decomposition of TEPs, visualize the
% outcome and run randomisation-based statistics on the results
% 
% Custom functions are included in the script, EEGLAB and letswave
% functions are being called from directories (paths need to be added)
% 
% Output:
%   --> figures are saved in a folder 'AG-SICI_P1_figures
% 
% 1) prepare single-trial TEP data
%       - loads data and saves them to a matlab structure (big!)
%       - downsamples and saves for letswave
%       - crops and exports to the EEGLAB .set format
% 
% 2) short-time time-frequency decomposition
%       - keeps only predefined channels {electrodes}
%       - computes STFT for both total and induced oscillations                 
%       - crops away epoch ends
%       - two-step normalization to z-score
%           (1) full-length single-trial normalization
%           (2) average across trials
%           (3) baseline correction
%       - saves for letswave
%
% 3) average across frequency bands
% 
% 
% 4) significant ERSP
%       - merges individual data to a single dataset per condition
%       - calculates t-test against 0 for baseline datasets
%           --> keeps only 'average' and C3 channels
%       - plots baseline group average data, marks significance
% 

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
folder_input = uigetdir(pwd, 'Choose the input folder');
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_figures = [folder_results '\AG-SICI_' study '_figures'];
output_file = [folder_results '\AG-SICI_' study '.mat'];

%% 1) single-trial TEP data
% ----- section input -----
prefix_old = 'bl icfilt-plus ica visual crop but fft-notchfilt prefilt prea P1'; 
placement = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity_old = {'stim_100' 'stim_120' 'stim_140'};
window_extract = [-0.5, 0.5];
suffix = '500Hz';
% -------------------------
% get new prefix
prefix = ['AGSICI ' study];
    
% loop through datasets
for s = 1:length(subject)
    cond_count = 1;
    for p = 1:length(placement)
        for c = 1:length(current)
            for i = 1:length(intensity_old)
                % identify dataset
                if subject(s) < 10
                    subj = ['0' num2str(subject(s))];
                else
                    subj = num2str(subject(s));
                end
                filename = sprintf('%s %s %s %s %s', prefix_old, subj, placement{p}, current{c}, intensity_old{i});

                % load dataset
                option = struct('filename', [folder_input '\' filename]);
                lwdata = FLW_load.get_lwdata(option);
                
                % crop raw data data & save to output variable
                x_start = (window_extract(1) - lwdata.header.xstart)/lwdata.header.xstep;
                x_end = (window_extract(2) - lwdata.header.xstart)/lwdata.header.xstep;
                data_save = lwdata.data(:, 1:length(electrodes), :, :, :, x_start:x_end);
                data_save = permute(single(data_save),[1,2,6,3,4,5]);
                AGSICI_data_individual(cond_count, i, s, 1:size(data_save, 1), 1:size(data_save, 2), 1:size(data_save, 3)) = data_save;
                
                % create new name
                lwdata.header.name = sprintf('%s %s %s %s', prefix, position{cond_count}, intensity{i}, subj);
                
                % downsample
                option = struct('x_dsratio', 4, 'suffix', suffix, 'is_save', 1);
                lwdata = FLW_downsample.get_lwdata(lwdata, option);
               
                % crop for EEGLAB
                x_start = (window_extract(1) - lwdata.header.xstart)/lwdata.header.xstep;
                x_end = (window_extract(2) - lwdata.header.xstart)/lwdata.header.xstep;
                lwdata.data = lwdata.data(:, 1:length(electrodes), :, :, :, x_start:x_end);
                data_extract = permute(single(lwdata.data),[2,6,1,3,4,5]);

                % create an EEGLAB structure & export
                load([folder_git '\EEGLAB_example.mat']);
                EEG = lw2eeglab(EEG, data_extract, lwdata.header);
                save(sprintf('%s\\AG-SICI_%s_export\\EEGLAB\\%s %s %s %s.set', folder_results, study, prefix, ...
                    position{cond_count}, intensity{i}, subj), 'EEG')
            end
            % update condition counter
            cond_count = cond_count + 1;
        end
    end
    % update 
    fprintf('Subject n.%d finished.\n', subject(s))
end

% save individual data 
save('C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Data\P1\AGSICI_data_individual.mat', 'AGSICI_data_individual', '-v7.3');

% update prefix
prefix = [suffix ' ' prefix];

% clean
clear s p c i subj filename lwdata data_save data_extract EEG x_start x_end AGSICI_data_individual ...
    cond_count prefix_old placement current intensity_old window_extract suffix 

%% 2) STFT
% ----- section input -----
oscillations = {'total' 'induced'};
hanning = 0.2;
slide = 1;
freqs = [1, 45, 1];     % [start, end, step]
output = 'power';
window_crop = [-0.5, 0.5];
baseline = [-0.3, -0.1];
suffix = {'stft_total' 'crop' 'full_z' 'avg' 'bl_z'};
% -------------------------
% loop through datasets
for s = 1:length(subject)
    for p = 1:length(position)
        for i = 1:length(intensity)
            % identify dataset
            if subject(s) < 10
                subj = ['0' num2str(subject(s))];
            else
                subj = num2str(subject(s));
            end
            filename = sprintf('%s %s %s %s', prefix, position{p}, intensity{i}, subj);
            fprintf('Subject n.%d: position %s, intensity %s%% rMT\n', subject(s), position{p}, intensity{i})

            % load the dataset
            option = struct('filename', [folder_input '\' filename]);                
            lwdata = FLW_load.get_lwdata(option);                

            % compute TF decomposition for both total and induced oscillations 
            for o = 1%:length(oscillations)
                % display the dataset info
                fprintf('%s oscillations:\n', oscillations{o})
                dataset = lwdata;
                dataset.header.name = [dataset.header.name ' ' oscillations{o}];

                % subtract evoked data if required
                if o == 2
                    % calculate mean evoked response 
                    evoked = mean(dataset.data, 1);
                    
                    % subtract
                    for c = 1:size(dataset.data, 1)
                        for e = 1:size(dataset.data, 2)
                            for t = 1:size(dataset.data, 6)
                                dataset.data(c, e, 1, 1, 1, t) = dataset.data(c, e, 1, 1, 1, t) - evoked(1, e, 1, 1, 1, t);
                            end
                        end
                    end
                end

                % SFFT 
                fprintf('...STFT')
                option = struct('hanning_width', hanning, 'sliding_step', slide, ...
                'low_frequency', freqs(1), 'high_frequency', freqs(2), ...
                'num_frequency_lines', (freqs(2) - freqs(1) + 1)/freqs(3), ...
                'output', output, 'show_progress', 0, 'suffix', suffix{1}, 'is_save', 0);
                dataset = FLW_STFT.get_lwdata(dataset, option);

                % crop data
                option = struct('xcrop_chk', 1, 'xstart', window_crop(1), 'xend', window_crop(2), ...
                    'suffix', suffix{2}, 'is_save', 0);
                dataset = FLW_crop_epochs.get_lwdata(dataset, option);

                % compute z-scores based on the full-length epoch
                fprintf('...full-length correction')
                option = struct('operation', 'zscore', 'xstart', window_crop(1), 'xend', window_crop(2), ...
                    'suffix', suffix{3}, 'is_save', 0);
                dataset = FLW_baseline.get_lwdata(dataset, option);

                % average across trials
                option = struct('operation', 'average', 'suffix', suffix{4}, 'is_save', 0);
                dataset = FLW_average_epochs.get_lwdata(dataset, option);

                % baseline correct
                fprintf('...baseline correction')
                option = struct('operation', 'zscore', 'xstart', baseline(1), 'xend', baseline(2),...
                    'suffix', suffix{5}, 'is_save', 1);
                dataset = FLW_baseline.get_lwdata(dataset, option);
                fprintf('...saved.\n')
            end
        end
    end
end

clear s p i o c e t subj filename option lwdata dataset hanning slide freqs ...
    output window_crop baseline evoked fig figure_name

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix 

%% 3) average across frequency bands
% ----- section input -----
fband.bands = {'delta' 'theta' 'alpha' 'beta1' 'beta2'};
fband.limits = {[1 4] [5 7] [8 12] [13 20] [21 30]};
suffix = {'merged_subj'};
% -------------------------
% loop through individual datasets 
for p = 3:length(position)
    for i = 1:length(intensity)
        % loop through frequency bands 
        for f = 1:length(fband.bands)
            % extract current fband, append subject to the dataset
            for s = 1:length(subject) 
                % identify dataset
                if subject(s) < 10
                    subj = ['0' num2str(subject(s))];
                else
                    subj = num2str(subject(s));
                end
                filename = sprintf('%s %s %s %s', prefix, position{p}, intensity{i}, subj);

                % load the dataset
                option = struct('filename', [folder_input '\' filename]);                
                lwdata = FLW_load.get_lwdata(option); 

                % create average channel
                lwdata.data(:, size(lwdata.data, 2) + 1, :, :, :, :) = mean(lwdata.data, 2);
                lwdata.header.datasize(2) = size(lwdata.data, 2);
                lwdata.header.chanlocs(size(lwdata.data, 2)).labels = 'average';                     

                % pool frequencies based on individual band limits
                lwdata_new.data(1, :, 1, 1, 1, :) = mean(lwdata.data(1, :, 1, 1, fband.limits{f}(1):fband.limits{f}(2), :), 5);
                        
                % update header
                lwdata_new.header = lwdata.header;
                lwdata_new.header.name = [fband.bands{f} ' ' lwdata_new.header.name];
                lwdata_new.header.datasize = size(lwdata_new.data);
                                                
                % append to the dataset
                lwdataset(s) = lwdata_new;                      
            end
                    
            % merge subjects & save to letswave - SAVES WITH '01'!
            option = struct('type', 'epoch', 'suffix', suffix, 'is_save', 1);
            lwdata = FLW_merge.get_lwdata(lwdataset, option); 
        end
    end
end
clear p i s f subj filename option lwdata lwdata_new lwdataset 

%% 4) test for significant ERSP & plot averages
channel = {'average'};
alpha = 0.05;
alpha_cl = 0.05;
n_perm = 1000;
axlim = {[-200, 500], [1, 45]};                     % limits of x (time) and y (frequency)
clim = [-2 2];                                      % limits of the colorscale {CS, TS}
suffix = {'merged_subj' 'signif'};
% -------------------------
for p = 1:length(position)
    for i = 1:length(intensity)
        % identify datasets
        for s = 1:length(subject) 
            if subject(s) < 10
                subj = ['0' num2str(subject(s))];
            else
                subj = num2str(subject(s));
            end
            filename{s} = sprintf('%s\\%s %s %s %s', folder_input, prefix, position{p}, intensity{i}, subj);
        end

        % load the datasets 
        option = struct('filename',{filename});
        lwdataset = FLW_load.get_lwdataset(option);

        % merge subjects & save to letswave - SAVES WITH '01'!
        option = struct('type', 'epoch', 'suffix', suffix{1}, 'is_save', 1);
        lwdata = FLW_merge.get_lwdata(lwdataset, option); 
                
        % create average channel
        lwdata.data(:, size(lwdata.data, 2) + 1, :, :, :, :) = mean(lwdata.data, 2);
        lwdata.header.datasize(2) = size(lwdata.data, 2);
        lwdata.header.chanlocs(size(lwdata.data, 2)).labels = 'average'; 
                    
        % identify channels(s) to plot
        labels = {lwdata.header.chanlocs.labels};
        for c = 1:length(channel)
            chanpos(c) = find(contains(labels, channel{c}) & cellfun('length', labels) == length(channel{c}));
        end
        
        % plot selected channel(s)
        for c = 1:length(chanpos)
            % plot TF map
            fig = figure(figure_counter);
            plot_map(lwdata, chanpos, 'axlim', axlim, 'clim', clim);
            figure_counter = figure_counter + 1; 
            
            % save the figure
            figure_name = sprintf('SFFT_%s_%s_induced', position{p}, intensity{i});
            savefig([folder_figures '\' figure_name '.fig'])
            saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        end
        
        % test for significant changes against 0
        option  = struct('constant', 0, 'tails', 'both', 'alpha', alpha, 'permutation', 1, 'cluster_threshold', alpha_cl, ...
            'num_permutations', n_perm, 'show_progress', 0, 'suffix', suffix{2}, 'is_save', 1);
        lwdata = FLW_ttest_constant.get_lwdata(lwdata, option);                        
    end
end
clear p i s c filename subj option lwdataset lwdata labels chanpos fig figure_name
clear channel window_crop alpha alpha_cl n_perm axlim clim

% update prefix
prefix = [suffix{1} ' ' prefix];
clear suffix 

%% test plot
% ----- section input -----
D = dataset;
channel = 'average';
axlim = {[-200, 500], [1, 45]};        
clim = [-5 5]; 
% -------------------------
fig = figure(figure_counter);
plot_map(D, channel, 'axlim', axlim, 'clim', clim);
figure_counter = figure_counter + 1; 
clear D channel axlim clim 

%% ) significant ERSPs
% ----- section input -----
window_crop = [0.01 0.5];
alpha = 0.05;
alpha_cl = 0.05;
n_perm = 1000;
suffix = 'signif';
% -------------------------
for s = 1:2
    for m = 1:length(medication)
        for t = 1:length(time)
            for o = 1:length(oscillations)
                for f = 1:length(fband.bands)
                    % identify dataset
                    filename = sprintf('merged_subj %s %s %s %s %s %s', ...
                        fband.bands{f}, prefix, medication{m}, time{t}, stimulus{s}, oscillations{o});
                    
                    % load dataset
                    option = struct('filename', [folder_input '\' filename]);
                    lwdata = FLW_load.get_lwdata(option);
                    
                    % crop time interval of interest
                    x_start = (window_crop(1) - lwdata.header.xstart)/lwdata.header.xstep;
                    x_end = (window_crop(2) - lwdata.header.xstart)/lwdata.header.xstep;
                    lwdata.data = lwdata.data(:, :, :, :, :, x_start:x_end);
                    
                    % update header
                    lwdata.header.datasize = size(lwdata.data);
                    lwdata.header.xstart = window_crop(1);
                    
                    % test for significant changes against 0
                    option  = struct('constant', 0, 'tails', 'both', 'alpha', alpha, 'permutation', 1, 'cluster_threshold', alpha_cl, ...
                        'num_permutations', n_perm, 'show_progress', 0, 'suffix', suffix, 'is_save', 1);
                    signif = FLW_ttest_constant.get_lwdata(lwdata, option);  
                end
            end
        end
    end
end
clear s m t o f filename option lwdata x_start x_end window_crop alpha alpha_cl n_perm 

% update prefix
prefix = [suffix ' ' prefix];
clear suffix 

%% save test plot
% ----- section input -----
figure_name = 'SFFT_';
% -------------------------
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
clear fig figure_name

%% functions
function EEG = lw2eeglab(EEG, data, header)
    % deal with datasize
    EEG.setname = header.name;
    EEG.filename = [];
    EEG.filepath = [];
    EEG.nbchan = size(data, 1);
    EEG.trials = size(data, 3);
    EEG.pnts = size(data, 2);
    EEG.srate = 1/header.xstep;
    EEG.times = header.xstart+(0:EEG.pnts-1)*header.xstep;
    EEG.xmin = EEG.times(1);
    EEG.xmax = EEG.times(end);    
    EEG.data = data; 
    
    % deal with channel locations
    EEG.chanlocs = header.chanlocs(1:size(data, 1));
    EEG.chanlocs = rmfield(EEG.chanlocs, {'SEEG_enabled' 'topo_enabled', 'sph_theta_besa', 'sph_phi_besa'});
    EEG.chanlocs = orderfields(EEG.chanlocs, [1, 2, 3, 6, 7, 8, 4, 5]);
    [EEG.chanlocs.sph_radius] = deal(85);
    [EEG.chanlocs.type] = deal('EEG');
    C = num2cell(1:length(EEG.chanlocs));
    [EEG.chanlocs.urchan] = C{:};
    [EEG.chanlocs.ref] = deal([]);
    EEG.urchanlocs = EEG.chanlocs;
    
    % deal with events
    EEG.event = header.events;
    if ~isempty(EEG.event)
        [EEG.event.type] = EEG.event.code;
        EEG.event = rmfield(EEG.event,'code');
        temp = num2cell([EEG.event.latency]/header.xstep);
        [EEG.event.latency] = deal(temp{:});     
        [EEG.event.urevent] = EEG.event.epoch;
        EEG.event = orderfields(EEG.event, [3, 1, 4, 2]);
        EEG.urevent = EEG.event;
    end
end
function plot_map(lwdata, channel, varargin)
    % prepare colormap 
    map = colorcet('D1','N', 201);
    
    % prepare data axes
    x = lwdata.header.xstart*1000 : lwdata.header.xstep*1000 : ((lwdata.header.datasize(6)-1)*lwdata.header.xstep + lwdata.header.xstart)*1000;
    y = lwdata.header.ystart : lwdata.header.ystep : (lwdata.header.datasize(5)-1)*lwdata.header.ystep + lwdata.header.ystart;
    
    % prepare data
    if strcmp(channel, 'average') & ~contains('average', {lwdata.header.chanlocs.labels})
        data = squeeze(mean(lwdata.data, [1, 2]));
    else
        data = squeeze(mean(lwdata.data(:, channel, :, :, :, :), 1));
    end
    
    % substitute with 0 where not significant
    s = find(strcmpi(varargin, 'signif'));
    if ~isempty(s)
        % prepare significance data
        signif = varargin{s + 1};
        map(101, :) = [1 1 1];
        for f = 1:size(data, 1)
            for i = 1:size(data, 2)
                if signif.data(1, channel, 3, 1, f, i)
                    data(f, i) = 0;
                end
            end
        end
    end
    
    % launch the figure
    set(gcf, 'units','centimeters','position',[10 10 15 10], 'color', 'w');
    hold on
    
    % plot the map
    imagesc(x,y,data)
    
    % set colorscale
    colorbar
    colormap(map)
    c = find(strcmpi(varargin, 'clim'));
    if ~isempty(c)
        clim = varargin{c + 1};
        set(gca, 'clim', clim)
    end 
    
    % set x and y limits
    a = find(strcmpi(varargin, 'axlim'));
    if ~isempty(a)
        axlim = varargin{a + 1};
        xlim(axlim{1})
        ylim(axlim{2})
    end
    
    % add TMS stimulus
    line([0, 0], axlim{2}, 'Color', [0 0 0], 'LineWidth', 2.5, 'LineStyle', '--')
    
    % other parameters
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
end


