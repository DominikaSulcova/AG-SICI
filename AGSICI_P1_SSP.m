%  processing with SSP-SIR

%% parameters
clear all; clc;

% dataset
subject = 13;
subjects = [1, 3:18, 20, 21];
sequence = {'along' 'across'; 'normal' 'reversed'};
block = 1:8;
orientation = {'along_normal' 'along_reversed' 'across_normal' 'across_reversed'};
intensity = {'100' '120' '140'};

% identify subject
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% choose relevant directories
folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');       % letswave + eeglab masterfiles
folder_output = uigetdir(pwd, 'Choose the output folder');         % processed data
cd(folder_output)

% load the finish sound
% load handel
load gong
soundwave = y; clear y Fs

%% import to letswave
% ----- section input -----
folders = [1:8; 1:8];
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% choose folder with raw NeurOne data
folder_input{1} = uigetdir('F:\AG-SICI\Raw data', 'S1: choose the input folder');
folder_input{2} = uigetdir('F:\AG-SICI\Raw data', 'S2: choose the input folder');

% import the datasets 
for s = 1:length(folder_input)
    fprintf('subject %d, session %d - %s\n', subject,  s,  sequence{1, s})
    for f = 1:size(folders, 2)        
        % determine current dirrection
        if mod(block(f), 2) == 1
            current = sequence{2, s};
        else
            current = sequence{2, ~ismember(sequence(2, :), sequence{2, s})};
        end
        
        % display current block
        fprintf('Block %d - %s\n', block(f), current)

        % import data
        [header, data] = EEG_import_MEGA(folder_input{s}, folders(s, f));

        % fill in letswave history entry
        header.history(1).configuration.parameters.input_folder  = folder_input{s};
        header.history(1).configuration.parameters.session_number  = folders(s, f);

        % remove extra 'Out' events
        disp('...reducing number of event categories')
        for i = 1:length(header.events)
            if ~strcmp(header.events(i).code, 'B - Out')
                index(i) = true;
            else
                index(i) = false;
            end
        end
        header.events = header.events(index);
        fprintf('...%d events found in the dataset\n', length(header.events))
        clear index

        % save dataset for letswave
        header.name = sprintf('AGSICI P1 S%s %s_%s b%d', subj, sequence{1, s}, current, block(f)); 
        CLW_save([], header, data);
    end
end
clear folder_input folders s f current header data i 
sound(soundwave)

%% preprocess 
% ----- section input -----
suffix = {'dc' 'ep' 'art-sup' 'ds'};
eventcode = 'B - Stimulation';
epoch = [-1 3];
interp = [-0.003, 0.003];
ds_ratio = 10;
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% cycle through sessions
for s = 1:2
    fprintf('subject %d, session %d - %s\n', subject,  s,  sequence{1, s})
    
    % cycle through blocks
    for b = block 
        fprintf('block %d:', b)

        % determine filename
        file = dir(['AGSICI P1 S' subj '*' sequence{1, s} '*b' num2str(b) '*.mat']);
        [~, filename, ~] = fileparts(file.name);

        % load data and header
        [header, data] = CLW_load(filename);
        
        % remove DC and apply linear detrend
        fprintf('...removing DC + detrending')
        [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);
        
        % segment 
        fprintf('...epoching')
        [header, data, ~] = RLW_segmentation(header, data, {{eventcode}}, 'x_start', epoch(1), 'x_duration', epoch(2));
        
        % interpolate TMS artifact
        fprintf('...interpolating the TMS artifact')
        [header, data, ~] = RLW_suppress_artifact_event(header, data,...
            'xstart', interp(1), 'xend',  interp(2), 'event_code', eventcode, 'interp_method', 'pchip');

        % downsample
        fprintf('...downsampling')
        [header, data, ~] = RLW_downsample(header, data, 'x_downsample_ratio', ds_ratio);
        
        % remove DC and apply linear detrend - per epoch
        fprintf('...detrending per epoch')
        [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);

        % save for letswave
        fprintf('...done.\n')
        header.name = [suffix{4} ' ' suffix{3} ' ' suffix{2} ' ' suffix{1} ' ' header.name];
        CLW_save([],header, data);
    end
end
fprintf('\n')

% create prefix
prefix = [suffix{4} ' ' suffix{3} ' ' suffix{2} ' ' suffix{1}];
clear suffix eventcode epoch interp ds_ratio s b file filename header data datasize 
sound(soundwave)

%% assign new eventcodes 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% choose folder with TMS sequence
folder_input = uigetdir('C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\AG-SICI\Subjects', 'TMS sequence: choose the input folder');

% load stim order matrices
load([folder_input '\stim_order_' subj '_N.mat'])
stim_order_N = outcome;
load([folder_input '\stim_order_' subj '_R.mat'])
stim_order_R = outcome;
clear outcome

% assign eventcodes
fprintf('assigning eventcodes:\n')
for s = 1:2
    fprintf('session %d - %s:', s, sequence{1, s})
    
    % identify stimulation order
    switch sequence{2, s}
        case 'normal'
            index = repmat('NR', 1, length(block)/2);
        case 'reversed'
            index = repmat('RN', 1, length(block)/2);
    end
    for b = block
        statement = ['stim_order(:, b) = stim_order_' index(b) '(:, b);'];
        eval(statement)
    end
    
    % rename eventcodes
    for b = block
        fprintf('block %d...', b)
        
        % determine filename
        file = dir([prefix ' AGSICI P1 S' subj '*' sequence{1, s} '*b' num2str(b) '*.mat']);
        [~, filename, ~] = fileparts(file.name);

        % load data and header
        [header, ~] = CLW_load(filename);
        
        % control for wrong number of events
        if length(header.events) < 75
            % ask for the right sequence
            fprintf('Found fewer than the expected number of events: %d', length(header.events))
            prompt = {'Use following TMS stimuli:'};
            dlgtitle = ['S' subj ' - block ' num2str(b)];
            dims = [1 50];
            definput = {'[1:75]'};
            stims2keep = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
            clear prompt dlgtitle dims definput  

            % replace the stimuli in the stim_order
            stim_order(1:length(stims2keep), b) = stim_order(stims2keep(1):stims2keep(end), b);
            
        elseif length(header.events) > 75
            fprintf('Found more than the expected number of events: %d', length(header.events))
        end
        
        % replace event codes
        for e = 1:length(header.events)
            header.events(e).code = num2str(stim_order(e, b));
        end
        
        % save dataset with the original name
        save([header.name '.lw6'], 'header');
    end
    fprintf('...done.\n')
end
fprintf('\n')

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));
clear folder_input stim_order_N stim_order_R stim_order s index b statement file filename header stims2keep e

%% parse by intensity
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% parse par eventcode & merge
fprintf('parsing by stimulation intensity:\n')
for o = 1:length(orientation)
    fprintf('%s...', orientation{o})
    
    % cycle through datasets
    files = dir([prefix ' AGSICI P1 S' subj '*' orientation{o} '*.mat']);
    for f = 1:length(files)
        % determine filenames
        [~, filename, ~] = fileparts(files(f).name);
        
      	% load the data
        option = struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);
        
        % split into datasets per eventcode
        for i = 1:length(intensity)
            % segment
            option = struct('event_labels', {intensity(i)}, 'x_start', -1, 'x_end', 2, 'x_duration', 3, ...
                'suffix', '', 'is_save', 0);
            lwdataset = FLW_segmentation_separate.get_lwdataset(lwdata, option); 
            
            % rename 
            lwdataset.header.name = sprintf('AGSICI P1 S%s %s %s', subj, orientation{o}, intensity{i});

            % save into a dataset for merging
            data2merge(i, f) = lwdataset;
        end
    end
    
    % merge epochs according to intensity
    for i = 1:length(intensity)
        % subset the dataset
        lwdataset = data2merge(i, :)';
        
        % merge subjects & save 
        option = struct('type', 'epoch', 'suffix', '', 'is_save', 1);
        lwdata = FLW_merge.get_lwdata(lwdataset, option); 
    end
end
fprintf('...done.\n')
fprintf('\n')
clear o files f filename option lwdata lwdataset data2merge i 

%% visual inspection, bad channels removal
% ----- section input -----
session = 'along'; %'along' 'across'
current = {'normal' 'reversed'};
prefix = 'visual';
epoch2remove = {[], ...
    [], ...
    []; ...
    [], ...
    [], ...
    []};
% -------------------------
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% cycle through current dirrections
fprintf('removing bad trials - %s STS:\n', session)
for c = 1:length(current)
    for i = 1:length(intensity)
        fprintf('%s current - %s %%rMT: ', current{c}, intensity{i})
        % determine the filename
        filename = sprintf('%s\\AGSICI P1 S%s %s_%s %s.lw6', folder_output, subj, session, current{c}, intensity{i});

        % load the data
        option = struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);
        
        % indentify epochs to keep
        epochs = 1:lwdata.header.datasize(1);
        epochs = epochs(~ismember(epochs, epoch2remove{c, i}));
        for e = 1:length(epochs)
            epoch2keep{e} = num2str(epochs(e));
        end
        
        % remove selected epochs
        option = struct('type', 'epoch', 'items', {epoch2keep}, 'suffix', prefix, 'is_save', 1);
        lwdata = FLW_selection.get_lwdata(lwdata, option);
        clear epoch2keep
        
        % indent if necessary
        if length(epoch2remove{c, i}) == 0;
            fprintf('No epochs removed.\n')
        end
    end
end
sound(soundwave)
fprintf('\n')
clear session epoch2remove c i filename epochs epoch2keep e option lwdata

%% export for EEGLAB
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export in .set format
fprintf('exporting to EEGLAB:\n')
for o = 1:length(orientation)
    fprintf('%s - ', strrep(orientation{o}, '_', ' '))
    for i = 1:length(intensity)
        fprintf('%s ... ', intensity{i})    

        % load the data
        name = sprintf('%s AGSICI P1 S%s %s %s', prefix, subj, orientation{o}, intensity{i});
        option = struct('filename', name);
        lwdata = FLW_load.get_lwdata(option);

        % export in .set format        
        export_EEGLAB(lwdata, name, subj);
    end
        fprintf('\n')
end
fprintf('done.\n')
fprintf('\n')

clear ref o i name option lwdata

%% bad channel inspection & SSP-SIR
% ----- section input -----
suffix = 'sspsir';
time_range = [-0.005, 0.3];
baseline = [-0.25 -0.005];
% ------------------------- 
% add eeglab and fastica to the of search path
addpath(fullfile(folder_toolbox,'eeglab2022.1'));
addpath(fullfile(folder_toolbox,'FastICA_25'));

% launch eeglab and generate an empty EEGLAB structure
fprintf('launching EEGLAB:\n')
eeglab 

% apply SSPSIR individually to each dataset and save the filters
counter = 1;
for o = 1:length(orientation)
    for i = 1:length(intensity)  
        % load the dataset
        name = sprintf('%s AGSICI P1 S%s %s %s.set', prefix, subj, orientation{o}, intensity{i}); 
        EEG = pop_loadset('filename', name, 'filepath', folder_output);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
        eeglab redraw  
        
        % add channel locations 
        EEG = pop_chanedit(EEG, 'lookup', [folder_toolbox '\eeglab2022.1\plugins\dipfit\standard_BEM\elec\standard_1005.elc']);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
        eeglab redraw;

        % visualize to identify possible bad channels
        figure; 
        pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);
        sgtitle(sprintf('S%s - %s, %s %%rMT', subj, strrep(orientation{o}, '_', ' '), intensity{i}))
        interp = inputdlg('Which channels do you want to interpolate?', 'Bad channels', [1 35], {'none'});
            
        % remove bad channels if necessary  
        interp = split(interp, ' ');
        if strcmp(interp{1}, 'none')
        else
            % identify electrodes to interpolate
            for t = 1:length(interp)
                for e = 1:length(EEG.chanlocs)
                    if strcmp(EEG.chanlocs(e).labels, interp{t})
                        interp_e(t) = e;
                    end
                end
            end
            EEG = pop_interp(EEG, interp_e, 'spherical'); 
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
            eeglab redraw;
        end

        % re-reference to common average & baseline correct
        EEG = pop_reref(EEG, []);
        EEG = pop_rmbase(EEG, baseline, []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, counter,'overwrite','on','gui','off'); 
        eeglab redraw;      

        % visual check - before SSP-SIR
        figure; 
        pop_timtopo(EEG, [-50  100], [5 10  25 45 75]);
        sgtitle(sprintf('S%s - %s, %s %%rMT', subj, strrep(orientation{o}, '_', ' '), intensity{i}))
        
        % ask if the original ssp filter should be applied
        answer{counter} = questdlg('Do you want to apply original SSP-SIR', 'SSP-SIR', 'YES', 'NO', 'YES'); 

        % SSP-SIR - spherical model 
        name = sprintf('%s AGSICI P1 S%s %s %s %s', prefix, subj, orientation{o}, intensity{i}, suffix); 
        [EEG, EEG_old] = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', time_range, 'PC', []);
        switch answer{counter}
            case 'YES'
                clear EEG_old
            case 'NO'
                EEG = EEG_old;
                clear EEG_old
        end    

        % baseline correct and save
        EEG = pop_rmbase(EEG, baseline, []);
        name = sprintf('%s AGSICI P1 S%s %s %s %s.set', prefix, subj, orientation{o}, intensity{i}, suffix); 
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, counter, 'setname', name, 'overwrite', 'on', 'gui', 'off');
        EEG.filename = name;
        eeglab redraw

        % visual check - after SSP-SIR
        figure; 
        pop_timtopo(EEG, [-50  100], [5 10  25 45 75]);
        sgtitle(sprintf('S%s - %s, %s %%rMT - FIRST FILTER', subj, strrep(orientation{o}, '_', ' '), intensity{i}))
        
        % save SSP-SIR parameters
        param(counter) = EEG.SSPSIR;

        % save dataset
        name = sprintf('%s AGSICI P1 S%s %s %s %s.set', prefix, subj, orientation{o}, intensity{i}, suffix);  
        pop_saveset(EEG, 'filename', name, 'filepath', folder_output);

        % update counter
        counter = counter + 1;
    end
end

% close all figures if required
answer_fig = questdlg('Do you want to close the figures?', 'figures', 'YES', 'NO', 'YES'); 
if strcmp(answer_fig, 'YES')
    figures = findall(0, 'Type', 'figure');
    for f = 1:length(figures)
        if strcmp(figures(f).Name, ' timtopo()')
            close(figures(f))
        elseif strcmp(figures(f).Name, '')
            close(figures(f))
        end
    end
end

% % apply ssp-sir filters of other compared datasets
% counter = 1;
% for o = 1:length(orientation)
%     for i = 1:length(intensity)  
%         % select the dataset
%         EEG = ALLEEG(counter);
%         [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
%         eeglab redraw
% 
%         % SSP-SIR 
%         other_datasets = [1:12];
%         other_datasets = other_datasets(strcmp(answer, 'YES'));
%         other_datasets = other_datasets(other_datasets ~= counter);
%         for a = other_datasets
%             [EEG_trials] = SSP_SIR_trials(EEG, param(counter).L_ave, ...
%                 param(a).topo, ...
%                 param(counter).kernel, []);
%             EEG.data =  EEG_trials.data;
%         end
% 
%         % baseline correct and save
%         EEG = pop_rmbase(EEG, baseline, []);
%         name = sprintf('%s AGSICI P1 S%s %s %s %s all_filt', prefix, subj, orientation{o}, intensity{i}, suffix); 
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, counter, 'setname', name, 'overwrite', 'on', 'gui', 'off'); 
%         eeglab redraw
% 
%         % visual check - after SSP-SIR
%         figure;
%         pop_timtopo(EEG, [-50  100], [5 10  25 45 75]);
%         sgtitle(sprintf('S%s - %s, %s %%rMT - FINAL', subj, strrep(orientation{o}, '_', ' '), intensity{i}))
% 
%         % save dataset
%         name = sprintf('%s AGSICI P1 S%s %s %s %s all_filt.set', prefix, subj, orientation{o}, intensity{i}, suffix);  
%         pop_saveset(EEG, 'filename', name, 'filepath', folder_output);
% 
%         % update counter
%         counter = counter + 1;
%     end
% end
clear suffix time_range baseline o i a counter name e t interp interp_e other_datasets param answer_fig f figures

%% export back to letswave
% ----- section input -----
suffix = 'sspsir';
% ------------------------- 
% add letswave 6 to the top of search path
addpath(genpath([folder_toolbox '\letswave6-master']));

% export data and header
counter = 1;
fprintf('exporting back to letswave: \n')
for o = 1:length(orientation)
    fprintf('%s: ', strrep(orientation{o}, '_', ' '))
    for i = 1:length(intensity)  
        fprintf('...%s ', intensity{i})
        % depending on whether EEGLAB is running:
        if ~isempty(findobj('Tag', 'EEGLAB'))
            if length(ALLEEG) >= counter
                % select the dataset in EEGLAB
                EEG = ALLEEG(counter);
                [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
                eeglab redraw   
            else
                % load the dataset
                name = sprintf('%s AGSICI P1 S%s %s %s %s.set', prefix, subj, orientation{o}, intensity{i}, suffix); 
                EEG = pop_loadset('filename', name, 'filepath', folder_output);
                [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
                eeglab redraw  
            end        
        else
            % add eeglab and fastica to the of search path
            addpath(fullfile(folder_toolbox,'eeglab2022.1'));
            addpath(fullfile(folder_toolbox,'FastICA_25'));
            
            % launch eeglab and generate an empty EEGLAB structure
            eeglab 
            
            % load the dataset
            name = sprintf('%s AGSICI P1 S%s %s %s %s.set', prefix, subj, orientation{o}, intensity{i}, suffix); 
            EEG = pop_loadset('filename', name, 'filepath', folder_output);
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, counter);
            eeglab redraw              
        end

        % load the original data
        [lwdata.header, lwdata.data] = CLW_load(sprintf('%s AGSICI P1 S%s %s %s', prefix, subj, orientation{o}, intensity{i}));
        
        % replace the data
        lwdata.data = [];
        for t = 1:size(EEG.data, 3)
            for e = 1:size(EEG.data, 1)
                for k = 1:size(EEG.data, 2)
                    lwdata.data(t, e, 1, 1, 1, k) = EEG.data(e, k, t);
                end
            end
        end
        
        % modify header
        lwdata.header.name = sprintf('%s %s AGSICI P1 S%s %s %s', suffix, prefix, subj, orientation{o}, intensity{i});
        lwdata.header.datasize = size(lwdata.data);
        lwdata.header.chanlocs = lwdata.header.chanlocs(1:size(lwdata.data, 2));
        lwdata.header.events = lwdata.header.events(1:size(lwdata.data, 1));
        
        % ssp-sir 
        lwdata.header.SSPSIR = EEG.SSPSIR;   
        
        % save 
        CLW_save([], lwdata.header, lwdata.data);

        % update counter
        counter = counter + 1;
    end
    fprintf('\n')
end
sound(soundwave)

% update prefix
prefix = [suffix  ' ' prefix];
clear name counter suffix ref o i t e k lwdata option ...
    ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG globalvars LASTCOM PLUGINLIST STUDY

%% preprocessing before ICA
% ----- section input -----
suffix = {'bandpass' 'notch' 'crop' 'bl' 'interp2'};
interp = [-0.005, 0.01];
bandpass = [0.1 80];
notch = 50;
crop = [-0.3 0.5];
baseline = [-0.25 -0.005];
% ------------------------- 
% cycle through datasets
for o = 1:length(orientation)
    fprintf('%s - ', strrep(orientation{o}, '_', ' '))
    for i = 1:length(intensity)
        fprintf('%s %%rMT:\n', intensity{i}) 
        % add letswave 7 to the top of search path
        addpath(genpath([folder_toolbox '\letswave7-master']));

        % load the data
        filename = sprintf('%s AGSICI P1 S%s %s %s', prefix, subj, orientation{o}, intensity{i});
        option = struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);

        % bandpass
        fprintf('...applying Butterworth bandpass filter')
        option = struct('filter_type', 'bandpass', 'high_cutoff', bandpass(2),'low_cutoff', bandpass(1),...
            'filter_order', 4, 'suffix', suffix{1}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);

        % notch
        fprintf('...applying FFT %dHz notch filter', notch)
        option = struct('filter_type', 'notch', 'notch_fre', notch, 'notch_width', 2, 'slope_width', 2,...
            'harmonic_num', 2, 'suffix', suffix{2},'is_save', 0);
        lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);

        % crop
        fprintf('...cropping')
        option = struct('xcrop_chk', 1, 'xstart', crop(1), 'xend', crop(2), 'suffix', suffix{3}, 'is_save', 0);
        lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);

        % baseline correction
        fprintf('...correcting for baseline\n')
        option = struct('operation','substract', 'xstart', baseline(1), 'xend', baseline(2), 'suffix', suffix{4}, 'is_save', 0);
        lwdata = FLW_baseline.get_lwdata(lwdata, option);

        % add letswave 6 to the top of search path
        addpath(genpath([folder_toolbox '\letswave6-master']));

        % interpolation
        fprintf('...interpolating\n')
        [header, data, ~] = RLW_suppress_artifact_event(lwdata.header, lwdata.data,...
            'xstart', interp(1), 'xend',  interp(2), 'interp_method', 'pchip');

        % save
        header.name = [suffix{5} ' ' header.name];
        CLW_save([], header, data);
    end
end
fprintf('done.\n')
fprintf('\n')
sound(soundwave)

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix interp bandpass notch crop baseline o i s filename option lwdata header data

 %% calculae ICA matrix 
% ----- section input -----
suffix = {'ica' 'sp_filter'};
% -------------------------
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% calculate ICA matrix separately for each session (along / across)
for s = 1:size(sequence, 1)
    fprintf('calculating ICA matrix: orientation %s STS\n', sequence{1, s}) 
    % identify the files
    filenames = dir(sprintf('%s AGSICI P1 S%s %s*.mat', prefix, subj, sequence{1, s}));
    for d = 1:length(filenames) 
        [~, filename, ~] = fileparts(filenames(d).name);
        dataset{d} = sprintf('%s\\%s.lw6', folder_output, filename); 
    end

    % load the dataset
    option = struct('filename', {dataset});
    lwdataset = FLW_load.get_lwdataset(option);

    % compute the ICA matrix
    option = struct('ICA_mode', 3, 'algorithm', 1, 'percentage_PICA', 100, 'criterion_PICA', 'LAP', 'suffix', suffix{1}, 'is_save', 1);
    lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    fprintf('\n')
end

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix counter s filenames filename d dataset option lwdataset

%% baseline correct and average
% ----- section input -----
suffix = {'bl' 'avg'};
baseline = [-0.25 -0.005];
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% cycle through datasets
for o = 1:length(orientation)
    fprintf('%s - ', strrep(orientation{o}, '_', ' '))
    for i = 1:length(intensity)
        fprintf('%s %%rMT:\n', intensity{i}) 
        % add letswave 7 to the top of search path
        addpath(genpath([folder_toolbox '\letswave7-master']));

        % load the data
        filename = sprintf('%s AGSICI P1 S%s %s %s', prefix, subj, orientation{o}, intensity{i});
        option = struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);

        % baseline correction
        fprintf('...correcting for baseline\n')
        option = struct('operation','substract', 'xstart', baseline(1), 'xend', baseline(2), 'suffix', suffix{1}, 'is_save', 0);
        lwdata = FLW_baseline.get_lwdata(lwdata, option);

        % average across epochs & save
        option = struct('operation', 'average', 'suffix', suffix{2}, 'is_save', 1);
        lwdata = FLW_average_epochs.get_lwdata(lwdata, option);
    end
end
fprintf('\n')

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear suffix baseline o i s option lwdata

%% export for Ragu
% ----- section input -----
time_window = [-0.05, 0.3];
% -------------------------
% calculate crop limits
load(sprintf('%s\\%s AGSICI P1 S%s %s %s.lw6', folder_output, prefix, subj, orientation{1}, intensity{1}), '-mat')                    
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;

% write text files for Ragu 
fprintf('exporting for Ragu:\n')
for o = 1:length(orientation)
    fprintf('%s ', strrep(orientation{o}, '_', ' '))
    for i = 1:length(intensity)
        fprintf('%s %%rMT... ', intensity{i}) 
        % load the data
        name_old = sprintf('%s\\%s AGSICI P1 S%s %s %s.mat', folder_output, prefix, subj, orientation{o}, intensity{i});
        load(name_old)
        data = squeeze(data(:, 1:32, :, :, :, x_start:x_end))';

        % save as .csv               
        name = sprintf('AGSICI P1 S%s %s %s.csv', subj, orientation{o}, intensity{i});
        writematrix(data, sprintf('%s\\export\\%s', folder_output, name));
    end
    fprintf('\n')
end
fprintf('done.\n')
fprintf('\n')

% create the montage file
% name = sprintf('%s\\export\\AGSICI_montage.xyz', folder_output);
% fileID = fopen(name, 'a');
% fprintf(fileID, '32\r\n');
% for a = 1:32
%     fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
%         header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
% end
% fclose(fileID);
clear o i name_old data x_start x_end name time_window fileID a header   

%% save individual data to the output structure
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% results folder
folder_results = uigetdir(pwd, 'Choose the Results folder');
output_file = [folder_results '\AG-SICI_P1_SSP.mat'];

% prepare empty structure
data_individual = [];
save(output_file, 'data_individual')

% loop through subjects
for s = 1:length(subjects)
    % identify subject
    if subjects(s) < 10
       subj = ['0' num2str(subjects(s))];
    else
       subj = num2str(subjects(s)); 
    end
    fprintf('subject %s ', subj)

    % prepare data
    for o = 1:length(orientation)    
        for i = 1:length(intensity)
            % load 
            option = struct('filename', sprintf('%s\\%s AGSICI P1 S%s %s %s.lw6', folder_output, prefix, subj, orientation{o}, intensity{i}));
            lwdata = FLW_load.get_lwdata(option);

            % crop
            option = struct('xcrop_chk', 1, 'xstart', time_window(1) - lwdata.header.xstep, 'xend', time_window(2), 'suffix', '', 'is_save', 0);
            lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);

            % save default header
            if s == 1 && o == 1 && i == 1
                header = lwdata.header;
                header.name = 'template';
                save(output_file, 'header', '-append');
            end
    
            % append to the data variable
            data_individual(o, i, s, :, :) = squeeze(lwdata.data(1, :, 1, 1, 1, :));
            save(output_file, 'data_individual', '-append');
        end
    end    
end
counter = 1;
for o = 1:2   
    for c = 1:2
        for i = 1:length(intensity)
            data_individual(counter, i, 16, :, :) = squeeze(AGSICI_data(o, c, i, 16, 1:32, :));
            save(output_file, 'data_individual', '-append');
        end
        counter = counter + 1;
    end
end
data_individual(o, i, 16, 1, 1:20)

%% muscle artifact extraction
% ----- section input -----
prefix = 'visual';
suffix = {'interp' 'reref' 'bl' 'muscles crop'};
session = {'along' 'across'};
time_range = [-0.05, 0.3];
baseline = [-0.25 -0.005];
electrodes = 1:2:31;
gfp_toi = [3, 10];
x_step = 0.5;
% ------------------------- 
% add eeglab and fastica to the of search path
addpath(fullfile(folder_toolbox,'eeglab2022.1'));
addpath(fullfile(folder_toolbox,'FastICA_25'));

% launch eeglab and generate an empty EEGLAB structure
fprintf('launching EEGLAB:\n')
eeglab 

% identify bad channels to interpolate
for s = 1:length(subjects)
    % identify subject
    if subjects(s) < 10
       subj = ['0' num2str(subjects(s))];
    else
       subj = num2str(subjects(s)); 
    end
    fprintf('subject %s: ', subj)

    for o = 1:length(session)
        fprintf('%s STS ... ', session{o})

        % load the dataset
        name = sprintf('%s AGSICI P1 S%s %s_normal %s.set', prefix, subj, session{o}, intensity{1}); 
        EEG = pop_loadset('filename', name, 'filepath', folder_output);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw  
        
        % add channel locations 
        EEG = pop_chanedit(EEG, 'lookup', [folder_toolbox '\eeglab2022.1\plugins\dipfit\standard_BEM\elec\standard_1005.elc']);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;

        % visualize to identify possible bad channels
        figure; 
        pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);
        sgtitle(sprintf('S%s - %s, %s %%rMT', subj, strrep(session{o}, '_', ' '), intensity{1}))
        interp = inputdlg('Which channels do you want to interpolate?', 'Bad channels', [1 35], {'none'});
            
        % identify electrodes to interpolate
        interp = split(interp, ' ');
        if strcmp(interp{1}, 'none')
            chan2interp{s, o} = [];
        else            
            for t = 1:length(interp)
                for e = 1:length(EEG.chanlocs)
                    if strcmp(EEG.chanlocs(e).labels, interp{t})
                        chan2interp{s, o}(t) = e;
                    end
                end
            end
        end
    end
    fprintf(' done.\n')
end
for s = 1:length(subjects)
    interp_manual(s).subject = subjects(s);
    for o = 1:length(session)
        if ~isempty(chan2interp{s, o})
            statement = ['interp_manual(s).' session{o} ' = lwdata.header.chanlocs(chan2interp{s, o}).labels;'];
            eval(statement)
        else
            statement = ['interp_manual(s).' session{o} ' = ''none'';'];
            eval(statement)
        end
    end
end
save(output_file, 'interp_manual', '-append');

% close all figures if required
answer_fig = questdlg('Do you want to close the figures?', 'figures', 'YES', 'NO', 'YES'); 
if strcmp(answer_fig, 'YES')
    figures = findall(0, 'Type', 'figure');
    for f = 1:length(figures)
        if strcmp(figures(f).Name, ' timtopo()')
            close(figures(f))
        elseif strcmp(figures(f).Name, '')
            close(figures(f))
        end
    end
end
clear ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG globalvars LASTCOM PLUGINLIST STUDY figures 

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% preprocess the data
for s = 1:length(subjects) 
    % identify subject
    if subjects(s) < 10
       subj = ['0' num2str(subjects(s))];
    else
       subj = num2str(subjects(s)); 
    end
    fprintf('subject %s:\n', subj)
    for o = 1:length(orientation)
        % determine the session
        if o <= 2
            a = 1;
        else
            a = 2;
        end
        for i = 1:length(intensity)
            fprintf('%s - %s %%rMT\n', strrep(orientation{o}, '_', ' '), intensity{i})
            % determine the name
            if ~isempty(chan2interp{s, a})
                name = sprintf('%s\\%s %s AGSICI P1 S%s %s %s.lw6', folder_output, suffix{1}, prefix, subj, orientation{o}, intensity{i});
            else
                name = sprintf('%s\\%s AGSICI P1 S%s %s %s.lw6', folder_output, prefix, subj, orientation{o}, intensity{i});
            end
    
            % load the data
            option = struct('filename', name);
            lwdata = FLW_load.get_lwdata(option);

            % re-reference to common average
            fprintf('re-referencing ...')
            option = struct('reference_list', {{lwdata.header.chanlocs(1:32).labels}}, 'apply_list', {{lwdata.header.chanlocs(1:32).labels}}, 'suffix', suffix{2}, 'is_save', 0);
            lwdata = FLW_rereference.get_lwdata(lwdata, option);

            % baseline correct
            fprintf('subtracting baseline ...')
            option = struct('operation', 'substract', 'xstart', baseline(1) ,'xend', baseline(2), 'suffix', suffix{3}, 'is_save', 0);
            lwdata = FLW_baseline.get_lwdata(lwdata, option);

            % crop
            fprintf('cropping ...')
            option=struct('xcrop_chk', 1, 'xstart', time_range(1),'xend', time_range(2), 'suffix', suffix{4},'is_save', 1);
            lwdata = FLW_crop_epochs.get_lwdata(lwdata,option);
            fprintf('done.\n')
        end
    end
end

% update the prefix
for s = 2:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end

% extract and save individual data
fprintf('subject ')
for s = 1:length(subjects) 
    % identify subject
    if subjects(s) < 10
       subj = ['0' num2str(subjects(s))];
    else
       subj = num2str(subjects(s)); 
    end
    fprintf(' %s ...', subj)
    for o = 1:length(orientation)
        for i = 1:length(intensity)
            % load the data
            name = sprintf('%s\\%s AGSICI P1 S%s %s %s.mat', folder_output, prefix, subj, orientation{o}, intensity{i});
            load(name)

            % append individual data
            data_muscles_trials(o, i, s, 1:size(data, 1), :, :) = squeeze(data);

            % calculate subject average
            data_muscles(o, i, s, :, :) = squeeze(mean(data, 1));

            % calculate artifact gfp
            data_muscles_gfp(o, i, s, :) = std(squeeze(data_muscles(o, i, s, electrodes, :)), 1);
        end
    end
end
fprintf('done.\n')
save(output_file, 'data_muscles_trials', '-append');
save(output_file, 'data_muscles', '-append');
save(output_file, 'data_muscles_gfp', '-append');

% fill in data for subject 17
data_temp = data_muscles(:, :, [16:18], :, :);
data_muscles(:, :, [16:18], :, :) = [];
data_muscles(:, :, [17:19], :, :) = data_temp;
clear data_temp
counter = 1; 
for a = 1:size(AGSICI_muscle_activity.contraction.data, 1)
    for b = 1:size(AGSICI_muscle_activity.contraction.data, 2)
        data_muscles(counter, :, 16, :, :) = AGSICI_muscle_activity.contraction.data(a, b, :, 16, :, 1:700);   
        data_muscles_gfp(counter, :, 16, :) = AGSICI_muscle_activity.contraction.GFP(a, b, :, 16, 1:700);
        counter = counter + 1;
    end
end
clear a b counter AGSICI_muscle_activity

% identify TOI limits
x_start = (toi(1) + 50)/x_step;
x_end = (toi(2) + 50)/x_step;

% export artifact GFP values in R-compatible table - subject mean values
AGSICI_P1_muscles_mean = table;
row_counter = height(AGSICI_P1_muscles_mean) + 1;
for s = 1:length(subjects) 
    for o = 1:length(orientation)  
        for i = 1:length(intensity)
            % fill mean values
            AGSICI_P1_muscles_mean.subject(row_counter) = s;
            AGSICI_P1_muscles_mean.subject_ID(row_counter) = subjects(s);
            AGSICI_P1_muscles_mean.orientation{row_counter} = orientation(o);
            AGSICI_P1_muscles_mean.intensity(row_counter) = str2num(intensity{i});
            AGSICI_P1_muscles_mean.gfp(row_counter) = mean(data_muscles_gfp(o, :, s, x_start:x_end), [2, 4]);

            % update the counter
            row_counter = row_counter + 1;
        end
    end
end
writetable(AGSICI_P1_muscles_mean, [folder_results '\AGSICI_P1_muscles_mean.csv'])

% export values in R-compatible table - trial values
AGSICI_P1_muscles_trials = table;
row_counter = height(AGSICI_P1_muscles_trials) + 1;
for s = 1:length(subjects) 
    for o = 1:length(orientation)  
        for i = 1:length(intensity)
            % determine number of trials in the category
            for e = 1:size(data_muscles_trials, 4)
                if sum(data_muscles_trials(o, i, s, e, 1, :)) == 0
                    idx(e) = false;
                else
                    idx(e) = true;
                end
            end
            
            for t = 1:sum(idx)
                % calculate GFP 
                gfp = std(squeeze(data_muscles_trials(o, i, s, t, electrodes, :)), 1);
                gfp = mean(gfp(x_start:x_end));

                % fill values
                AGSICI_P1_muscles_trials.subject(row_counter) = s;
                AGSICI_P1_muscles_trials.subject_ID(row_counter) = subjects(s);
                AGSICI_P1_muscles_trials.orientation{row_counter} = orientation(o);
                AGSICI_P1_muscles_trials.intensity(row_counter) = str2num(intensity{i});
                AGSICI_P1_muscles_trials.trial(row_counter) = t;
                AGSICI_P1_muscles_trials.gfp(row_counter) = gfp;

                % update the counter
                row_counter = row_counter + 1;
            end
        end
    end
end
writetable(AGSICI_P1_muscles_trials, [folder_results '\AGSICI_P1_muscles_trials.csv'])

clear session e f o s t 

%% functions
function export_EEGLAB(lwdata, filename, subj)
    % dataset
    EEG.setname = filename;
    EEG.filename = [];
    EEG.filepath = [];
    EEG.subject = subj; 
    EEG.session = 1;
    
    % time properties
    EEG.nbchan = lwdata.header.datasize(2);
    EEG.trials = lwdata.header.datasize(1);
    EEG.pnts = lwdata.header.datasize(6);
    EEG.srate = 1/lwdata.header.xstep;
    EEG.times = lwdata.header.xstart + (0:EEG.pnts-1)*lwdata.header.xstep;
    EEG.xmin = EEG.times(1);
    EEG.xmax = EEG.times(end);
    EEG.data = permute(single(lwdata.data),[2,6,1,3,4,5]);
    EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'SEEG_enabled');
    EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'topo_enabled');
    
    % create events with appropriate latencies
    EEG.event = lwdata.header.events;
    if ~isempty(EEG.event)
        [EEG.event.type] = EEG.event.code;
        for e = 1:length(EEG.event)
            EEG.event(e).latency = (e-1)*EEG.pnts + 2001;
        end
        EEG.event = rmfield(EEG.event,'code');
    end
    
    % create required empty fields
    EEG.icawinv = [];
    EEG.icaweights = [];
    EEG.icasphere = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    EEG.icaweights = [];
    save([filename,'.set'], 'EEG');
end
function [EEG_out] = SSP_SIR_trials(EEG_in, L, art_topographies, filt_ker, M)
%     L = param(c).L_ave;
%     EEG_in = EEG;
%     art_topographies = rmfield(param, {'PC', 'filter', 'L_ave'});
%     filt_ker = param(c).filter;
%     M = [];
%     EEG_out = EEG_in;
%     
%     % prepare individual filters
%     for t = 1:length(art_topographies)
%         statement = sprintf('P%d = eye(size(EEG_in.data,1)) - art_topographies(t).topographies*art_topographies(t).topographies'';', t);
%         eval(statement)
%     end
%     
%     % create final composite filter
%     statement = 'P = P1';
%     if length(art_topographies) > 1
%         for t = 2:length(art_topographies)
%             statement = [statement sprintf(' * P%d', t)];
%         end
%     end
%     statement = [statement ';'];
%     eval(statement);

    EEG_out = EEG_in;
    P = eye(size(EEG_in.data,1)) - art_topographies*art_topographies';
    
    % apply filter
    for i = 1:size(EEG_in.data,3)

        data = EEG_in.data(:,:,i);

        %Suppressing the artifacts:
        data_clean = P*data;

        %Performing SIR for the suppressed data:
        PL = P*L;

        if isempty (M)
            M = rank(data_clean) -  1 ;
        end

        tau_proj = PL*PL';
        [U,S,V] = svd(tau_proj);
        S_inv = zeros(size(S));
        S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
        tau_inv = V*S_inv*U';
        suppr_data_SIR = L*(PL)'*tau_inv*data_clean;

        %Performing SIR for the original data:
        tau_proj = L*L';
        [U,S,V] = svd(tau_proj);
        S_inv = zeros(size(S));
        S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
        tau_inv = V*S_inv*U';
        orig_data_SIR = L*(L)'*tau_inv*data;

        if isempty(filt_ker)
            data_correct = suppr_data_SIR;
        else
            filt_ker_B = repmat(filt_ker,[size(suppr_data_SIR,1),1]);
            data_correct = filt_ker_B.*suppr_data_SIR + orig_data_SIR - filt_ker_B.*orig_data_SIR;
            filt_ker_B = [];
        end

        EEG_out.data(:,:,i) = data_correct;

    end
end









