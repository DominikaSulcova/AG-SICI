%  processing with SSP-SIR

%% parameters
clear all; clc;

% dataset
subject = 1;
sequence = {'along' 'across'; 'reversed' 'normal'};
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
folder_toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave + eeglab masterfiles
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
suffix = {'reref' 'dc' 'ep' 'art-sup' 'ds'};
eventcode = 'B - Stimulation';
epoch = [-0.2 0.5];
interp = [-0.005, 0.01];
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
        
        % rereference to common average
        fprintf('...rereferencing to average')
        [header, data, ~] = RLW_rereference(header, data, 'apply_list', {header.chanlocs.labels}, 'reference_list', {header.chanlocs.labels});

        % remove DC and apply linear detrend
        fprintf('...removing DC + detrending')
        [header, data, ~] = RLW_dc_removal(header, data, 'linear_detrend', 1);
        
        % segment 
        fprintf('...epoching')
        [header, data, ~] = RLW_segmentation(header, data, {{eventcode}}, 'x_start', -1, 'x_duration', 3);
        
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
        header.name = [suffix{5} ' ' suffix{4} ' ' suffix{3} ' ' suffix{2} ' ' suffix{1} ' ' header.name];
        CLW_save([],header, data);
    end
end

% create prefix
prefix = [suffix{5} ' ' suffix{4} ' ' suffix{3} ' ' suffix{2} ' ' suffix{1}];
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
    fprintf('Session %d - %s:', s, sequence{1, s})
    
    % identify stimulation order
    switch sequence{2, s}
        case 'normal'
            index = repmat(['NR'], 1, length(block)/2);
        case 'reversed'
            index = repmat(['RN'], 1, length(block)/2);
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
clear o files f filename option lwdata lwdataset data2merge i 

%% visual inspection, bad channels removal
% ----- section input -----
subject = 21;
session = 'across';
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

% identify subject
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% cycle through current dirrections
fprintf('Subject %s - %s:\n', subj, session)
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
% ----- section input -----
subject = 1;
% -------------------------
% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% identify subject
if subject < 10
   subj = ['0' num2str(subject)];
else
   subj = num2str(subject); 
end

% export in .set format
fprintf('exporting to EEGLAB:\n')
for o = 1:length(orientation)
    fprintf('%s - ', orientation{o})
    for i = 1:length(intensity)
        fprintf('%s ... ', intensity{i})
        % identify the file
        file = dir([prefix ' AGSICI P1 S' subj ' ' orientation{o} ' ' intensity{i} '.lw6']);
        
        % export
        FLW_export_EEGLAB.get_lwdata('filename', file.name, 'pathname', file.folder);
    end
    fprintf('\n')
end
fprintf('done.\n')
clear subj o i file

%% launch EEG lab
% add eeglab to the of search path
addpath(fullfile(folder_toolbox,'eeglab2022.1'));

% add fastica to the of search path
addpath(fullfile(folder_toolbox,'FastICA_25'));

% launch eeglab and generate an empty EEGLAB structure
eeglab

%% process the data: ICA + SOUND + SSP-SIR
% ----- section input -----
subj = '01';
ori = 'along_normal';
int = '100';
baseline = [-0.25 -0.005];
ICA_comp = 29;
lambda = 0.5;
% ------------------------- 
% load the dataset
filename = sprintf('visual AGSICI P1 S%s %s %s.set', subj, ori, int);
EEG = pop_loadset('filename', filename, 'filepath', folder_output);
eeglab redraw   
        
%% baseline correct - subtraction
EEG = pop_rmbase(EEG, baseline,[]);

% visual check
figure; 
pop_timtopo(EEG, [-100  300], [25  45  75  100 180], 'Original data');
        
%% bad channels
% save the original channels locations 
pop_saveset(pop_select(EEG, 'trial', 1), 'filename', [EEG.setname ' all_channels.set'], 'filepath', folder_output);

% visualize the average response to identify possible bad channels
figure; 
pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);

% remove bad channels 
pop_eegplot(pop_select(EEG, 'time', [-0.1 0.3]), 1, 1, 1);
EEG = pop_select(EEG);
        
%% remove bad trials
% % verify if the final number of trials is acceptable --> SNR should be less than 10%
% n_original = 100;
% n_accepted = 100 - 15;
% SNR = (sqrt(n_accepted)-sqrt(n_original))/sqrt(n_original);

% % remove bad trials
% pop_eegplot(EEG, 1, 1, 1);
% pop_saveset(EEG, 'filename', [EEG.setname ' trial_select.set'], 'filepath', folder_output);

%% ICA --> remove ocular artifacts and TMS-independent muscular activity
% determine number of components
if EEG.nbchan < 32
    n_comp = ICA_comp - (32 - EEG.nbchan);
else
    n_comp = ICA_comp;
end
EEG = pop_tesa_pcacompress(EEG, 'compVal', n_comp, 'plot', 'on');

% run ICA 
EEG = pop_tesa_fastica(EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off');
EEG = pop_tesa_compplot(EEG,'figSize', 'large', 'plotTimeX', [-0.5 0.5], 'plotFreqX', [1 100],...
    'freqScale', 'log', 'saveWeights','off');

% baseline correct & save
EEG = pop_rmbase(EEG, baseline, []);
pop_saveset(EEG, 'filename', EEG.setname, 'filepath', folder_output);

% visualize
figure; 
pop_timtopo(EEG, [-100  300], [25  45  75  100 180], 'Data after ICA');

%% SSP-SIR --> leftover muscular artifact
% - use spherical model 
EEG = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', [0,12], 'PC', []);

clear baseline ICA_comp n_comp o i filename










