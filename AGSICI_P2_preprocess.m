%% AG SICI: EEG DATA PRE-PROCESSING
% Written by Dominika for the AG-SICI project, part 2 (2022)
% 
% - each section of the script performs one data pre-processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - when the manual input is needed (ICA, visual inspection), the process
%   is performed in the letswave GUI, then the corresponding section of the
%   script takes care of the logfile update and encodes the output
%   parameters in a structure 'AGSICI_info'
% 
% Output:
%   --> variables are saved in a global MATLAB output file 'AG-SICI_P2.mat'
%   --> Ragu datasets are exported to the 'export' subfolder in results 
% 
% 1) INFO FILE
%       - launch info structure 'AGSICI_info'
%       - encode individual rMT
% 
% 2) PRE-PROCESSING BLOCK 1
%       - load raw TEP data 
%           ! musct be in lw format, import from MEGA using 'AGSICI_import.m'
%       - correct channel location info 
%       - re-reference to common average
%       - segment to [-1 1]s epochs
%       - remove DC shift + linear trend 
%       - remove TMS artifact at [-5 10]ms --> cubic interpolation
%       - downsample to 2kHz
% 
% 3) PARSE EVENTS + MERGE PER CONDITION
%       - load matrix with stimulation order 'SICI_stim_order_subject'
%       - rename event codes accordingly
%       - subset blocks per category (7 x 11 events)
%       - merge subset in one dataset per category
% 
% 4) ICA 1 - COMPUTE MATRIX
%       - calculate squared ICA matrix
% 
% 5) ENCODE ICA 1  
%       - extracts indexes of removed components
%       - encodes to the info structure and to the logfile
% 
% 6) PRE-PROCESSING BLOCK 1
%       - bandpass Butterworth filter at [0.1 80]Hz
%       - FFT notch filter at 50Hz - width 2 Hz, slope 2 Hz
%       - crop at [-0.75 0.75]s
% 
% 8) ENCODE ICA 2  
%       !!! second round of ICA is run manually from letswave, it is also
%       encoded manually to the logfile !!!
%       - extracts the indexes of removed ICs, encodes to the info
%       structure
% 
% 9) PRE-PROCESSING BLOCK 3
%       - baseline correct by subtracting tha average of [-0.2 -0.005]s
%       - average across trials
% 
% 10) GROUP AVERAGE & EXPORT
%       - merge subject data per category, average and save for letswave
%       - write in .csv formate accepted by Ragu
%       - create a .xyz file with montage coordinates

%% parameters
clear all; clc

% ----------- dataset -----------
study = 'P2';
<<<<<<< HEAD:AGSICI_P2_preprocess.m
subject = [1:18];
=======
subject = [15];
>>>>>>> main:AGSICI_P2_process.m
block = 1:9;
protocol = {'spTMS', 'ppTMS'};
intensity = {'CS1', 'CS2', 'CS3'};
electrodes = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2',...
    'CP1','CP2','FC5','FC6','CP5','CP6','P5', 'P6', 'C1','C2'};
% --------------------------------

% define conditions
condition{1} = [protocol{1} ' TS'];
for p = 1:length(protocol)
    for i = 1:length(intensity)
        condition{(p-1)*length(intensity) + i + 1} = [protocol{p} ' ' intensity{i}];
    end
end
clear p i

% choose the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');
addpath(folder_git)

% choose the input folder 
folder_input = uigetdir(pwd, 'Coose the input folder');

% choose the results folders + output
folder_results = uigetdir(pwd, 'Coose the results folder');
output_file = [folder_results '\AG-SICI_' study '.mat' ];

% choose the folder with logfiles
folder_logfiles = uigetdir(pwd, 'Coose the logfile folder');

%% 1) INFO FILE
% load if available
if exist(output_file) 
    load(output_file, 'AGSICI_info');
end

% loop through datasets
for s = 1:length(subject)
    % encode subject number
    AGSICI_info(subject(s)).subject_ID = subject(s);
    
    % identify logfile
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
        
    % extract & encode rMT value
    rMT = get_rMT(filename);
    AGSICI_info(subject(s)).rMT = rMT;    
end
clear s subj filename rMT EMG_ID statement

% save the output file
save(output_file, 'AGSICI_info', '-append');

%% 2) PRE-PROCESSING BLOCK 1
% ----- section input -----
suffix = {'reref' 'ep' 'dc' 'art-sup' 'ds'};
% -------------------------
% load default header
load([folder_git '\header_example.mat'])
        
% loop through subjects
for s = 1:length(subject)
    fprintf('subject %d:\n', subject(s))
    
    % loop through blocks    
    for b = block
        fprintf('block %d:\n', b)
        
        % identify subject
        if subject(s) < 10
           subj = ['0' num2str(subject(s))];
        else
           subj = num2str(subject(s)); 
        end
        
        % get the dataset        
        dataset = sprintf('%s\\%s %s b%d.lw6', folder_input, study, subj, b);
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);
        
        % correct channel locations
        lwdata.header.chanlocs = header.chanlocs(1:length(electrodes));
                
        % re-reference to common average
        fprintf('re-referencing to common average ... ')
        option = struct('reference_list', {electrodes}, 'apply_list', {electrodes},...
            'suffix', suffix{1}, 'is_save', 0);
        lwdata = FLW_rereference.get_lwdata(lwdata, option);
                
        % segment
        fprintf('epoching from -1 to 1 s relative to the event ... ')
        option = struct('event_labels', {{'Stimulation'}}, 'x_start', -1, 'x_end', 1, 'x_duration', 2, ...
            'suffix', suffix{2}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);

        % remove DC + linear detrend
        fprintf('removing DC and applying linear detrend...')
        option = struct('linear_detrend', 1, 'suffix', suffix{3}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata,option);
                
        % artifact removal, interpolation
        fprintf('interpolating the TMS artifact from -5 to 10 ms...')
        [lwdata.header, lwdata.data, message_string] = RLW_suppress_artifact_event(lwdata.header, lwdata.data,...
            'xstart', -0.005, 'xend', 0.01, 'event_code', 'Stimulation', 'interp_method', 'pchip');
                                
        % downsample
        fprintf('downsampling...\n')
        option = struct('x_dsratio', 10, 'suffix', [suffix{5} ' ' suffix{4}], 'is_save', 1);
        lwdata = FLW_downsample.get_lwdata(lwdata, option);
    end     

    % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    logfile_entry('block 1', filename)
end
clear header s b subj dataset filename lwdata message_string
disp('Finished! Sendwich?')

% update prefix
prefix = study;
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix

%% 3) PARSE EVENTS + MERGE PER CONDITION
% loop through subjects
for s = 1:length(subject)
    fprintf('subject %d:\n', subject(s))
    
    % identify subject
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end
    
    % load order of events
    load([folder_input '\SICI_stim_order_' subj '.mat'])
    stim_order = outcome;
    clear outcome
    
    % rename events
    fprintf('assigning event labels...')
    for b = block       
        % load dataset
        load([prefix ' ' subj ' b' num2str(b) '.lw6'], '-mat'); 
    
        % control for wrong number of events
        if length(header.events) < 77
            % ask for the right sequence
            fprintf('\nFound fewer than the expected number of events: %d\n', length(header.events))
            prompt = {'Use following TMS stimuli:'};
            dlgtitle = ['Subject ' subj ' - block ' num2str(b)];
            dims = [1 50];
            definput = {'[1:77]'};
            stims2keep = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
            clear prompt dlgtitle dims definput  

            % replace the stimuli in the stim_order
            stim_order(1:length(stims2keep), b) = stim_order(stims2keep, b);

        elseif length(header.events) > 77
            fprintf('\nFound more than the expected number of events: %d\n', length(header.events))
            continue
        end
    
        % replace event codes, save header
        for e = 1:length(header.events)
            header.events(e).code = condition{stim_order(e, block(b))};
        end
        save([prefix ' ' subj ' b' num2str(b) '.lw6'], 'header');
    end
    clear b stims2keep e header

    % split per event code
    fprintf('parsing per condition...')
    for b = block
        % get the dataset        
        dataset = sprintf('%s %s b%d.lw6', prefix, subj, b);
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);

        % split into datasets per eventcode
        for c = 1:length(condition)
            % segment
            option = struct('event_labels', {condition(c)}, 'x_start', -1, 'x_end', 1, 'x_duration', 2, ...
                'suffix', '', 'is_save', 0);
            lwdataset = FLW_segmentation_separate.get_lwdataset(lwdata, option); 
            
            % rename 
            lwdataset.header.name = sprintf('%s %s %s', study, subj, condition{c});

            % save into a dataset for merging
            data2merge(c, b) = lwdataset;
        end
    end
    clear b dataset option lwdata lwdataset 

    % merge epochs according to conditions
    fprintf('merging datasets...')
    for c = 1:length(condition)
        % subset the dataset
        lwdataset = data2merge(c, :)';
        
        % merge subjects & save 
        option = struct('type', 'epoch', 'suffix', '', 'is_save', 1);
        lwdata = FLW_merge.get_lwdata(lwdataset, option); 
    end
    clear c 

     % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    logfile_entry('parse', filename)   
    fprintf('done.\n')
end
clear s subj lwdataset stim_order lwdata filename option data2merge

% update prefix
prefix = study;

%% 4) ICA 1 - COMPUTE MATRIX
% ----- section input -----
% prefix = study;
suffix = 'prea';
components = 31;
% -------------------------
% define ICA mode
if components == length(electrodes)
    option_ica = struct('ICA_mode', 1, 'algorithm', 1, 'suffix', suffix, 'is_save', 1);
else
    option_ica = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', components, 'suffix', suffix, 'is_save', 1);
end

% loop through datasets
for s = 1:length(subject) 
    fprintf('subject %d:\n', subject(s))
    
    % identify subject
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end

    % load the dataset
    for c = 1:length(condition)
        dataset{c} = sprintf('%s %s %s.lw6', prefix, subj, condition{c});
    end
    option = struct('filename', {dataset});
    lwdataset = FLW_load.get_lwdataset(option);

    % compute the ICA matrix
    fprintf('computing ICA matrix...')
    option = option_ica;
    lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
    fprintf('done.\n')
end
clear components option_ica s subj c subject dataset option lwdataset

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

%% 5) ENCODE ICA 1  
% ----- section input -----
% prefix = ['prea ' study];
suffix = 'sp_filter';
% -------------------------
% load info file
load(output_file, 'AGSICI_info');

% loop through datasets
for s = 1:length(subject)
    fprintf('subject %d:\n', subject(s))
    
    % identify subject
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end
           
    % load one of the datasets
    dataset = sprintf('%s %s %s spTMS TS.lw6', suffix, prefix, subj);
    option = struct('filename', dataset);
    lwdata = FLW_load.get_lwdata(option);
            
    % extract numbers of removed components
    components = lwdata.header.history(length(lwdata.header.history)).option.remove_idx;  
    disp(['Components removed: ' num2str(components)])
        
    % encode to the output structure        
    AGSICI_info(subject(s)).ICA_1 = components;
        
    % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    logfile_entry('ICA 1', filename, 'components', components)    
end
clear s filename dataset option lwdata subj components
disp('Finished! Ad Hoc?!')

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

% save the output file
save(output_file, 'AGSICI_info', '-append');

%% 6) PRE-PROCESSING BLOCK 2
% ----- section input -----
% prefix = ['sp_filter prea ' study];
suffix = {'bandpass' 'notch' 'crop'};
bandpass = [0.1 80];
crop_window = [-0.75 0.75];
% -------------------------
% loop through datasets
for s = 1:length(subject)
    fprintf('subject %d:\n', subject(s))
    
    for c = 1:length(condition)
        fprintf('%s:\n', condition{c})
        
        % identify subject
        if subject(s) < 10
           subj = ['0' num2str(subject(s))];
        else
           subj = num2str(subject(s)); 
        end
           
        % get the dataset
        dataset = sprintf('%s %s %s.lw6', prefix, subj, condition{c});
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);

        % bandpass
        fprintf('Applying Butterworth bandpass filter...')
        option = struct('filter_type', 'bandpass', 'high_cutoff', bandpass(2),'low_cutoff', bandpass(1),...
            'filter_order', 4, 'suffix', suffix{1}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);

        % notch
        fprintf('Applying FFT notch filter...')
        option = struct('filter_type', 'notch', 'notch_fre', 50, 'notch_width', 2, 'slope_width', 2,...
            'harmonic_num', 2,'suffix', suffix{2},'is_save', 0);
        lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);

        % crop
        fprintf('Cropping the data...')
        option = struct('xcrop_chk', 1, 'xstart', crop_window(1), 'xend', crop_window(2), 'suffix', suffix{3}, 'is_save', 1);
        lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
        fprintf('Done.\n')
    end
    % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    logfile_entry('block 2', filename)
end
clear s c subj dataset option lwdata bandpass crop_window
disp('Finished! Sendwich?')

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix 

%% 7) VISUAL INSPECTION 
% ----- section input -----
% prefix = ['crop notch bandpass sp_filter prea ' study];
suffix = 'visual';
% -------------------------
% load info file
load(output_file, 'AGSICI_info');

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

% loop through datasets
for s = 2:length(subject)
    % identify subject
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end
    
    % loop through conditions
    for c = 1:length(condition)       
        % load the dataset
        dataset = sprintf('%s %s %s.lw6', prefix, subj, condition{c}); 
        load(dataset, '-mat')
            
        % extract removed epochs
        epochs = header.history.configuration.parameters.rejected_epochs;  
            
        % extract number of retained epochs
        epochs_n = header.datasize(1);  
        
        % encode to the output structure    
        AGSICI_info(s).epochs(c) = epochs_n;
        AGSICI_info(s).epochs_removed(c, 1:length(epochs)) = epochs;
    end
    
    % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    epochs_n = AGSICI_info(s).epochs;
    epochs = AGSICI_info(s).epochs_removed; 
    logfile_entry('visual', filename, 'epochs_n', epochs_n, 'epochs', epochs)   
end
clear s c subj dataset header epochs epochs_n filename

% save the output file
save(output_file, 'AGSICI_info', '-append');

%% 8) ENCODE ICA 2  
% ----- section input -----
% prefix = ['visual crop notch bandpass sp_filter prea ' study];
suffix = 'icfilt ica';
% -------------------------
% load info file
load(output_file, 'AGSICI_info');

% loop through datasets
for s = 1:length(subject)   
    % identify subject
    if subject(s) < 10
       subj = ['0' num2str(subject(s))];
    else
       subj = num2str(subject(s)); 
    end
           
    % load one of the datasets
    dataset = sprintf('%s %s %s spTMS TS.lw6', suffix, prefix, subj);
    load(dataset, '-mat')
            
    % extract numbers of removed components
    IC_n = size(header.history(length(header.history)).configuration.parameters.ICA_um, 1);
    IC = 1:IC_n;
    IC = IC(~ismember(IC, header.history(length(header.history)).configuration.parameters.IC_list));  
        
    % encode to the output structure        
    AGSICI_info(subject(s)).ICA_2 = IC; 
end
clear s subj dataset header IC_n IC

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

% save the output file
save(output_file, 'AGSICI_info', '-append');

%% 9) PRE-PROCESSING BLOCK 3
% ----- section input -----
% prefix = ['icfilt ica visual crop notch bandpass sp_filter prea ' study];
suffix = {'bl_correct' 'avg'};
baseline = [-0.2 -0.005];
% -------------------------
% loop through datasets
for s = 1:length(subject)  
    fprintf('subject %d ... ', subject(s))
    for c = 1:length(condition)        
        % identify subject
        if subject(s) < 10
           subj = ['0' num2str(subject(s))];
        else
           subj = num2str(subject(s)); 
        end
           
        % get the dataset
        dataset = sprintf('%s %s %s.lw6', prefix, subj, condition{c});
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);

        % baseline correct
        option = struct('operation', 'substract', 'xstart', baseline(1), 'xend', baseline(2), 'suffix', suffix{1}, 'is_save', 0);
        lwdata = FLW_baseline.get_lwdata(lwdata, option);

        % average
        option = struct('operation', 'average', 'suffix', suffix{2}, 'is_save', 1);
        lwdata = FLW_average_epochs.get_lwdata(lwdata, option);
    end
    % update the logfile
    filename = [folder_logfiles '\AGSICI_' study '_' subj '.txt'];
    logfile_entry('block 3', filename)
end
clear s c subj dataset option lwdata baseline
disp('Finished! Sendwich?')

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix 

%% 10) GROUP AVERAGE & EXPORT
% ----- section input -----
time_window = [-0.05, 0.3];
% -------------------------
% average across subjects by condition
for c = 1:length(condition)
    % load the dataset
    for s = 1:length(subject)
        % identify subject
        if subject(s) < 10
           subj = ['0' num2str(subject(s))];
        else
           subj = num2str(subject(s)); 
        end
            
        % add to the dataset
        dataset{s} = [folder_input '\' prefix ' ' subj ' ' condition{c} '.lw6'];
    end
    option = struct('filename', {dataset});
    lwdataset = FLW_load.get_lwdataset(option);
        
    % merge epochs 
    option = struct('type', 'epoch', 'suffix', '', 'is_save', 0);
    lwdata = FLW_merge.get_lwdata(lwdataset, option);
            
    % rename, save
    lwdata.header.name = ['merged ' condition{c}];
    CLW_save(lwdata)
    
    % average, save
    option = struct('operation', 'average', 'suffix', 'avg', 'is_save', 1);
    lwdata = FLW_average_epochs.get_lwdata(lwdata, option);
end
clear c s subj dataset option lwdata lwdataset

% calculate crop limits
load([folder_input '\' prefix ' 01 spTMS TS.lw6'], '-mat')                    
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;

% write text files for Ragu 
for c = 1:length(condition)
    for s = 1:length(subject)
        % identify subject
        if subject(s) < 10
           subj = ['0' num2str(subject(s))];
        else
           subj = num2str(subject(s)); 
        end

        % load the data
        name_old = [folder_input '\' prefix ' ' subj ' ' condition{c} '.mat'];
        load(name_old)
        data = squeeze(data(:, 1:32, :, :, :, x_start:x_end))';

        % save as .csv               
        name = ['AGSICI ' study ' ' condition{c} ' ' subj '.csv']; 
        writematrix(data, [folder_results '\AG-SICI_' study '_export\' name])
    end
end
clear c s subj name_old data x_start x_end name time_window    

% create the montage file
name = [folder_results '\AG-SICI_' study '_export\AGSICI_montage.xyz'];
fileID = fopen(name, 'a');
fprintf(fileID, '32\r\n');
for a = 1:32
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear name fileID a header   

%% functions
function rMT = get_rMT(filename)
    fileID = fopen(filename, 'r');
    text = textscan(fileID, '%s%s%d', 'Headerlines', 27);
    rMT = text{3}(1);
    fclose(fileID);
end
function logfile_entry(entry, filename, varargin)
    switch entry
        case 'block 1'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, 'DATA CLEANING\r\n');
            fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '1	dataset import: used a matlab script ''AGSICI_import.m''\r\n');
            fprintf(fileID, '		- data automatically loaded from MEGA output\r\n');
            fprintf(fileID, '		- ''Out'' category discarded, ''Stimulation'' category kept\r\n');
            fprintf(fileID, '		- checked for missing/extra events, first blind event discarded\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '2	first step preprocessing: used a matlab script ''AGSICI_P2_process.m''\r\n');
            fprintf(fileID, '		- electrodes were assigned standardized 10/20 labels\r\n');
            fprintf(fileID, '		- rereferenced to the common average\r\n');
            fprintf(fileID, '		- segmentation relative to events --> epochs [-1 1]s, epoch size 60000 bins\r\n');
            fprintf(fileID, '		- DC removal + linear detrend\r\n');
            fprintf(fileID, '		- TMS artifact suppression --> [-0.005 0.01]s\r\n');
            fprintf(fileID, '		- downsample by 1/10\r\n');
            fprintf(fileID, '	--> file prefix: ds art-sup dc ep reref P2\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);

        case 'parse'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '3  	events parsed and merged in one dataset per condition\r\n');
            fprintf(fileID, '	--> new name: P2 %s spTMS TS/spTMS CS1/spTMS CS2/spTMS CS3/ppTMS CS1/ppTMS CS2/ppTMS CS3\r\n', filename(end-5 : end-4)); 
            fprintf(fileID, '\r\n');
            fclose(fileID);

        case 'ICA 1'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'components'));
                if ~isempty(a)
                    components = varargin{a + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '4	preliminary ICA was computed to get rid of the major decay artifact\r\n');
            fprintf(fileID, '		- squared matrix --> 32 ICs\r\n');
            fprintf(fileID, '		- removed component(s): IC %s\r\n', num2str(components));
            fprintf(fileID, '	--> file prefix: sp_filter prea\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 2'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '5	frequency filters applied\r\n');
            fprintf(fileID, '		- FFT notch filter --> 50 Hz, width 2 + slope 2\r\n');
            fprintf(fileID, '		- Butterworth bandpass filter --> [0.1 80]Hz, 4th order\r\n');
            fprintf(fileID, '	--> file prefix: bandpass notch\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '6 	signal cropped to get rid of edge artifacts\r\n');
            fprintf(fileID, '		- [-0.75 0.75]s\r\n');
            fprintf(fileID, '	--> file prefix: crop\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'visual'
            if ~isempty(varargin)
                e = find(strcmpi(varargin, 'epochs'));
                if ~isempty(e)
                    epochs = varargin{e + 1};
                end
                
                n = find(strcmpi(varargin, 'epochs_n'));
                if ~isempty(n)
                    epochs_n = varargin{n + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '7	noisy epochs were removed based on visual inspection\r\n');
            fprintf(fileID, '		- spTMS TS: %s --> %d epochs retained\r\n', num2str(epochs(1, epochs(1, :) > 0)), epochs_n(1));
            fprintf(fileID, '		- spTMS CS1: %s --> %d epochs retained\r\n', num2str(epochs(2, epochs(2, :) > 0)), epochs_n(2));
            fprintf(fileID, '		- spTMS CS2: %s --> %d epochs retained\r\n', num2str(epochs(3, epochs(3, :) > 0)), epochs_n(3));
            fprintf(fileID, '		- spTMS CS3: %s --> %d epochs retained\r\n', num2str(epochs(4, epochs(4, :) > 0)), epochs_n(4));
            fprintf(fileID, '		- ppTMS CS1: %s --> %d epochs retained\r\n', num2str(epochs(5, epochs(5, :) > 0)), epochs_n(5));
            fprintf(fileID, '		- ppTMS CS2: %s --> %d epochs retained\r\n', num2str(epochs(6, epochs(6, :) > 0)), epochs_n(6));
            fprintf(fileID, '		- ppTMS CS3: %s --> %d epochs retained\r\n', num2str(epochs(7, epochs(7, :) > 0)), epochs_n(7));
            fprintf(fileID, '	--> file prefix: visual\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 3'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '9  baseline corrected - subtracted mean of [-0.2 -0.005]s\r\n');
            fprintf(fileID, '   --> file prefix: bl_correct\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '10 signal averaged across epochs\r\n');
            fprintf(fileID, '   --> file prefix: avg\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
    end
end