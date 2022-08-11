%% AG SICI: PROCESS EEG DATA
% Written by Dominika for the AG-SICI project, part 2 (2022)
% 
% - each section of the script performs one data processing step
% - the information about each finished process is automatically encoded in
%   a .txt logfile and saved to the shared folder
% - when the manual input is needed (ICA, visual inspection), the process
%   is performed in the letswave 7 GUI and the corresponding section of the
%   script takes care of the logfile update and encodes the output
%   parameters in a structure 'AGSICI_info'
% 
% Output:
%   --> figures are saved in a folder 'AG-SICI_P2_figures'
%   --> variables are saved in a global MATLAB output file 'AG-SICI_P2.mat'
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
% 4) PRELIMINARY ICA MATRIX
%       - calculate preliminary ICA matrix with chosen number of
%       components


%% parameters
clear all; clc

% ----------- dataset -----------
study = 'P2';
subject = [1:5];
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

% create the results folders
folder_results = uigetdir(pwd, 'Coose the results folder');
folder_figures = [folder_results '\AG-SICI_' study '_figures'];
output_file = [folder_results '\AG-SICI_' study '.mat' ];

% choose the folder with logfiles
folder_logfiles = uigetdir(pwd, 'Coose the logfile folder');

%% 1) INFO FILE
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
if ~exist(output_file) 
    save(output_file, 'AGSICI_info');
else
    save(output_file, 'AGSICI_info', '-append');
end

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
clear header s b subj dataset filename lwdata
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

%% 4) PRELIMINARY ICA MATRIX
% ----- section input -----
suffix = 'prea';
components = 32;
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

%% -------------------------------- FROM HERE NOT YET MODIFIED -------------------------------
% preliminary ICA 
% ----- section input -----
suffix = 'prefilt';
% -------------------------
% load output file
load(output_file);

% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        disp([subject ' ' condition{c}])
        
        % load one of the datasets
        dataset = [folder_input '\sp_filter ' prefix ' ' subject ' ' condition{c} ' ' time{1} '.lw6'];
        option = struct('filename', dataset);
        lwdata = FLW_load.get_lwdata(option);
            
        % extract numbers of removed components
        components = lwdata.header.history(length(lwdata.header.history)).option.remove_idx;  
        disp(['Components removed: ' num2str(components)])
        
        % encode to the output structure        
        statement = ['CAPSTEP_info(participant(p)).preICA.' condition{c} ' = components;'];
        eval(statement)
        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('preICA', filename, 'components', components)    
    end
end
clear p c subject filename dataset option lwdata 
disp('Finished! Ad Hoc?!')

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

% save the output file
save(output_file, 'CAPSTEP_info', '-append');

%% 7) processing block 2
% ----- section input -----
suffix = {'bandpass' 'notch' 'crop'};
bandpass = [0.1 80];
crop_window = [-0.75 0.75];
% participant = [];
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        
        for t = 1:length(time)     
            % get the dataset
            disp([subject ' ' condition{c} ' ' time{t}])
            dataset = [folder_input '\' prefix ' ' subject ' ' condition{c} ' ' time{t} '.lw6'];
            option = struct('filename', dataset);
            lwdata = FLW_load.get_lwdata(option);

            % bandpass
            disp('Applying Butterworth bandpass filter...')
            option = struct('filter_type', 'bandpass', 'high_cutoff', bandpass(2),'low_cutoff', bandpass(1),...
                'filter_order', 4, 'suffix', suffix{1}, 'is_save', 0);
            lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);

            % notch
            disp('Applying FFT notch filter...')
            option = struct('filter_type', 'notch', 'notch_fre', 50, 'notch_width', 2, 'slope_width', 2,...
                'harmonic_num', 2,'suffix', suffix{2},'is_save', 0);
            lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);

            % crop
            disp('Cropping the data [0.75 0.75]s...')
            option = struct('xcrop_chk', 1, 'xstart', crop_window(1), 'xend', crop_window(2), 'suffix', suffix{3}, 'is_save', 1);
            lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);

        end        
        % update the logfile
        filename = [folder_logfiles '\CAPS-TEP_' subject '_' condition{c} '.txt'];
        logfile_entry('block 2', filename)
    end
end
clear p c subject t dataset option lwdata bandpass crop_window
disp('Finished! Sendwich?')

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix 

%% 8) visual inspection 
% ----- section input -----
suffix = 'visual';
% participant = [];
% -------------------------
% load output file
load(output_file);

% update prefix
prefix = [suffix ' ' prefix];
clear suffix

% loop through datasets
for p = 1:length(participant)
    for c = 1:length(condition)
        % determine the subject
        if participant(p) < 10
           subject = ['0' num2str(participant(p))];
        else
           subject = num2str(participant(p)); 
        end
        disp([subject ' ' condition{c}])
        
        for t = 1:length(time)
            % load the dataset
            dataset = [folder_input '\' prefix ' ' subject ' ' condition{c} ' ' time{t} '.lw6'];
            option = struct('filename', dataset);
            lwdata = FLW_load.get_lwdata(option);
            
            % extract numbers of removed epochs
            epochs = lwdata.header.history(length(lwdata.header.history)).configuration.parameters.rejected_epochs;  
            disp([num2str(length(epochs)) ' epochs removed'])
            
            % extract final number of retained epochs
            n_epochs = lwdata.header.datasize(1);  
            disp([num2str(n_epochs) ' epochs retained'])
        
            % encode to the output structure        
            statement = ['CAPSTEP_info(participant(p)).epochs_removed.' condition{c} ' = epochs;'];
            eval(statement)
            statement = ['CAPSTEP_info(participant(p)).epochs_n.' condition{c} ' = n_epochs;'];
            eval(statement)   
        end
    end
end
clear p c t subject filename dataset option lwdata epochs n_epochs
disp('Finished! Beer?!')

% save the output file
save(output_file, 'CAPSTEP_info', '-append');

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

        case 'bad events'
            a = find(strcmpi(varargin, 'removed'));
            if ~isempty(a)
                removed = varargin{a + 1};
            end
            
            fileID = fopen(filename, 'a');
            if size(removed) == 0
                fprintf(fileID, '4	no missed stimuli --> no faulty events discarded\r\n');
                fprintf(fileID, '\r\n');
                fclose(fileID);
            else
                fprintf(fileID, '4	faulty (missed) events were removed:\r\n');
                for r = 1:size(removed, 1)
                    fprintf(fileID, '		- n. %s in %s\r\n', num2str(removed{r, 2}), time{removed{r, 1}});
                end
                fprintf(fileID, '\r\n');
                fclose(fileID);
            end

        case 'preICA'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'components'));
                if ~isempty(a)
                    components = varargin{a + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, '5	preliminary ICA was computed to get rid of the major decay artifact\r\n');
            fprintf(fileID, '		- squared matrix --> 30 ICs\r\n');
            fprintf(fileID, '		- removed component(s): IC %s\r\n', num2str(components));
            fprintf(fileID, '	--> file prefix: prefilt prea\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 2'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '6	frequency filters applied\r\n');
            fprintf(fileID, '		- FFT notch filter --> 50 Hz, width 2 + slope 2\r\n');
            fprintf(fileID, '		- Butterworth bandpass filter --> [0.1 80]Hz, 4th order\r\n');
            fprintf(fileID, '	--> file prefix: bandpass notch\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '7 	signal cropped to get rid of edge artifacts\r\n');
            fprintf(fileID, '		- [-0.75 0.75]s\r\n');
            fprintf(fileID, '	--> file prefix: crop\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'visual'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '	--> file prefix: visual\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'ICA'
            if ~isempty(varargin)
                a = find(strcmpi(varargin, 'components'));
                if ~isempty(a)
                    components = varargin{a + 1};
                end
            end
            
            fileID = fopen(filename, 'a');
            fprintf(fileID, ['9	second ICA was computed to get rid of all remaining artifacts (' ' independent components)\r\n']);
            fprintf(fileID, '		- squared matrix --> 30 ICs\r\n');
            fprintf(fileID, '		- removed component(s): IC %s\r\n', num2str(components));
            fprintf(fileID, '	--> file prefix: prefilt prea\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
            
        case 'block 3'
            fileID = fopen(filename, 'a');
            fprintf(fileID, '10 baseline corrected - subtracted mean of [-0.2 -0.005]s\r\n');
            fprintf(fileID, '   --> file prefix: bl\r\n');
            fprintf(fileID, '\r\n');
            fprintf(fileID, '11 signal was averaged across epochs\r\n');
            fprintf(fileID, '	--> file prefix: avg\r\n');
            fprintf(fileID, '\r\n');
            fclose(fileID);
    end
end