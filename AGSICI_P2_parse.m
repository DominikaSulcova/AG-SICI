%% PARSE EVENTS ACCORDING TO EXPERIMENTAL CONDITIONS
% Written by Dominika for AG-SICI project (2021)
% 
% 1) Identifies the order of events for current session
%       - loads conditions matrix --> 'SICI_stim_order_XX.mat'
%       ! check that the matlab file is in the working directory
% 2) Renames events according to stimulation sequence
%       - loads header for each dataset
%       - assigns corresponding codes to events, saves back to letswave 
% 3) Splits datasets per event code
%       - each block is split into 7 datasets per 11 stimuli
%       - data and header saved back to letswave
% 4) Merges epochs according to experimental conditions
%       - loads data from letswave and merges blocks of each condition
%       together 
%       --> 'condition P2 XX'
% 
% Important: check whether all the files are disponible in the current directory 

%% session info
clear all
clc

% dataset
block = 1:9;
condition = {'spTMS TS' 'spTMS CS1' 'spTMS CS2' 'spTMS CS3' 'ppTMS CS1' 'ppTMS CS2' 'ppTMS CS3'};
prefix = 'dc ds art-sup ep dc reref';

% session 
prompt = {'Project:' 'Subject:'};
dlgtitle = 'Subject';
dims = [1 35];
definput = {'P2' '00'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% create a prefix
name = [prefix ' ' session_info{1} ' ' session_info{2}];

%% 1) load order of events
load(['SICI_stim_order_' session_info{2} '.mat'])
stim_order = outcome;
clear outcome

%% 2) rename events 
% loop through datasets
for b = 1:length(block)
    load([name ' b' num2str(block(b)) '.lw6'],'-mat'); 
    
    % control for wrong number of events
    if length(header.events) < 77
        % ask for the right sequence
        disp([session_info{1} ' ' session_info{2} ' block ' num2str(block(b))])
        disp(['Found fewer than the expected number of events: ' num2str(length(header.events))])
        prompt = {'Use following TMS stimuli:'};
        dlgtitle = [session_info{1} ' ' session_info{2} ' block ' num2str(block(b))];
        dims = [1 50];
        definput = {'[1:77]'};
        stims2keep = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
        clear prompt dlgtitle dims definput  
        
        % replace the stimuli in the stim_order
        stim_order(1:length(stims2keep), block(b)) = stim_order(stims2keep(1):stims2keep(end), block(b));
        
    elseif length(header.events) > 77
        disp([session_info{1} ' ' session_info{2} ' block ' num2str(block(b))])
        disp(['Found more than the expected number of events: ' num2str(length(header.events))])
        continue
    end
    
    % replace event codes, save header
    for e = 1:length(header.events)
        header.events(e).code = condition{stim_order(e, block(b))};
    end
    save([name ' b' num2str(block(b)) '.lw6'],'header');
end
clear b e c header

%% 3) split per event code
% loop through datasets
for b = 1:length(block)
    % load dataset
    load([name ' b' num2str(block(b)) '.mat']); 
    load([name ' b' num2str(block(b)) '.lw6'],'-mat'); 
    
    % split into datasets per eventcode
    [out_datasets, message_string] = RLW_segmentation2(header, data, condition, 'x_start', -1, 'x_duration', 3);    
    
    % save datasets 
    for c = 1:length(condition)
        % header
        header = out_datasets(c).header;
        header.name = [condition{c} ' ' session_info{1} ' ' session_info{2} ' b' num2str(block(b))];
        save([header.name '.lw6'], 'header');
        
        % data
        data = out_datasets(c).data;
        save([header.name '.mat'], 'data');
    end
end
clear b c data header message_string

%% 4) merge epochs according to conditions
for c = 1:length(condition)
    % merge data
    merge_idx = 1:length(block);
    datasets = struct;
    for b = 1:length(block)
        load([condition{c} ' ' session_info{1} ' ' session_info{2} ' b' num2str(block(b)) '.lw6'], '-mat');
        datasets(b).header = header;
        load([condition{c} ' ' session_info{1} ' ' session_info{2} ' b' num2str(block(b)) '.mat']);
        datasets(b).data = data;
    end
    [header, data, message_string] = RLW_merge_epochs(datasets, merge_idx); 
    
    % save datasets
    header.name = [condition{c} ' ' session_info{1} ' ' session_info{2}];
    save([header.name '.lw6'], 'header')  
    save([header.name '.mat'], 'data') 
 end
clear c b merge_idx data header message_string 

