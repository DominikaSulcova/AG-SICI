%% PARSE EVENTS ACCORDING TO EXPERIMENTAL CONDITIONS
% Written by Dominika for AG-SICI project (2021)
% 
% 1) Identifies the order of events for current session
%       - loads matrices for both types of TMS scripts:    
%         normal or reversed current direction --> 'stim_order_00_N/R.mat'
%       - creates true sequence of stiumuli based on session info
% 2) Renames events according to stimulation sequence
%       - loads header for each dataset
%       - assigns corresponding codes to events, saves back to letswave 
% 3) Merges epochs according to experimental conditions
%       - loads 
% 
% Important: check whether all the files are disponible in the current directory 

%% session info
clear all
clc

% dataset
block = 1:8;
intensity = {'stim_100' 'stim_120' 'stim_140'};
prefix = 'dc ds art-sup ep dc reref';

% session 
prompt = {'Project:' 'Subject:' 'Session:'};
dlgtitle = 'Subject';
dims = [1 35];
definput = {'P1' '00' 'S1'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% coil placement
session_info{4} = questdlg('Coil placement?', 'Experimental condition',...
    'along', 'across', 'along');

% create a prefix
name = [prefix ' ' session_info{1} ' ' session_info{2} ' ' session_info{4}];

% choose real current direction
session_info{5} = questdlg('Real current direction?', 'Experimental condition',...
    '1 - normal  2 - reversed', '1 - reversed  2 - normal', '1 - normal  2 - reversed');

% choose real current direction
index = questdlg('TMS protocol sequence?', 'Parsing sequence',...
    '1 - normal  2 - reversed', '1 - reversed  2 - normal', '1 - normal  2 - reversed');

%% 1) load the order of events
% identify sequence of TMS protocols
switch index
    case '1 - normal  2 - reversed'
        index = repmat(['NR'], 1, numel(block)/2);
    case '1 - reversed  2 - normal'
        index = repmat(['RN'], 1, numel(block)/2);
end

% load stim order matrices
load(['stim_order_' session_info{2} '_N.mat'])
stim_order_N = outcome;
load(['stim_order_' session_info{2} '_R.mat'])
stim_order_R = outcome;
clear outcome

% create final stim sequence
for s = 1:length(block)
    statement = ['stim_order(:, s) = stim_order_' index(s) '(:, s);'];
    eval(statement)
end
clear statement stim_order_N stim_order_R index s

%% 2) rename events
% loop through datasets
for b = 1:length(block)
    load([name ' b' num2str(block(b)) '.lw6'],'-mat'); 
    for e = 1:length(header.events)
        header.events(e).code = num2str(stim_order(e, b));
    end
    save([name ' b' num2str(block(b)) '.lw6'],'header');
end
clear b e header

%% 3) merge epochs according to conditions
for i = 1:length(intensity)
    for b = 1:length(block)
    end
end



