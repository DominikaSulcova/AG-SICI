%% EEG DATA IMPORT
% Written by Dominika for AG-SICI project (2021)
% 
% 1) Loads raw EEG datasets recorded by the MEGA system
%    and saves them in the current directory as .mat data + .lw6 header
%    ATTENTION: requires a letswave history file EEG_history_import.mat 
% 2) Discards 'Out' events, keeps 'Stimulation' and counts them
% 
%% session info
clear all; clc;

% dataset
block = 1:9;
prompt = {'Project:' 'Subject:'};
dlgtitle = 'Subject';
dims = [1 35];
definput = {'P2' '00'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% create a prefix
prefix = [session_info{1} ' ' session_info{2}];

% choose the folder with raw data
input_folder = uigetdir(pwd, 'Choose folder with raw data');

% specify EEG block subfolders in the imput folder
prompt = {'Admit following EEG blocks:'};
dlgtitle = 'EEG blocks';
dims = [1 35];
definput = {'[1:9]'};
answer = cell2mat(inputdlg(prompt,dlgtitle,dims,definput));
eval(['folder = ' answer ';']); 
clear answer prompt dlgtitle dims definput            

% choose the Git folder 
git_folder = uigetdir(pwd, 'Choose Git folder');

%% import MEGA datasets
% import the datasets - blocks indicated by folders vector
for b = 1:length(block)
    % 1) import the appropriate dataset
    [header, data] = EEG_import_MEGA(input_folder, folder(b));
                
    % create the name for the dataset
    dataset_name = [prefix ' b' num2str(block(b))];
    
    % create the first letswave history entry
    load([git_folder '\EEG_history_import.mat'])
    EEG_history_import.configuration.parameters.input_folder = input_folder;
    EEG_history_import.configuration.parameters.session_number = folder(b);   
    header.history(1) =  EEG_history_import;
    
    % 2) ditch extra 'out' category
    T = struct2table(header.events);                        % creates a table of the structure 'events'
    T.code = categorical(T.code);                           % in order to work with strings --> convert to categorical
    sortedT = T(T.code == 'Stimulation', :); 
    sortedT.code = cellstr(sortedT.code);                   % turns the categorical data back to string cells
    header.events = table2struct(sortedT); 
    header.events = header.events';
    
    % verify the number of events
    message = [dataset_name ': ' num2str(size(sortedT, 1)) ' events found.'];
    disp(message)
    
    % save the data and the header as letswave files
    header.name = dataset_name;
    save([dataset_name '.mat'], 'data');
    save([dataset_name '.lw6'], 'header');
    disp('Datasets saved for letswave.')
end

%% functions
function [out_header,out_data]=EEG_import_MEGA(input_folder,session_number);
% Author : Andre Mouraux
% Institute of Neurosciences (IONS)
% Universite catholique de louvain (UCL)
% Belgium
% 
% modified by Dominika (13/08/2020)

out_header=[];
out_data=[];

message_string{1}=['Loading : ' input_folder];

%recording
recording=module_read_neurone(input_folder,session_number);

%prepare header
message_string{end+1}='Generating header.';
out_header.filetype='time_amplitude';
out_header.name=[input_folder,'_',num2str(session_number)];
out_header.tags={};
out_header.history(1).configuration=[];
out_header.datasize=double([1 length(recording.signalTypes) 1 1 1 recording.properties.length*recording.properties.samplingRate]);
out_header.xstart=1/recording.properties.samplingRate;
out_header.ystart=0;
out_header.zstart=0;
out_header.xstep=1/recording.properties.samplingRate;
out_header.ystep=1;
out_header.zstep=1;

%dummy chanloc
chanloc.labels='';
chanloc.topo_enabled=0;
chanloc.SEEG_enabled=0;
%set chanlocs
message_string{end+1}=['Importing ',num2str(length(recording.signalTypes)),' channel labels'];
for chanpos=1:length(recording.signalTypes);
    chanloc.labels=recording.signalTypes{chanpos};
    out_header.chanlocs(chanpos)=chanloc;
end;
            
%set events
out_header.events=[];
if isempty(recording.markers.index);
else
    numevents=length(recording.markers.index);
    message_string{end+1}=['Importing ',num2str(numevents),' events'];
    for eventpos=1:numevents;
        event.code='unknown';
        if ~isempty(recording.markers.type(eventpos));
            event.code=recording.markers.type{eventpos};
        end;
        if isnumeric(event.code);
            event.code=num2str(event.code);
        end;
        event.latency=recording.markers.time(eventpos);
        event.epoch=1;
        events(eventpos)=event;
    end;
    out_header.events=events;
end;

%data
message_string{end+1}=['Importing data (',num2str(out_header.datasize(6)),' samples, ',num2str(out_header.datasize(1)),' epoch(s))'];
out_data=zeros(round(out_header.datasize));
for k=1:length(recording.signalTypes)
    eval(['out_data(1,k,1,1,:)=squeeze(recording.signal.',recording.signalTypes{k},'.data);']);
end

out_header.datasize=size(out_data);
end
