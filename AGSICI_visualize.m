%% AG-SICI: TEP AMPLITUDE & LATENCY GROUP PLOTS 
% Written by Dominika for AG-SICI project (2021)
% 
% 1) load data 
%       - in a matlab structure 'AGSICI_TEP_subject.mat' saved from the
%       previous processing step --> see script 'AGSICI_process.m'
% 2) plot mean amplitude        
%       - 3 repeated measures (intensities), 4 conditions (position x current)
%       - all in one figure, saves
% 3) plot mean latency
%       - as a boxplot 

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
subject = [1, 3:8, 10, 12];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = [100, 120, 140];
filename = 'AG-SICI.mat';

% statistics 
z = 1.96;
alpha = 0.05;
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

%% 1) load the data