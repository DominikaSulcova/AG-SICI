%% EXPERIMENTAL SESSION AG-SICI
% Written by Dominika for AG-SICI project (2021)
% 
% The script runs along the experimental session
%   1) initializes a logfile for current participant, asks to fill in 
%   all relevant information
%   2) calculates stimulation intensities based on rMT

%% session info
clear all; clc

% session
prompt = {'Study:', 'Subject:', 'Session:'};
dlgtitle = 'Session information';
dims = [1 45];
definput = {'AG-SICI', 'P1 - 01', '1'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% coil placement
session_info{4} = questdlg('Coil placement?', 'Experimental condition',...
    'along axis', 'across axis', 'along axis');

% stimulation
intensity = [100 120 140];

%% stimulation intensity
% M1 - rMT
session_info{5} = cell2num(inputdlg('rMT at the beginning of the sesion:', 'rMT'));

% calculate intensities
intensity = round(intensity/100 * session_info{5});
session_info{6} = intensity;
disp('Stimulation intensities:')
disp(['100 %rMT --> ' num2str(intensity(1)) ' %MSO'])
disp(['120 %rMT --> ' num2str(intensity(2)) ' %MSO'])
disp(['140 %rMT --> ' num2str(intensity(3)) ' %MSO'])

% closest electrodes 
prompt = {'Closest electrodes:'};
dlgtitle = 'AG target';
dims = [1 50];
definput = {'27, 7, 15'};
answer = cell2num(inputdlg(prompt,dlgtitle,dims,definput));
clear prompt dlgtitle dims definput

% associate numbers with electrode names
labels = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F7' 'F8' 'T7' 'T8' 'P7' 'P8'...
    'Fz' 'Cz' 'Pz' 'Iz' 'FC1' 'FC2' 'CP1' 'CP2' 'FC5' 'FC6' 'CP5' 'CP6' 'P5' 'P6' 'C1' 'C2'};
for e = 1:length(answer)
    session_info{end + 1} = [num2str(answer(e)) ' (' labels{answer(e)} ')'];
end

% create logfile
filename = initialize_logfile(session_info);

%% functions
function  filename = initialize_logfile(session_info)
% get the date
c = fix(clock);
date = [num2str(c(3)) '/' num2str(c(2)) '/' num2str(c(1))];
clear c

% define the filename
filename = [session_info{1} '_' session_info{2} '_S' session_info{3} '.txt'];

% count eletrodes
for e = 1:numel(session_info) - 6
    electrodes{e} = session_info{6 + e};
end

% write the file 
fileID = fopen(filename,'w');
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, ['study: ', session_info{1}, '\r\n']); 
fprintf(fileID, ['subject: ', session_info{2}, '\r\n']);
fprintf(fileID, ['session: S', session_info{3}, '\r\n']);
fprintf(fileID, ['coil placement: ', session_info{4}, '\r\n']);
fprintf(fileID, ['recorded: ', date, '\r\n']);
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, 'DATA ACQUISITION\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '- TMS: MagVenture MagPro X100 + Visor2 neuronavigation based on individual T1 MRI image\r\n');
fprintf(fileID, '     - over left AG\r\n')
fprintf(fileID, '     - biphasic sin pulse with 2 different current directins\r\n');
fprintf(fileID, '     - 3 tested intensities: 100 - 120 - 140 %%rMT\r\n');
fprintf(fileID, '- EEG: Bittium system - TESLA amplifiers, NeurOne recording software, 20kHz sampling rate, 3500Hz LP device filter\r\n');
fprintf(fileID, '     - 32 electrodes (10-20), M1/M2 not recorded - replaced by P5/P6 \r\n');
fprintf(fileID, '     - referenced to AFz, ground at F6\r\n');
fprintf(fileID, '- TMS-EEG is recorded in 8 blocks, each consisting of 75 stimuli with variable ISI 4 - 6s (~6.5 min)\r\n');
fprintf(fileID, '- every block consists of 25 stimuli of each tested intensity\r\n');
fprintf(fileID, '- current direction in the coil changes every other block\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, ['rMT : ' num2str(session_info{5}) ' %%MSO\r\n']);
fprintf(fileID, ['stimulation intensities: ' num2str(session_info{6}(1)) ' - ' num2str(session_info{6}(2)) ' - ' num2str(session_info{6}(3)) ' %%MSO\r\n']);
fprintf(fileID, '\r\n');
fprintf(fileID, ['AG closest electrodes: ' electrodes{:} '\r\n']);
fprintf(fileID, '\r\n');
fprintf(fileID, 'notes:\r\n');
fprintf(fileID, '\r\n');
fclose(fileID);
end
