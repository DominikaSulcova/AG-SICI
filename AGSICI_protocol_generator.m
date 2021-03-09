%% STIMULATION PROTOCOL GENERATOR
% Written by Dominika for AG-SICI project (2021)
% 
% The script creates text files compatible with MagPro X100 TMS stimulator. 
% According to the experimental design, each file consists of
%       - 75 single-pulse biphasic TMS stimuli 
%       - 3 intensities (3 x 25 stimuli) in semi-random order
%       - pseudo-random ISI <4, 6>s --> mean ISI 5s
% 
%% parameters
clear all, clc

target_cortex = 'AG';  
participant = {'03'};
n_protocol = 8;
amplitude = [100 120 140; 25 25 25];        % disponible output stim amplitudes 
                                            % --> r1: intensities (%rMT); r2: repetitions
ISI = [4000, 5000, 6000; 25, 25, 24];       % disponible ISI times
                                            % --> r1: times (ms); r2: repetitions
max_rep = 3;                                % maximum repetition of both ISI and intensity

%% write the file
% cycle through participants
repetitions = sum(amplitude(2,:));
for p = 1:length(participant)
    % outcome file
    outcome = [];
    
    % cycle through blocks
    for n = 1:n_protocol
        waiting_time = randomize(ISI, max_rep);
        intensity = randomize(amplitude, max_rep);
        filename = [target_cortex '_' participant{p} '_' num2str(n) '.CG3'];
        MagPro_initialize_bi(filename,repetitions);
        for i= 1:repetitions
            if i==1
                MagPro_single_bi(filename, i, 10000, intensity(i));
            else
                MagPro_single_bi(filename, i, waiting_time(i - 1), intensity(i));
            end        
        end
        % append to aoutcome
        outcome = cat(2, outcome, intensity); 
    end
    
    % save stim order
    filename_stim = ['stim_order_' participant{p}];
    save([filename_stim '.mat'], 'outcome')
end
clear p n 

%% functions
function seq = randomize(matrix, max_rep)
% creates a randomized vertical vector of N elements 
%   - elements specified by r1 of input matrix
%   - respective ns of repetition specified by r2 of inputmatrix 
%   --> N = r1(1)*r2(1) + r1(2)*r2(2)... + r1(n)*r2(n)
%   - maximum n of contiguous repetitions specified integer max_rep

seq = [];
for i = 1:length(matrix)
    single_isi = repelem(matrix(1, i), matrix(2, i)); 
    seq = [seq single_isi];
end

while true
     seq = seq(randperm(numel(seq)));           % shuffle randomly, may contain repetitions
     
     test_rep = [];
     for a = 1:(numel(seq) - (max_rep + 1))     % split the shuffled vector into sliding chunks                     
         b = seq(a:(a + max_rep));              % of n = max_rep + 1, and verify for each chunk 
         if range(b) == 0                       % that the elements are not all the same
             c = 1;                             % --> each chunk adds 1 element to a test_rep 
         else                                   % vector: 1 = all elements are the same, 0 = not the same
             c = 0;                                                          
         end
         test_rep = [test_rep c];
     end
     
     if all(test_rep == 0)                      % check if all elements in test_rep = 0        
        break;                                  % if the condition holds, exit loop
     end
end
seq = seq';                                     % the output vector should be in vertical form
end

function  [fileID] = MagPro_initialize_bi(filename,n_trials)
% opens and writes in a text file specified by filename
% = the first bloc of instructions required by the MagPro X100
% in order to read the experimental protocol file.
% NOTE : the extension given in the filename should be CG3

fileID = fopen(filename,'w');
fprintf(fileID, '[Model Option]\r\n');
fprintf(fileID, 'ModelOption=3\r\n');
fprintf(fileID, '[Main Menu]\r\n');
fprintf(fileID, 'Mode=0\r\n');
fprintf(fileID, 'Current Direction=0\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Inter Pulse Interval=100\r\n');
fprintf(fileID, 'Burst Pulses=0\r\n');
fprintf(fileID, 'Pulse BA Ratio=100\r\n');
fprintf(fileID, '[Timing Menu]\r\n');
fprintf(fileID, 'Timing Control=0\r\n');
fprintf(fileID, 'Rep Rate=3\r\n');
fprintf(fileID, 'Pulses in train=1\r\n');
fprintf(fileID, 'Number of Trains=1\r\n');
fprintf(fileID, 'Inter Train Interval=30\r\n');
fprintf(fileID, '[Trigger Menu]\r\n');
fprintf(fileID, 'Trig Output=1\r\n');
fprintf(fileID, 'Twin Trig output=1\r\n');
fprintf(fileID, 'Twin Trig Input=0\r\n');
fprintf(fileID, 'Polarity Input=1\r\n');
fprintf(fileID, 'Polarity output=1\r\n');
fprintf(fileID, 'Delay Input Trig=0\r\n');
fprintf(fileID, 'Delay Output Trig=0\r\n');
fprintf(fileID, '[Configuration Menu]\r\n');
fprintf(fileID, 'Charge Delay=300\r\n');
fprintf(fileID, 'Auto Discharge Time=60\r\n');
fprintf(fileID, 'Prior Warning Sound=0\r\n');
fprintf(fileID, 'Coil Type Display=1\r\n');
fprintf(fileID, '[MEP Menu]\r\n');
fprintf(fileID, 'Time=5 ms/div\r\n');
fprintf(fileID, 'Sensitivity=500 uV/div\r\n');
fprintf(fileID, 'Scroll=0 ms\r\n');
fprintf(fileID, 'Curve No=0\r\n');
fprintf(fileID, 'Baseline=1\r\n');
fprintf(fileID, 'Lower Freq=100 Hz\r\n');
fprintf(fileID, 'Upper Freq=5 kHz\r\n');
fprintf(fileID, 'Trigger mode=0\r\n');
fprintf(fileID, 'Size=0\r\n');
fprintf(fileID, 'On Top=1\r\n');
fprintf(fileID, '[Protocol Setup]\r\n');
fprintf(fileID,'Number of Lines=%d\r\n',n_trials);
fclose(fileID);
end

function [fileID]=MagPro_single_bi(filename, trial_number, isi, amplitude)
fileID = fopen(filename,'a');
fprintf(fileID, '[protocol Line %d]\r\n',trial_number);
fprintf(fileID, 'Delay=%d\r\n',isi);
fprintf(fileID, 'Amplitude A Gain=%d\r\n',amplitude/10);
fprintf(fileID, 'Mode=0\r\n');
fprintf(fileID, 'Current Direction=1\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Burst Pulses=2\r\n');
fprintf(fileID, 'Inter Pulse Interval=10000\r\n');
fprintf(fileID, 'BA Ratio=100\r\n');
fprintf(fileID, 'Repetition Rate=10\r\n');
fprintf(fileID, 'Train Pulses=1\r\n');
fclose(fileID);
end

