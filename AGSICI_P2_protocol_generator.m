%% STIMULATION PROTOCOL GENERATOR
% Written by Dominika for AG-SICI project (2021)
% 
% The script creates text files compatible with MagPro X100 TMS stimulator. 
% According to the experimental design, each file consists of
%       - 77 biphasic TMS stimuli of 7 conditions, each per 10 stimuli
%       - in semi-random order, ITI <4, 6>s --> mean ISI 5s 
%       - 4 single-pulse conditions:
%               - TS: 120% rMT
%               - CS1: 60% rMT
%               - CS2: 80% rMT
%               - CS3: 100% rMT
%       - 3 paired-pulse conditions, ISI 2.5ms:
%               - CS1-TS, CS2-TS, CS3-TS

%% parameters
clear all, clc

target_cortex = 'AG';  
participant = {'15'};
n_block = 9;
rep_cond = 11;
amp_TS = 120;
amp_CS = [60, 80, 100]; 
ITI = [4000, 5000, 6000; 26, 25, 25];       % disponible ISI times
                                            % --> r1: times (ms); r2: repetitions
ISI = 2500;
max_rep = 3;                                % maximum repetition of both ISI and intensity

%% write the file
% prepre input
n_cond = length(amp_TS) + length(amp_CS) + length(amp_TS)*length(amp_CS);
rep_block = rep_cond * n_cond;
cond = [1:n_cond; repelem(rep_cond, n_cond)];
amp = {amp_TS, amp_CS(1), amp_CS(2), amp_CS(3), [amp_CS(1), amp_TS], [amp_CS(2), amp_TS], [amp_CS(3), amp_TS]};

% cycle through participants
for p = 1:length(participant)
    % set up the outcome file
    outcome = [];

    % cycle through blocks
    for n = 1:n_block
        % initiate the text file
        filename = [target_cortex '_SICI_' participant{p} '_'  num2str(n) '.CG3'];
        MagPro_initialize_bi(filename, rep_block);
        
        % randomize ITI and conditions        
        waiting_time = randomize(ITI, max_rep);
        order_cond = randomize(cond, max_rep);
        
        % append condition order to outcome
        outcome = cat(2, outcome, order_cond); 

        % cycle through trials
        for t = 1:rep_block 
            % identify condition
            amp_t = amp{order_cond(t)};
            
            % single pulse
            if length(amp_t) == 1
                if t == 1
                    MagPro_single_bi_R(filename, t, 10000, amp_t);
                else
                    MagPro_single_bi_R(filename, t, waiting_time(t - 1), amp_t);
                end
                
            % paired pulse
            else
                if t == 1
                    MagPro_paired_bi_R(filename, t, 10000, amp_t(1), amp_t(2), ISI);
                else
                    MagPro_paired_bi_R(filename, t, waiting_time(t - 1), amp_t(1), amp_t(2), ISI);
                end
            end           
        end
    end

    % save stim order
    filename_stim = ['SICI_stim_order_' participant{p}];
    save([filename_stim '.mat'], 'outcome')
end
clear p n t 

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
function  [fileID] = MagPro_initialize_bi(filename, n_trials)
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
function [fileID] = MagPro_single_bi_R(filename, trial_number, ITI, amplitude)
fileID = fopen(filename,'a');
fprintf(fileID, '[protocol Line %d]\r\n',trial_number);
fprintf(fileID, 'Delay=%d\r\n',ITI);
fprintf(fileID, 'Amplitude A Gain=%d\r\n',amplitude/10);
fprintf(fileID, 'Mode=0\r\n');
fprintf(fileID, 'Current Direction=2\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Burst Pulses=2\r\n');
fprintf(fileID, 'Inter Pulse Interval=10000\r\n');
fprintf(fileID, 'BA Ratio=100\r\n');
fprintf(fileID, 'Repetition Rate=10\r\n');
fprintf(fileID, 'Train Pulses=1\r\n');
fclose(fileID);
end
function [fileID] = MagPro_paired_bi_R(filename, trial_number, ITI, amplitude_CS, amplitude_TS, ISI)
fileID = fopen(filename,'a');
BAratio = round((amplitude_TS/amplitude_CS)*100);
if mod(BAratio,5) > 0
    if mod(BAratio,5) < 3
        BAratio = BAratio - mod(BAratio,5);
    end
    if mod(BAratio,5) >= 3
        BAratio = BAratio + (5 - mod(BAratio,5));
    end
end 
BAratio = round((amplitude_TS/amplitude_CS)*100);
fprintf(fileID, '[protocol Line %d]\r\n',trial_number);
fprintf(fileID, 'Delay=%d\r\n',ITI);
fprintf(fileID, 'Amplitude A Gain=%d\r\n',amplitude_CS/10);
fprintf(fileID, 'Mode=2\r\n');
fprintf(fileID, 'Current Direction=2\r\n');
fprintf(fileID, 'Wave Form=1\r\n');
fprintf(fileID, 'Burst Pulses=2\r\n');
fprintf(fileID, 'Inter Pulse Interval=%d\r\n',ISI);
fprintf(fileID, 'BA Ratio=%d\r\n', BAratio);
fprintf(fileID, 'Repetition Rate=10\r\n');
fprintf(fileID, 'Train Pulses=1\r\n');
fclose(fileID)
end
