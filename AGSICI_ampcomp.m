%% AG-SICI: AMPLITUDE COMPARISON
% written by Dominika for AG-SICI project (2021)
% 
% 1) loads amplitudes of the chosen peak 
%   - from M1 TEPs: GABA-AD YC both sessions, baseline --> 40 entries
%   - from AG TEPs: AG-SICI P1 --> 19 entries
% 2) M1 - 80 / 120 %rMT - boxplot
% 3) AG - 100 / 120 / 140 %rMT - boxplot

%% parameters
clear all, clc

% datafile location
datafile_M1 = 'E:\UCL\O365G-NOCIONS - People\dsulcova\GABA-AD_results\Data\GABA_results';
datafile_AG = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\AG-SICI_plus';

% choose peak to compare
peak = 1;
peaks = [{'P30' 'N45' 'P60' 'N100' 'P180'};...
    {'P25' 'N45' 'P75' 'N120' 'P200'}];
intensities = [100 120 140];

% graphics
figure_counter = 1;
load('E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\colours')

%% 1) data extraction
% ----- M1 data -----
load(datafile_M1)

% launch outcome variables
amp_M1_80 = []; amp_M1_120 = [];

% extract amplitude
for a = 1:length(GABA_results)
    for b = 1:length(GABA_results(a).TEP)
        % sub-threshold stimulus --> CS = 80%rMT
        amp_M1_80 = [amp_M1_80; GABA_results(a).TEP(b).CS.pre(1, peak)]; 
        
        % sub-threshold stimulus --> TS = 120%rMT
        amp_M1_120 = [amp_M1_120; GABA_results(a).TEP(b).TS.pre(1, peak)]; 
    end
end

% ----- AG data -----
load(datafile_AG, 'AGSICI_TEP_subject')

% launch outcome variables 
amp_AG_100 = []; amp_AG_120 = []; amp_AG_140 = [];

% extract amplitude
for a = 1:length(AGSICI_TEP_subject)
    for b = [1 2]
        for c = [1 2]
            % 100%rMT
            amp_AG_100(a, (b - 1)*2 + c) = AGSICI_TEP_subject(a).amplitude(b, c, 1, peak); 
            
            % 120%rMT
            amp_AG_120(a, (b - 1)*2 + c) = AGSICI_TEP_subject(a).amplitude(b, c, 2, peak); 
            
            % 140%rMT
            amp_AG_140(a, (b - 1)*2 + c) = AGSICI_TEP_subject(a).amplitude(b, c, 3, peak); 
        end
    end
end

clear AGSICI_TEP_subject GABA_results a b c

%% 2) plot mean values - M1
% get data
data_visual = [amp_M1_80, amp_M1_120];

% launch the figure
fig = figure(figure_counter);
boxplot(data_visual, 'colors', [0 0.45 0.75; 0.88 0.2 0.3])
set(gca, 'xtick', 1:2, 'xticklabel', {'80 %rMT' '120 %rMT'})
set(gca, 'Fontsize', 16)
ylim([-2, 12])
title(['M1 ' peaks{1, peak}], 'FontWeight', 'bold', 'FontSize', 18)
text(0.75, 11, sprintf('mean = %0.3f ?V', mean(data_visual(:, 1))), 'fontsize', 12)
text(1.75, 11, sprintf('mean = %0.3f ?V', mean(data_visual(:, 2))), 'fontsize', 12)

% save the figure
figure_name = 'ampcomp_M1';
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;

%% 3) plot mean values - AG
for i = 1:length(intensities)
    % get data
    statement = ['data_visual = amp_AG_' num2str(intensities(i)) ';'];
    eval(statement)

    % launch the figure
    fig = figure(figure_counter);
    boxplot(data_visual, 'colors', colours)
    set(gca, 'xtick', 1:4, 'xticklabel', {'along N' 'along R' 'across N' 'across R'})
    set(gca, 'Fontsize', 16)
    ylim([-2, 6])
    title(['AG ' peaks{2, peak} ' - ' num2str(intensities(i)) ' %rMT'], 'FontWeight', 'bold', 'FontSize', 18)
    for a = 1:size(data_visual, 2)
        text(a - 0.25, 5.5, sprintf('%0.3f ?V', mean(data_visual(:, a))), 'fontsize', 12)
    end

    % save the figure
    figure_name = ['ampcomp_AG' num2str(intensities(i))];
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;
end


