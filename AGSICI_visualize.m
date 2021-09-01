%% AG-SICI: TEP AMPLITUDE & LATENCY GROUP PLOTS 
% Written by Dominika for AG-SICI project (2021)
% 
% 1) load data 
%       - in a matlab structure 'AGSICI_TEP_subject.mat' saved from the
%       previous processing step --> see script 'AGSICI_process.m'
% 2) plot mean absolute TEP amplitude        
% 3) plot mean TEP latency
%       - 3 repeated measures (intensities), 4 conditions (position x current)
%       - all in one figure, saves in .fig and .png
% 4) peak-2-peak N40 amplitude
%       - calculates the amplitude for each subject/condition as absolute
%       N40 - absolute P25
%       - plots all conditions together as previous figures


%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
subject = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = [100, 120, 140];
peak = {'N45'};
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

% launch the figure counter
figure_counter = 1;

%% 1) prepare the data
% load data structure
load(filename, 'AGSICI_TEP_subject');

% extract data
% p = 1; c = 1; i = 1; k = 1; s = 1;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            for k = length(peak)
                % extract individual data 
                for s = 1:length(subject)
                    amplitude(p, c, i, k, s) = AGSICI_TEP_subject(s).amplitude(p, c, i, k);
                    latency(p, c, i, k, s) = AGSICI_TEP_subject(s).latency_real(p, c, i, k); 
                end
            end
        end
    end
end
clear p c i k s

%% 2) plot absolute TEP amplitude 
% set x
x = 1:length(intensity);

% plot separate figures for each peak
for k = 1:length(peak)
    % decide names
    fig_name = ['AGSICI_amp_abs_' peak{k}];
    fig_title = ['Mean absolute TEP amplitude : ' peak{k}];
    
    % launch the figure
    fig = figure(figure_counter); 
    hold on

    % plot the data
    counter = 1;
    for p = 1:length(position)
        for c = 1:length(current)
            for i = 1:length(intensity)
                % calculate the data and 95% CI
                y(i) = mean(amplitude(p, c, i, k, :));
                SEM(i) = std(amplitude(p, c, i, k, :)) / sqrt(length(subject));
                CI(i) = SEM(i) * z;
            end
            
            % plot
            perr(counter) = errorbar(x, y, SEM);

            % adjust parameters
            perr(counter).Color = colours(counter, :);
            perr(counter).LineWidth = 1.5;
            perr(counter).Marker = 'o';
            perr(counter).MarkerFaceColor = colours(counter, :);
            perr(counter).MarkerSize = 10;
            
            % update counter
            counter = counter + 1;
        end
    end
    
    % add features, adjust parameters
    set(gca, 'xtick', 1:length(intensity), 'xticklabel', intensity)
    set(gca, 'Fontsize', 14)
    title(fig_title, 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('stimulation intensity (% rMT)'); ylabel('TEP amplitude (\muV \pm SEM)');
    xlim([0.75, length(intensity) + 0.25])
    leg = legend(perr, {'along STS - normal' 'along STS - reversed' 'across STS - normal' 'across STS - reversed'});
    set(leg,'NumColumns',2,'Location','southeast');

    % save the figure       
    savefig([fig_name '.fig'])
    saveas(fig, [fig_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;   
end
clear x k fig p c i y SEM CI perr leg

%% 3) plot mean TEP latency
% set x
x = 1:length(intensity);

% plot separate figures for each peak
for k = 1:length(peak)
    % decide names
    fig_name = ['AGSICI_lat_' peak{k}];
    fig_title = ['Mean TEP latency : ' peak{k}];
    
    % launch the figure
    fig = figure(figure_counter); 
    hold on

    % plot the data
    counter = 1;
    for p = 1:length(position)
        for c = 1:length(current)
            for i = 1:length(intensity)
                % calculate the data and 95% CI
                y(i) = mean(latency(p, c, i, k, :));
                SEM(i) = std(latency(p, c, i, k, :)) / sqrt(length(subject));
                CI(i) = SEM(i) * z;
            end
            
            % plot
            perr(counter) = errorbar(x, y, SEM);

            % adjust parameters
            perr(counter).Color = colours(counter, :);
            perr(counter).LineWidth = 1.5;
            perr(counter).Marker = 'o';
            perr(counter).MarkerFaceColor = colours(counter, :);
            perr(counter).MarkerSize = 10;
            
            % update counter
            counter = counter + 1;
        end
    end
    
    % add features, adjust parameters
    set(gca, 'xtick', 1:length(intensity), 'xticklabel', intensity)
    set(gca, 'Fontsize', 14)
    title(fig_title, 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('stimulation intensity (% rMT)'); ylabel('peak latency (ms \pm SEM)');
    xlim([0.75, length(intensity) + 0.25])
    leg = legend(perr, {'along STS - normal' 'along STS - reversed' 'across STS - normal' 'across STS - reversed'});
    set(leg,'NumColumns',2,'Location','southeast');

    % save the figure       
    savefig([fig_name '.fig'])
    saveas(fig, [fig_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;   
end
clear x k fig p c i y SEM CI perr leg

%% 4) peak-2-peak N40
% calculate amplitudes for all conditions
% p = 1; c = 1; i = 1; s = 1;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)            
            for s = 1:length(subject)
                % p2p N40 = abs N40 - abs P25
                N40_p2p(p, c, i, s) = amplitude(p, c, i, 2, s) - amplitude(p, c, i, 1, s);
            end
        end
    end
end
clear p c i s

% set parameters
x = 1:length(intensity);
fig_name = 'AGSICI_N40_p2p';
fig_title = 'Peak-to-peak TEP amplitude : N40';


% launch the figure
fig = figure(figure_counter); 
hold on

% plot the data
counter = 1;
for p = 1:length(position)
    for c = 1:length(current)
        for i = 1:length(intensity)
            % calculate the data and 95% CI
            y(i) = mean(N40_p2p(p, c, i, :));
            SEM(i) = std(N40_p2p(p, c, i, :)) / sqrt(length(subject));
            CI(i) = SEM(i) * z;
        end

        % plot
        perr(counter) = errorbar(x, y, SEM);

        % adjust parameters
        perr(counter).Color = colours(counter, :);
        perr(counter).LineWidth = 1.5;
        perr(counter).Marker = 'o';
        perr(counter).MarkerFaceColor = colours(counter, :);
        perr(counter).MarkerSize = 10;

        % update counter
        counter = counter + 1;
    end
end

% add features, adjust parameters
set(gca, 'xtick', 1:length(intensity), 'xticklabel', intensity)
set(gca, 'Fontsize', 14)
title(fig_title, 'FontWeight', 'bold', 'FontSize', 16)
xlabel('stimulation intensity (% rMT)'); ylabel('TEP amplitude (\muV \pm SEM)');
xlim([0.75, length(intensity) + 0.25])
leg = legend(perr, {'along STS - normal' 'along STS - reversed' 'across STS - normal' 'across STS - reversed'});
set(leg,'NumColumns',2,'Location','southeast');

% save the figure       
savefig([fig_name '.fig'])
saveas(fig, [fig_name '.png'])

% update the counter
figure_counter = figure_counter + 1;   

clear x k fig p c i y SEM CI perr leg







