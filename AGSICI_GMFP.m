%% AG-SICI: GMFP - STRENGTH OF THE RESPONSE
% written by Dominika for AG-SICI project (2021)
% 
% 1) calculate GMFP for each subject/condition
%       - pools intensities together, averages
%       - calculates GMFP as standard deviation across electrodes
%       - pools together both datasets of the same coil position 
%           --> to contrast the differences in global activation
% 2) plot mean GMFP, visualize AUC 
%       - plots mean GMFP for both coil positions
%       - marks differences in AUC
%       - adds mean topoplots for chosen peaks
%       - saves the figure in the target folder
% 3) calculate AUC
%       - for each subject/coil position
%       - time window of interest --> see parameters
%       - plots the AUC values - paired boxplot, marks outliers
%       - saves the figure in the target folder

%% parameters
clear all, clc

% datafile location
data_path = 'E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Results\AG-SICI_plus.mat';
data_var = 'AGSICI_data';

% dataset
subject = [1, 3:18, 20, 21];
position = {'along' 'across'}; 
current = {'normal' 'reversed'};
intensity = {'stim_100' 'stim_120' 'stim_140'};

% load a header
load('E:\UCL\O365G-NOCIONS - People\dsulcova\AG-SICI\Data\P1\Processed data\avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1 03 along reversed stim_120.lw6', '-mat')

% times of interest
TOI = [0.015, 0.085];
TOI_peaks = [0.023 0.047 0.075];
peaks = {'P25' 'N45' 'P75'};

% visualization 
figure_counter = 8;
colours = [0 0.48 0.74; 0.99 0.18 0.18; 0.96 0.68 0.68];
lines = {':' '-'};
time_window = [-0.05, 0.3]; 
x = [time_window(1):header.xstep:time_window(2)];

% create output folders
foldername = 'AG-SICI_plus_figs';
folderpath = [pwd '\' foldername];
if ~exist(folderpath) 
    mkdir(folderpath)
end     

%% 1) calculate GMFP 
% load the data
load(data_path, data_var)

% calculate GMFP for each subject/condition
for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:size(AGSICI_data, 4)    
            % average across intensities
            for t = 1:size(AGSICI_data, 6)
                gmfp_data(p, c, s, :, t) = mean(squeeze(AGSICI_data(p, c, :, s, :, t)), 1);
            end

            % calculate GMFP
            GMFP(p, c, s, :) = std(squeeze(gmfp_data(p, c, s, 1:size(gmfp_data, 4) - 1, :)), 1);
        end
    end
end
clear p c i s t gmfp_data

% merge 'current' conditions 
for p = 1:length(position)
    for s = 1:size(AGSICI_data, 4)   
        GMFP_position(p, s, :) = mean(squeeze(GMFP(p, :, s, :)), 1);
    end
end
clear p s 

%% 2) plot - contrast coil position
% average data across subjects                
data_visual = squeeze(mean(squeeze(GMFP_position(:, :, :)), 2));

% launch the figure
fig = figure(figure_counter);
h_axis(1) = subplot(4, length(TOI_peaks) + 1 , 1 : 2*(length(TOI_peaks) + 1));
title(sprintf('GMFP: different coil position %d - %dms', TOI(1)*1000, TOI(2)*1000), 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12)
xlabel('time (s)')
ylabel('GMFP (\muV)')
hold on

% set limits of the figure
plot(x, data_visual(2, :))
yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
xl = get(gca, 'xlim');
cla(h_axis(1))
ylim(yl); xlim(xl);

% plot the TOI
rectangle('Position', [TOI(1), yl(1)+0.01, TOI(2) - TOI(1), yl(2) - yl(1)], 'FaceColor', colours(3, :), 'EdgeColor', 'none')

% plot interpolated part
rectangle('Position', [0, yl(1)+0.01, 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none') 

% plot the shading
F = fill([x fliplr(x)],[data_visual(1, :) fliplr(data_visual(2, :))], ...
    colours(2, :) , 'linestyle', 'none');

% plot GMFP
for p = 1:length(position)
    % plot signal
    P(p) = plot(x, data_visual(p, :), 'Parent', h_axis(1), ...
        'Color', colours(1, :), 'LineStyle', lines{p}, 'LineWidth', 2.5);
end

% mark TMS stimulus
line([0, 0], yl, 'Parent', h_axis(1), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)  

% add legend
lgd = legend(P, {'along STS' 'across STS'}, 'Location', 'southeast');
lgd.FontSize = 12; title(lgd, 'coil position')

% add topoplots
for p = 1:length(position)
    % choose data for topoplots 
    for e = 1:size(AGSICI_data, 5)
        for i = 1:size(AGSICI_data, 6)
            data_topoplot(1, e, 1, 1, 1, i) = mean(squeeze(AGSICI_data(p, :, :, :, e, i)), 'all');
        end
    end

    % add topoplots for all peaks
    for k = 1:length(TOI_peaks)
        n = 1 + (p-1)*length(TOI_peaks) + k;
        
        % plot the topoplot
        h_axis(n) = subplot(4, length(TOI_peaks) + 1, 2*(length(TOI_peaks) + 1) + (p-1)*(length(TOI_peaks)+1) + k);
        topo_plot(header, data_topoplot, TOI_peaks(k), time_window(1), [-2, 2])

        % modifies the layout
        if p == 1
            pos = get(h_axis(n), 'Position');
            pos(2) = pos(2) - 0.04;
            set(h_axis(n), 'Position', pos);
        else
            text(-0.2, -0.7, peaks{k}, 'FontWeight', 'bold', 'FontSize', 12)
        end

        if k == length(TOI_peaks)
            if p == 1
                descript = 'along STS';
            else
                descript = 'across STS';
            end                    
            text(0.8, 0, descript, 'Parent', h_axis(n), ...
                'Color', [0 0 0], 'FontSize', 12)
        end
    end
end
hold off

% change figure position
fig.Position = [-1200 200 680 600];

% save figure
if length(peaks)>1
    figure_name = sprintf('AGSICI_GMFP_%s-%s', peaks{1}, peaks{end});
else
    figure_name = sprintf('AGSICI_GMFP_%s', peaks{1});
end
savefig([folderpath '\' figure_name '.fig'])
saveas(fig, [folderpath '\' figure_name '.png'])

% update figure counteer
figure_counter = figure_counter + 1;

clear fig h_axis yl xl F P lgd p e i k n pos descript data_topoplot figure_name

%% 3) calculate and plot AUC 
% define TOI limits
x_start = (TOI(1) - time_window(1))/header.xstep;
x_end = (TOI(2) - time_window(1))/header.xstep;

% calculate AUC for each subject/coil position
for p = 1:length(position)
    for s = 1:size(AGSICI_data, 4)    
        AUC_position(p, s) = trapz(squeeze(GMFP_position(p, s, x_start:x_end)));
    end
end
AUC_position = AUC_position';
clear p s 

% plot mean values
fig = figure(figure_counter);        
boxplot(AUC_position, 'color', colours([1 2], :))
hold on

% plot the lines 
for s = 1:size(AUC_position, 1)
    P(s) = plot([1 2], AUC_position(s, [1 2]), '-o',...
        'Color', [0.75, 0.75, 0.75],...
        'MarkerSize', 10,...
        'MArkerEdge', 'none');
    hold on
end

% plot the markers
for b = [1 2]
    scat(b) = scatter(repelem(b, size(AUC_position, 1)), AUC_position(:, b),...
        75, colours(b, :), 'filled');
    hold on
end

% add parameters
set(gca, 'xtick', [1 2], 'xticklabel', {'along STS' 'across STS'})
set(gca, 'Fontsize', 12)
title(sprintf('GMFP: area under the curve %d - %dms', TOI(1)*1000, TOI(2)*1000), 'FontWeight', 'bold', 'FontSize', 16)
xlabel('coil position'); ylabel('AUC');

% mark outliers
h_out = flipud(findobj(gcf,'tag','Outliers'));
for h = 1:length(h_out)
    x_out =  get(h_out(h), 'XData');
    y_out =  get(h_out(h), 'YData');
    for i = 1:length(x_out)
        if ~(isnan(x_out(i)))
            index_out(h, i) = find(AUC_position(:, h) == y_out(i));
            text(x_out(i) + 0.05, double(y_out(i)), sprintf('subject %d', subject(index_out(h, i))))
        end
    end
end
hold off

% save figure
if length(peaks)>1
    figure_name = sprintf('AGSICI_AUC_%s-%s', peaks{1}, peaks{end});
else
    figure_name = sprintf('AGSICI_AUC_%s', peaks{1});
end
savefig([folderpath '\' figure_name '.fig'])
saveas(fig, [folderpath '\' figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;    

clear x_start x_end x_out y_out h_out P scat fig p s b h i  

%% functions
function topo_plot(header, data, x_pos, x_start, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - x_start)/header.xstep);
vector = data(1, :, 1, 1, 1, x_visual);

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i=1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i=i+1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end