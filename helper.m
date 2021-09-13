%% merge epochs
pfx = 'AGSICI TEP'; % tracked across normal N100 centered
for p = 1:length(position)
    for c = 1:length(current)
        for k = 1:length(AGSICI_TEP_avg.peak)-1
            for g = 1:length(category)  
                % load original .mat data
                load([[pfx ' tracked ' position{p} ' ' current{c} ' ' AGSICI_TEP_avg.peak{k} ' ' category{g}] '.mat']) 
                data1 = data;
                
                % load new data
                load([pfx ' tracked_2 ' position{p} ' ' current{c} ' ' AGSICI_TEP_avg.peak{k} ' ' category{g}]) 
                data2 = data;
                
                % concatenate
                data_cat = cat(1, data1, data2);
                
                % average
                data = mean(data_cat);
                
                % save data 
                save([pfx ' tracked_all ' position{p} ' ' current{c} ' ' AGSICI_TEP_avg.peak{k} ' ' category{g} '.mat'], 'data')
            end
        end
    end
end

%% replace values by default
for s = 1:length(subject)
    for p = 1:length(position)
        for c = 1:length(current)
            s_index = length(subject_order) - length(subject) + s;
            AGSICI_TEP_subject(s_index).latency(p, c, 6) = 0.200;  
        end
    end
end

%% table work
table_N45 = AGSICI_outcome(457:end, :);
n_entries = length(position) * length(current) * length(intensity) * length(subject_order);
N45_final = table;
N45_final.subject = cell(n_entries, 1); 
N45_final.position = cell(n_entries, 1); 
N45_final.current = cell(n_entries, 1); 
N45_final.intensity = cell(n_entries, 1); 
N45_final.peak = cell(n_entries, 1); 
N45_final.amplitude = cell(n_entries, 1); 
N45_final.latency = cell(n_entries, 1); 
row_counter = 0;
row_counter_final = 1;
% a = 1; b = 1;
for a = 1:n_entries
    for b = 1:width(table_N45)
        N45_final{row_counter_final, b} = {table_N45{row_counter + a, b}};
        row_counter = row_counter + a;
    end
    row_counter_final = row_counter_final + 1;
end
AGSICI_outcome(457:end, :) = [];
AGSICI_outcome = [AGSICI_outcome; N45_final];
N45_final(1, 1) = cell2mat(N45_final(1, 1));
AGSICI_outcome_N45 = N45_final;
save(filename, 'AGSICI_outcome_N45', 'AGSICI_outcome', '-append');