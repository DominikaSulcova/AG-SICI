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