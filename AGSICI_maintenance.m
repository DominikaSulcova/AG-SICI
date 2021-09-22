for p = 1:length(position)
    for c = 1:length(current)
        for s = 1:length(subject)
            % define subject
            if subject(s) < 10
                subj = ['0' num2str(subject(s))];
            else
                subj = num2str(subject(s));
            end
            
            % average across intensities
            for e = 1:size(AGSICI_data, 5)
                data(1, e, 1, 1, 1, :) = squeeze(mean(AGSICI_data(p, c, :, s, e, :), 3));
            end
            
            % adapt and save header
            header.name = [subj ' ' position{p} ' ' current{c} ' all_intensities'];
            save([header.name '.lw6'], 'header')

            % save for LW              
            save([header.name '.mat'], 'data')
        end
    end
end
clear p c s data subj filename     