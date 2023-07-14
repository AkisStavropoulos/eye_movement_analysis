function [RT_avg,RT_sig] = jerk_onset_hist(subject,N)
%% Plot distribution of cumulative jerk of every trial and every subject
dt = 1/60;
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        
        if ~isnan(subject(i).trials(j).continuous.jerk)
            count = count +1;
            indx = find(subject(i).trials(j).continuous.jerk > 0 ,1);
            onset(count) = indx*dt;
            
            if onset(count) > 5
                onset(count) = nan;
            end
            
        end
    end
end
RT_avg = nanmean(onset - 1.13);
RT_sig = nanstd(onset - 1.13);

figure;
hist(onset,N);xlabel('time since target appearance [s]');title('Distribution of time of subjects motion start');hold on;
xlim([0,5]);ylim([0,250]);
if length(subject) == 1
    suptitle(subject.name) 
end


%%
% strategies = 2;
% jerkbins = [0 thresh 10^4];
% leg_loc = {'northeast' , 'northwest'};
% for s = 1:strategies
%     subplot(strategies,1,s);hold on;
%     colr = parula(length(taubinnames));
%     for n = 1:length(taubinnames)
%         jerk = [];
%         count = 0;
%         for i = 1:length(subject)
%             for j = 1:length(subject(i).tau.bins.(taubinnames{n}))
%                 indx = subject(i).tau.bins.(taubinnames{n})(j);
%                 if (subject(i).trials(indx).continuous.jerk(end) >= jerkbins(s)) && (subject(i).trials(indx).continuous.jerk(end) < jerkbins(s+1))
%                     count = count + 1;
%                     jerk(count) = subject(i).trials(indx).continuous.jerk(end);
%                 end
%             end
%         end
%         [y(n,:),x(n,:)] = hist(jerk,N);
%         
%         plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;vline(thresh);
%         
%         leg_input{n} = ['\tau ' taubinnames{n} ' - [' num2str(subject(i).tau.lims.(taulimnames{n})) ']'];
%     end
%     legend(leg_input{:},'Location',leg_loc{s});
%     xlim([0 150]);
%     xlabel('cumulative jerk');title(['Histogram of J|\tau for strategy ' num2str(s) ', J_t_h_r_e_s_h = ' num2str(thresh)]);
% end
% if length(subject) == 1
%     suptitle(subject.name) 
% end
