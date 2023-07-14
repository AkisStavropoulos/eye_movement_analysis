function tau_jerk_strat_hist(subject,N,thresh)
%% Plot distribution of cumulative jerk of every trial and every subject

figure;
taubinnames = fieldnames(subject(1).tau.bins);
taulimnames = fieldnames(subject(1).tau.lims);

strategies = 2;
jerkbins = [0 thresh 10^4];
leg_loc = {'northeast' , 'northwest'};
for s = 1:strategies
    subplot(strategies,1,s);hold on;
    colr = parula(length(taubinnames));
    for n = 1:length(taubinnames)
        jerk = [];
        count = 0;
        for i = 1:length(subject)
            for j = 1:length(subject(i).tau.bins.(taubinnames{n}))
                indx = subject(i).tau.bins.(taubinnames{n})(j);
                if (subject(i).trials(indx).continuous.jerk(end) >= jerkbins(s)) && (subject(i).trials(indx).continuous.jerk(end) < jerkbins(s+1))
                    count = count + 1;
                    jerk(count) = subject(i).trials(indx).continuous.jerk(end);
                end
            end
        end
        [y(n,:),x(n,:)] = hist(jerk,N);
        
        plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;vline(thresh);
        
        leg_input{n} = ['\tau ' taubinnames{n} ' - [' num2str(subject(i).tau.lims.(taulimnames{n})) ']'];
    end
    legend(leg_input{:},'Location',leg_loc{s});
    xlim([0 150]);
    xlabel('cumulative jerk');title(['Histogram of J|\tau for strategy ' num2str(s) ', J_t_h_r_e_s_h = ' num2str(thresh)]);
end
if length(subject) == 1
    suptitle(subject.name) 
end
