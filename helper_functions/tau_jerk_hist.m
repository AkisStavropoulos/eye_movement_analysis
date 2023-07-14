function tau_jerk_hist(subject,N)
%% Plot distribution of cumulative jerk of every trial and every subject

figure;hold on;
whitebg([.5 .5 .5]);
taubinnames = fieldnames(subject(1).tau.bins);
taulimnames = fieldnames(subject(1).tau.lims);
colr = parula(length(taubinnames));
for n = 1:length(taubinnames)
    jerk = [];
    count = 0;
    for i = 1:length(subject)
        for j = 1:length(subject(i).tau.bins.(taubinnames{n}))
            count = count + 1;
            indx = subject(i).tau.bins.(taubinnames{n})(j);
            jerk(count) = subject(i).trials(indx).continuous.cumjerk(end);
        end
    end
    [y(n,:),x(n,:)] = hist(jerk,N);
    
    plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;
    
    leg_input{n} = ['\tau ' taubinnames{n} ' - [' num2str(subject(i).tau.lims.(taulimnames{n})) ']'];
end
xlim([0 150*102]);
legend(leg_input{:});xlabel('cumulative jerk^2');ylabel('P(J^2|\tau)');
title('Histogram of cumulative jerk^2 for different tau bins');
if length(subject) == 1
    suptitle(subject.name) 
end
