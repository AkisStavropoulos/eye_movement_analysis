function ang_jerk_hist(subject,N,M)
%% Plot distribution of cumulative jerk of every trial and every subject conditioned on angle
% N: histogram bins
% M: angle bins

% bin target angles
angbins = linspace(0,40,M+1);

figure;hold on;
whitebg([.5 .5 .5]);
colr = parula(M);
for n = 1:M
    jerk = [];
    count = 0;
    for i = 1:length(subject)
        for j = 1:length(subject(i).trials)
            if (abs(subject(i).trials(j).prs.th_tar) >= angbins(n)) && (abs(subject(i).trials(j).prs.th_tar) < angbins(n+1))
            count = count + 1;
            jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
            end
        end
    end
    [y(n,:),x(n,:)] = hist(jerk,N);
    
    plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;
    
    leg_input{n} = ['\theta_t_a_r bin' num2str(n) ' - [' num2str(angbins(n)) ' - ' num2str(angbins(n+1)) ']'];
end
legend(leg_input{:});xlabel('cumulative jerk^2');ylabel('P(J^2|\theta_t_a_r)');
title('Histogram of cumulative jerk^2 for different \theta_t_a_r bins');
if length(subject) == 1
    suptitle(subject.name) 
end
