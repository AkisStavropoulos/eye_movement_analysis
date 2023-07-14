function dist_jerk_hist(subject,N,M)
%% Plot distribution of cumulative jerk of every trial and every subject conditioned on target distance
% N: histogram bins
% M: angle bins

% bin target angles
distbins = linspace(250,550,M+1);

figure;hold on;

colr = parula(M);
for n = 1:M
    jerk = [];
    count = 0;
    for i = 1:length(subject)
        for j = 1:length(subject(i).trials)
            if (subject(i).trials(j).prs.r_tar >= distbins(n)) && (subject(i).trials(j).prs.r_tar < distbins(n+1))
            count = count + 1;
            jerk(count) = subject(i).trials(j).continuous.jerk(end);
            end
        end
    end
    [y(n,:),x(n,:)] = hist(jerk,N);
    
    plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;
    
    leg_input{n} = ['d_t_a_r bin' num2str(n) ' - [' num2str(distbins(n)) ' - ' num2str(distbins(n+1)) ']'];
end
legend(leg_input{:});xlabel('cumulative jerk');ylabel('P(J^2|d_t_a_r)');
title('Histogram of cumulative jerk for different d_t_a_r bins');
if length(subject) == 1
    suptitle(subject.name) 
end
