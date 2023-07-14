function [rho,pval] = jerk_vs_tau_all(subject,regr)
%% heatmap of the joint distribution of jerk and tau: P(tau,J)
% regr: choose what type of regression you want to plot
figure;
whitebg([1 1 1]);

jerk = [];
tau = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        tau(count) = subject(i).trials(j).prs.tau;
        jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
    end
end
[rho,pval] = nancorr(tau',jerk');


plot(tau,jerk,'.','MarkerSize',4.5);
xlim([0 8]);ylim([0 150*10^2]);

title(['\tau vs Jerk^2, all trials - \rho = ' num2str(rho) ', pval = ' num2str(pval)]);ylabel('J^2');xlabel('\tau [s]');

choose_regress(tau,jerk,regr);

if length(subject) == 1
    suptitle(subject.name) 
end

%
%     x = [ceil(min(jerk)) ceil(max(jerk))];
%     y = [ceil(min(tau)) ceil(max(tau))];
%     plot(jerk,tau,'k.');grid on;
%     xlabel('Cumulative jerk');ylabel('\tau [s]');
