function [rho,pval] = jerk_vs_ang_all(subject,regr)
%% jerk and target angle for all subjects
% regr: choose what type of regression you want to plot
figure;hold on;
whitebg([1 1 1]);

jerk = [];
theta = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        theta(count) = abs(subject(i).trials(j).prs.th_tar);
        jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
        if theta(count) > 40
            theta(count) = nan;
        end
    end
end

[rho,pval] = nancorr(theta',jerk');

plot(theta,jerk,'.','MarkerSize',4.5),
xlim([0 40]);ylim([0 150*10^2]);

title(['\theta_t_a_r_g_e_t vs Jerk^2, all trials, \rho = ' num2str(rho) ', pval = ' num2str(pval)]);ylabel('J^2');xlabel('\theta_t_a_r [deg]');

choose_regress(theta,jerk,regr);

% calculate moving mean to find slope that divides strategies
% N = 10000;
% [theta,S] = sort(theta);
% jerk = jerk(S);
% M = movmedian(jerk,N,'omitnan');
% nanind = find(isnan(jerk));
% X = theta;
% X(nanind) = nan;
% plot(X,M-1000,'r','LineWidth',2.5);
% 
if length(subject) == 1
    suptitle(subject.name) 
end

%
%     x = [ceil(min(jerk)) ceil(max(jerk))];
%     y = [ceil(min(tau)) ceil(max(tau))];
%     plot(jerk,tau,'k.');grid on;
%     xlabel('Cumulative jerk');ylabel('\tau [s]');
