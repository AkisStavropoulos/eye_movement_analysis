function [rho,pval] = jerk_vs_mod_tau_all(subject,regr)
%% heatmap of the joint distribution of jerk and tau: P(tau,J)
% regr: choose what type of regression you want to plot
stimtype = fieldnames(subject(1).stimtype);
binnames = fieldnames(subject(1).tau.bins);
limnames = fieldnames(subject(1).tau.lims);
stimnames = {'vestibular' , 'visual' , 'combined'};
leg_inp = [stimnames ; stimnames];
leg_inp = leg_inp(:);

figure;
whitebg([1 1 1]);
colr = brewermap(length(stimnames),'Dark2');

for s = 1:length(stimtype)-1
    
    jerk = [];
    tau = [];
    count = 0;
    for i = 1:length(subject)
        
        for j = 1:length(subject(i).stimtype.(stimtype{s}))
            count = count + 1;
            indx = subject(i).stimtype.(stimtype{s})(j);
            tau(count) = subject(i).trials(indx).prs.tau;
            jerk(count) = subject(i).trials(indx).continuous.cumjerk(end);
        end
    end
    [rho(s),pval(s)] = nancorr(tau',jerk');
    
    
    plot(tau,jerk,'.','Color',colr(s,:),'MarkerSize',4.5);
    xlim([0 6]);ylim([0 150*10^2]);
    
    
    [XX,YY] = choose_regress(tau,jerk,regr);
    plot(XX,YY,'Color',colr(s,:),'LineWidth',2);
%     leg_inp{s} = stimnames{s};
end
corrstring = ['\rho = [' num2str(rho) '], pval = [' num2str(pval) ']'];
title({'\tau vs Jerk^2, all trials',corrstring});ylabel('J^2');xlabel('\tau [s]');
legend(leg_inp,'Location','northeast');

if length(subject) == 1
    suptitle(subject.name) 
end

%
%     x = [ceil(min(jerk)) ceil(max(jerk))];
%     y = [ceil(min(tau)) ceil(max(tau))];
%     plot(jerk,tau,'k.');grid on;
%     xlabel('Cumulative jerk');ylabel('\tau [s]');
