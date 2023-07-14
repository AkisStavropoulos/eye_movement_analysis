function taus_vs_jerk(trials,stimtype,regr)
% regr: type of regression, 'poly' for 2nd order polynomial, 'lin' for linear regression
    
stnames = fieldnames(stimtype);
dt = 1/60;
N = 50;
%% see distance error as a function of tau
figure;
colr = brewermap(length(stnames),'Dark2');
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    prs = [];
    taus = [];
    mc = [];
    jerk = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        indx = stimtype.(stnames{n})(k);
        jerk(k) = trials(indx).continuous.cumjerk(end);
    end
    [jerk,S] = sort(jerk);
    taus = taus(S);
    
    plot(taus,jerk,'.','Color',colr(n,:),'MarkerSize',5);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(taus,jerk,2);
        f = polyval(p,taus);
        [~,I] = sort(taus);
        hold on;plot(taus(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(taus,jerk);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'movavg')
        hold on;
        [X,I] = sort(taus);
        M = movmean(jerk(I),N,'omitnan');
        nanind = find(isnan(jerk(I)));
        X(nanind) = nan;
        plot(X,M,'Color',colr(n,:),'LineWidth',2.5);

    end
end
grid on;xlabel('\tau [s]');ylabel('jerkiness');xlim([0 8]);ylim([0 150*10^2]);
title('total jerk as a function \tau');legend(leg_input);

suptitle(trials(1).prs.subject);
