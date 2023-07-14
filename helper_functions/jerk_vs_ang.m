function jerk_vs_ang(trials,stimtype,regr)
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
    ang = [];
    mc = [];
    jerk = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    ang = abs([prs.th_tar]);
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        indx = stimtype.(stnames{n})(k);
        jerk(k) = trials(indx).continuous.cumjerk(end);
    end
    [ang,S] = sort(ang);
    jerk = jerk(S);
    
    plot(ang,jerk,'.','Color',colr(n,:),'MarkerSize',5);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(ang,jerk,2);
        f = polyval(p,ang);
        [~,I] = sort(ang);
        hold on;plot(ang(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(ang,jerk);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'movavg')
        hold on;
        M = movmean(jerk,N,'omitnan');
        nanind = find(isnan(jerk));
        X = ang;
        X(nanind) = nan;
        plot(X,M,'Color',colr(n,:),'LineWidth',2.5);
    end
end
grid on;xlabel('target angle [deg]');ylabel('jerkiness');xlim([0 40]);ylim([0 100*10^2]);
title('total jerk as a function target angle');legend(leg_input);

suptitle(trials(1).prs.subject);
