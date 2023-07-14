function errors_vs_jerk(trials,stimtype,regr)
% errors vs cum squared jerk
% regr: type of regression, 'poly' for 2nd order polynomial, 'lin' for linear regression
    
stnames = fieldnames(stimtype);
dt = 1/60;
%% see distance error as a function of tau
figure;subplot(2,1,1);
colr = brewermap(length(stnames),'Dark2');
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
%% Compute cumulative squared jerk
for n = 1:length(stnames)-1
    stats = [];
    d_err = [];
    mc = [];
    jerk = [];
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        if ~isnan(mc(k).JS_X_Raw)
        jerk(k) = sum((diff(mc(k).JS_X_Raw)./dt).^2 + (diff(mc(k).JS_Yaw_Raw)./dt).^2);
        else
            jerk(k) = nan;
        end
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    d_err = [stats.d_err];
    [jerk,S] = sort(jerk);
    d_err = d_err(S);
    
    plot(jerk,d_err,'.','Color',colr(n,:),'MarkerSize',8);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(jerk,d_err,2);
        f = polyval(p,jerk);
        [~,I] = sort(jerk);
        hold on;plot(jerk(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(jerk,d_err);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
end
grid on;xlabel('jerkiness');ylabel('error (cm)');ylim([-500 500]);xlim([0 10^8]);
title('Distance error as a function total jerk');legend(leg_input);
% absolute distance error
% subplot(2,2,2);
% for n = 1:length(stnames)-1
%     stats = [];
%     d_err = [];
%     prs = [];
%     taus = [];
%     prs = [trials(stimtype.(stnames{n})).prs];
%     taus = [prs.tau];
%     stats = [trials(stimtype.(stnames{n})).stats];
%     d_err = abs([stats.d_err]);
%     [x,y,b,c] = find_regress(taus,d_err);
%     plot(taus,d_err,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
% end
% grid on;xlabel('\tau (s)');ylabel('ABS error (cm)');
% title('ABSOLUTE distance error as a function \tau');legend(leg_input);

suptitle(trials(1).prs.subject);
%% see angular error as a function of tau
p = [];
f = [];
subplot(2,1,2);
colr = brewermap(length(stnames),'Dark2');
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    stats = [];
    th_err = [];
    mc = [];
    jerk = [];
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        jerk(k) = sum((diff(mc(k).JS_X_Raw)./dt).^2 + (diff(mc(k).JS_Yaw_Raw)./dt).^2);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err = [stats.th_err];
    
    [jerk,S] = sort(jerk);
    th_err = th_err(S);
    plot(jerk,th_err,'.','Color',colr(n,:),'MarkerSize',8);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(jerk,th_err,2);
        f = polyval(p,jerk);
        [~,I] = sort(jerk);
        hold on;plot(jerk(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(jerk,th_err);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
    
end
grid on;xlabel('jerkiness');ylabel('error (deg)');ylim([-30 30]);xlim([0 10^8]);
title('Angular error as a function of total jerk');legend(leg_input);
% abolute angular error
% subplot(2,2,4);
% for n = 1:length(stnames)-1
%     stats = [];
%     th_err = [];
%     prs = [];
%     taus = [];
%     prs = [trials(stimtype.(stnames{n})).prs];
%     taus = [prs.tau];
%     stats = [trials(stimtype.(stnames{n})).stats];
%     th_err = abs([stats.th_err]);
%     [x,y,b,c] = find_regress(taus,th_err);
%     plot(taus,th_err,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
% end
% grid on;xlabel('\tau (s)');ylabel('ABS error (deg)');ylim([-30 30]);
% title('ABSOLUTE angular error as a function \tau');legend(leg_input);
% 
suptitle(trials(1).prs.subject);
