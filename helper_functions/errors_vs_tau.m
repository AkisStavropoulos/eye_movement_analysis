function errors_vs_tau(trials,stimtype,regr)
% regr: type of regression, 'poly' for 2nd order polynomial, 'lin' for linear regression
    
stnames = fieldnames(stimtype);
%% see distance error as a function of tau
figure;subplot(2,1,1);
colr = brewermap(length(stnames),'Dark2');
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    stats = [];
    d_err = [];
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    d_err = [stats.d_err];
    [taus,S] = sort(taus);
    d_err = d_err(S);
    
    plot(taus,d_err,'.','Color',colr(n,:),'MarkerSize',8);xlim([0 8]);ylim([-500 500]);
    
    if strcmp(regr,'poly')
        p = polyfit(taus,d_err,2);
        f = polyval(p,taus);
        [~,I] = sort(taus);
        hold on;plot(taus(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(taus,d_err);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
end
grid on;xlabel('\tau (s)');ylabel('error (cm)');
title('Distance error as a function \tau');legend(leg_input);
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
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err = [stats.th_err];
    [taus,S] = sort(taus);
    th_err = th_err(S);
    plot(taus,th_err,'.','Color',colr(n,:),'MarkerSize',8);xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(taus,th_err,2);
        f = polyval(p,taus);
        [~,I] = sort(taus);
        hold on;plot(taus(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(taus,th_err);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
    
end
grid on;xlabel('\tau (s)');ylabel('error (deg)');ylim([-30 30]);
title('Angular error as a function \tau');legend(leg_input);
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
