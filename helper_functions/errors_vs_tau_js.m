function errors_vs_tau_js(trials,stimtype)
stnames = fieldnames(stimtype);
%% see distance error as a function of tau
figure;subplot(2,2,1);
colr = colormap(jet);
colr = colr(round(linspace(1,length(colr),length(stnames))),:);
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    stats = [];
    d_err2 = [];
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    d_err2 = [stats.d_err2];
    [x,y,b,c] = find_regress(taus,d_err2);
    plot(taus,d_err2,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('\tau (s)');ylabel('error (cm)');
title('Distance error as a function \tau');legend(leg_input);
% absolute distance error
subplot(2,2,2);
for n = 1:length(stnames)-1
    stats = [];
    d_err2 = [];
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    d_err2 = abs([stats.d_err2]);
    [x,y,b,c] = find_regress(taus,d_err2);
    plot(taus,d_err2,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('\tau (s)');ylabel('ABS error (cm)');
title('ABSOLUTE distance error as a function \tau');legend(leg_input);

suptitle(trials(1).prs.subject);
%% see angular error as a function of tau
subplot(2,2,3);
colr = colormap(jet);
colr = colr(round(linspace(1,length(colr),length(stnames))),:);
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    stats = [];
    th_err2 = [];
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err2 = [stats.th_err2];
    [x,y,b,c] = find_regress(taus,th_err2);
    plot(taus,th_err2,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('\tau (s)');ylabel('error (deg)');ylim([-30 30]);
title('Angular error as a function \tau');legend(leg_input);
% abolute angular error
subplot(2,2,4);
for n = 1:length(stnames)-1
    stats = [];
    th_err2 = [];
    prs = [];
    taus = [];
    prs = [trials(stimtype.(stnames{n})).prs];
    taus = [prs.tau];
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err2 = abs([stats.th_err2]);
    [x,y,b,c] = find_regress(taus,th_err2);
    plot(taus,th_err2,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('\tau (s)');ylabel('ABS error (deg)');ylim([-30 30]);
title('ABSOLUTE angular error as a function \tau');legend(leg_input);

suptitle(trials(1).prs.subject);
