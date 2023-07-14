function eucl_errors_vs_jerk(trials,stimtype,regr)
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
    eucl_err = [];
    mc = [];
    jerk = [];
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        jerk(k) = trials(k).continuous.cumjerk(end);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    eucl_err = [stats.eucl_err];
    [jerk,S] = sort(jerk);
    eucl_err = eucl_err(S);
    
    plot(eucl_err,jerk,'.','Color',colr(n,:),'MarkerSize',8);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(eucl_err,jerk,2);
        f = polyval(p,eucl_err);
        [~,I] = sort(eucl_err);
        hold on;plot(eucl_err(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(eucl_err,jerk);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
end
grid on;ylabel('jerkiness');xlabel('Euclideian error [cm]');xlim([0 500]);ylim([0 100*10^2]);
title('Euclideian error as a function of total jerk');legend(leg_input);

%% Compute acceleration
for n = 1:length(stnames)-1
    stats = [];
    eucl_err = [];
    mc = [];
    acc = [];
    mc = [trials(stimtype.(stnames{n})).mc];
    for k = 1:length(mc)
        acc(k) = trials(k).continuous.cumacc(end);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    eucl_err = [stats.eucl_err];
    [acc,S] = sort(acc);
    eucl_err = eucl_err(S);
    
    subplot(2,1,2);
    plot(eucl_err,acc,'.','Color',colr(n,:),'MarkerSize',8);%xlim([0 8]);
    
    if strcmp(regr,'poly')
        p = polyfit(eucl_err,acc,2);
        f = polyval(p,eucl_err);
        [~,I] = sort(eucl_err);
        hold on;plot(eucl_err(I),f(I),'Color',colr(n,:),'LineWidth',2.5);
    elseif strcmp(regr,'lin')
        [x,y,b,c] = find_regress(eucl_err,acc);
        hold on;plot(x,b*x + c,'Color',colr(n,:),'LineWidth',2.5);
    end
end
grid on;ylabel('acceleration');xlabel('Euclideian error [cm]');xlim([0 500]);ylim([0 1000]);
title('Euclideian error as a function of total acceleration');legend(leg_input);



suptitle(trials(1).prs.subject);
