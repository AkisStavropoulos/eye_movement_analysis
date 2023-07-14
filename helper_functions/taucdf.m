function taucdf(subject,distr_flag)
%% Distribution of taus for each subject
% distr_flag: choose if you want to plot CDF or PDF
% distr_flag: 'pdf' or 'cdf'

if strcmp(distr_flag,'pdf')
    plt = 1;
elseif strcmp(distr_flag,'cdf')
    plt = 0;
else
    error('2nd input must be a string, ''pdf'' or ''cdf''. ');
end

figure;hold on;
colr = brewermap(length(subject),'Paired');
for i = 1:length(subject)
    
    for j = 1:length(subject(i).trials)
        tau_temp(j) = subject(i).trials(j).prs.tau;
    end
    if plt
        figure;hist(tau_temp,50);xlabel('\tau [sec]');title([subject(i).name ' - \tau distribution']);
        xlim([0 8]);
    else
    taus{i} = tau_temp;
    [f,x] = ecdf(taus{i});
    plot(x,f,'Color',colr(i,:));xlabel('\tau [sec]');title('CDF of \tau distributions for all subjects');
    leg_input{i} = subject(i).name;
    end
end
legend(leg_input,'Location','southeast');
