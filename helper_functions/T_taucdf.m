function T_taucdf(subject,distr_flag)
%% distribution of T/tau
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
colr = jet(length(subject));
for i = 1:length(subject)
    T = [];
    tau = [];
    for j = 1:length(subject(i).trials)
        T(j) = subject(i).trials(j).continuous.ts(end);
        tau(j) = subject(i).trials(j).prs.tau;
    end
    
    if plt
        figure;hist(T./tau,50);xlabel('T/\tau');title([subject(i).name ' - T/\tau distribution']);
        xlim([0 25]);
    else
        [f,x] = ecdf(T./tau);
        plot(x,f,'Color',colr(i,:));xlabel('T/\tau');title('CDF of T/\tau distribution for all subjects');
        xlim([0 30]);
        leg_input{i} = subject(i).name;
    end
end
legend(leg_input{:},'Location','southeast');