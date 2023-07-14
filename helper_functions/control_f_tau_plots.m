function control_f_tau_plots(T,s,T_trial,s_trial,condition,trials)

%% Control differences as a function of tau

for i = 1:length(trials)
taus(i) = trials(i).prs.tau;
end
condnames = fieldnames(condition);
if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2')) % stimtype
    cond = [{'vestibular'} {'visual'} {'combined'} {'all trials'}];
end
for n = 1:length(condnames)*3
    h(n) = figure;
end
for k = 1:length(condnames)
    for i = 1:length(condition.(condnames{k}))
        indx = condition.(condnames{k})(i);
        colr = parula(length(condition.(condnames{k})));
        figure(h(k));
        hold on;
        plot(taus(indx),s_trial(indx) - s(indx),'.','Color',colr(i,:)); % switchtime
        hline(0);title([cond{k} ' - switchtime difference as a function of \tau']);xlabel('\tau (s)');ylabel('switchtime difference (s)');
        g = colorbar('Ticks',linspace(0,1,5),...
            'TickLabels',linspace(0,length(condition.(condnames{k})),5));
        ylabel(g,'Experience in trials')
        
        figure(h(k+length(condnames)));
        hold on;
        plot(taus(indx),T_trial(indx) - T(indx),'.','Color',colr(i,:)); % Trial duration
        hline(0);title([cond{k} ' - Trial duration difference as a function of \tau']);xlabel('\tau (s)');ylabel('Trial duration difference (s)');
        g = colorbar('Ticks',linspace(0,1,5),...
            'TickLabels',linspace(0,length(condition.(condnames{k})),5));
        ylabel(g,'Experience in trials')
        
        figure(h(k+2*length(condnames)));
        hold on;
        plot(taus(indx),T_trial(indx) - s_trial(indx) - (T(indx)-s(indx)),'.','Color',colr(i,:)); % brake time
        hline(0);title([cond{k} ' - Brake time difference as a function of \tau']);xlabel('\tau (s)');ylabel('Brake time difference (s)');
        g = colorbar('Ticks',linspace(0,1,5),...
            'TickLabels',linspace(0,length(condition.(condnames{k})),5));
        ylabel(g,'Experience in trials');
        
    end
end
for n = 1:length(h)
    figure(h(n));suptitle(trials(1).prs.subject)
end