function [bias,variab] = TrackersVsNonTrackersPerformance(subject,subindx)

params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
[bias,skipsub] = get_bias(subject,params,0);

Nsubs = numel(subject);
Nstim = numel(subindx);
colr = brewermap(Nstim,'Dark2');

%% difference in performance between vertical trackers and non-trackers
% bias
figure('name','Bias differences between trackers and non-trackers','numbertitle','off','position',[0 545 950 450]);
for s = 1:Nstim
    trackers = subindx{s};
    nontrack = ~trackers;
    
    subplot(2,Nstim,s); hold on;
    plot(ones(sum(trackers),1)+0.1*randn(sum(trackers),1),bias.r(trackers,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')
    plot(2*ones(sum(nontrack),1)+0.1*randn(sum(nontrack),1),bias.r(nontrack,s),'o','markeredgecolor',colr(s,:),'markerfacecolor','none')
    ylabel('linear bias'); xticks([1 2]); xticklabels({'trackers','non-trackers'}); axis([0 3 0 2]);hline(1,'k--');
    
    subplot(2,Nstim,s+Nstim); hold on;
    plot(ones(sum(trackers),1)+0.1*randn(sum(trackers),1),bias.th(trackers,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')
    plot(2*ones(sum(nontrack),1)+0.1*randn(sum(nontrack),1),bias.th(nontrack,s),'o','markeredgecolor',colr(s,:),'markerfacecolor','none')
    ylabel('angular bias'); xticks([1 2]); xticklabels({'trackers','non-trackers'}); axis([0 3 0 2]);hline(1,'k--');

end
% variability
for i = 1:Nsubs
    for s = 1:Nstim
        indx = poolindx{i,s};
        r_tar = arrayfun(@(x) x.prs.r_tar, subject(i).trials(indx));      th_tar = arrayfun(@(x) x.prs.th_tar, subject(i).trials(indx));
        r_sub = arrayfun(@(x) x.prs.r_sub, subject(i).trials(indx));      th_sub = arrayfun(@(x) x.prs.th_sub, subject(i).trials(indx));
        r_var(i,s) = std(r_sub - bias.r(i,s)*r_tar);    th_var(i,s) = std(th_sub - bias.th(i,s)*th_tar);
    end
end
figure('name','Variability differences between trackers and non-trackers','numbertitle','off','position',[0 25 950 450]);
for s = 1:Nstim
    trackers = subindx{s};
    nontrack = ~trackers;
    
    subplot(2,Nstim,s); hold on;
    plot(ones(sum(trackers),1)+0.1*randn(sum(trackers),1),r_var(trackers,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')
    plot(2*ones(sum(nontrack),1)+0.1*randn(sum(nontrack),1),r_var(nontrack,s),'o','markeredgecolor',colr(s,:),'markerfacecolor','none')
    ylabel('linear variability'); xticks([1 2]); xticklabels({'trackers','non-trackers'}); axis([0 3 0 200]);
    
    subplot(2,Nstim,s+Nstim); hold on;
    plot(ones(sum(trackers),1)+0.1*randn(sum(trackers),1),th_var(trackers,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')
    plot(2*ones(sum(nontrack),1)+0.1*randn(sum(nontrack),1),th_var(nontrack,s),'o','markeredgecolor',colr(s,:),'markerfacecolor','none')
    ylabel('angular variability'); xticks([1 2]); xticklabels({'trackers','non-trackers'}); axis([0 3 0 15]);

end

variab.r = r_var;
variab.th = th_var;