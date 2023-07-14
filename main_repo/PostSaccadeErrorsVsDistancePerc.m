function [rho,pval,beta] = PostSaccadeErrorsVsDistancePerc(saccade)

%% Tracking errors after saccade as a function of distance percentage
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'ver','hor','all'};

figure('name','Tracking error as a function of distance percentage','numbertitle','off');
for k = 1:numel(fldnames)
for s = 1:Nstim
    % from target
    beta.targ(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.targ.(fldnames{k}).beta,saccade(:,s),'un',0));
    rho.targ(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.targ.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.targ(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.targ.(fldnames{k}).pval,saccade(:,s),'un',0));
    % from stop position
    beta.stop(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.stop.(fldnames{k}).beta,saccade(:,s),'un',0));
    rho.stop(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.stop.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.stop(:,s) = cell2mat(cellfun(@(x) x.error.postsac.distperc.stop.(fldnames{k}).pval,saccade(:,s),'un',0));
    
    subplot(2,numel(fldnames),k); hold on;
    plot(s-0.2*rand(1,Nsubs),rho.targ(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:));
    plot(s+0.2*rand(1,Nsubs),rho.stop(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:)*0.5);
    title(fldnames{k}); ylabel('corr. coefficient'); axis([0 4 -1 1]); xticks([1 2 3]); xticklabels({'vestibular','visual','combined'});
    hline(0,'k');
    
    subplot(2,numel(fldnames),k+numel(fldnames)); hold on;
    plot(s-0.2*rand(1,Nsubs),beta.targ(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:));
    plot(s+0.2*rand(1,Nsubs),beta.stop(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:)*0.5);
    title(fldnames{k}); ylabel('regression slope'); axis([0 4 -40 40]);xticks([1 2 3]); xticklabels({'vestibular','visual','combined'});
    hline(0,'k');
end
end
