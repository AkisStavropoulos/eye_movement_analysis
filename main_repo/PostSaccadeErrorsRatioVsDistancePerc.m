function [rho,pval,beta,mu,sd] = PostSaccadeErrorsRatioVsDistancePerc(saccade)

%% Errors Ratio change after saccade as a function of distance percentage
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'ver','hor','all'};

figure('name','Errors Ratio change as a function of distance percentage','numbertitle','off');
for k = 1:numel(fldnames)
for s = 1:Nstim
    % Errors ratio difference before and after saccade
    beta.diff.(fldnames{k})(:,s) = cell2mat(cellfun(@(x) x.error.ratio.(fldnames{k}).diff.beta,saccade(:,s),'un',0));
    rho.diff.(fldnames{k})(:,s) = cell2mat(cellfun(@(x) x.error.ratio.(fldnames{k}).diff.rho,saccade(:,s),'un',0));
    pval.diff.(fldnames{k})(:,s) = cell2mat(cellfun(@(x) x.error.ratio.(fldnames{k}).diff.pval,saccade(:,s),'un',0));
    
    subplot(2,numel(fldnames),k); hold on;
    plot(s-0.2*rand(1,Nsubs),rho.diff.(fldnames{k})(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:));
    plot(s+0.2*rand(1,Nsubs),rho.diff.(fldnames{k})(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:)*0.5);
    title(fldnames{k}); ylabel('corr. coefficient'); axis([0 4 -1 1]); xticks([1 2 3]); xticklabels({'vestibular','visual','combined'});
    hline(0,'k');
    
    subplot(2,numel(fldnames),k+numel(fldnames)); hold on;
    plot(s-0.2*rand(1,Nsubs),beta.diff.(fldnames{k})(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:));
    plot(s+0.2*rand(1,Nsubs),beta.diff.(fldnames{k})(:,s),'o','markeredgecolor','none','markerfacecolor',colr(s,:)*0.5);
    title(fldnames{k}); ylabel('regression slope'); axis([0 4 -40 40]);xticks([1 2 3]); xticklabels({'vestibular','visual','combined'});
    hline(0,'k');
end
end



    
figure('name','Errors Ratio change: Over distance %','numbertitle','off');
for k = 1:numel(fldnames)
for s = 1:Nstim
    
    trl_prog = saccade{1}.error.ratio_vs_distperc.groups.diff.val;
    % error change
    mu.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.error.ratio_vs_distperc.(fldnames{k}).diff.mu,saccade(:,s),'un',0));
    sd.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.error.ratio_vs_distperc.(fldnames{k}).diff.sd,saccade(:,s),'un',0));
    
    subplot(3,Nstim,s+(k-1)*Nstim); 
    shadedErrorBar(trl_prog,mean(mu.(fldnames{k}){s}),std(mu.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); hold on;
    hline(0,'k'); title(fldnames{k}); xlabel('distance %'); ylabel('ratio change [pre-post]');


end
end
