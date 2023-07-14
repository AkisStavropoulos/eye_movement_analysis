function [rho,pval,beta,multi] = SaccadeAmplitudeRegression(saccade)

%% Regress saccade amplitude against TTE and SPTE before saccade and plot weights over distance %
% play a bit here, seems promising
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'hor','ver','all'};

figure('name','Saccade amplitude vs TTE/SPTE over distance %','numbertitle','off');
for k = 1:numel(fldnames)
for s = 1:Nstim
    
    trl_prog = saccade{1}.eyepos.expsac_vs_distperc.groups.val;
    % from target
    rho.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.targ.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.targ.(fldnames{k}).pval,saccade(:,s),'un',0));
    beta.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.targ.(fldnames{k}).beta,saccade(:,s),'un',0));
    multi.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.targ.(fldnames{k}).multi,saccade(:,s),'un',0));
    
    % from stop position
    rho.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.stop.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.stop.(fldnames{k}).pval,saccade(:,s),'un',0));
    beta.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.stop.(fldnames{k}).beta,saccade(:,s),'un',0));
    multi.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.stop.(fldnames{k}).multi,saccade(:,s),'un',0));

    % from belief position
    rho.blv.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.blv.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.blv.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.blv.(fldnames{k}).pval,saccade(:,s),'un',0));
    beta.blv.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.blv.(fldnames{k}).beta,saccade(:,s),'un',0));
    multi.blv.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.expsac_vs_distperc.blv.(fldnames{k}).multi,saccade(:,s),'un',0));

    ngroups = size(rho.stop.(fldnames{k}){s},2);
    
    selection = multi;
    subplot(3,Nstim,s+(k-1)*Nstim); hold on;
    bar(-0.75,mean(selection.targ.(fldnames{k}){s}(:,1)),'facecolor',colr(s,:));errorbar(-0.75,mean(selection.targ.(fldnames{k}){s}(:,1)),std(selection.targ.(fldnames{k}){s}(:,1))./sqrt(Nsubs),'k','capsize',0);
    bar(0,mean(selection.stop.(fldnames{k}){s}(:,1)),'facecolor',colr(s,:)*0.5);errorbar(0,mean(selection.stop.(fldnames{k}){s}(:,1)),std(selection.stop.(fldnames{k}){s}(:,1))./sqrt(Nsubs),'k','capsize',0);
    bar(+0.75,mean(selection.blv.(fldnames{k}){s}(:,1)),'facecolor',[.5 .5 .5]);errorbar(+0.75,mean(selection.blv.(fldnames{k}){s}(:,1)),std(selection.blv.(fldnames{k}){s}(:,1))./sqrt(Nsubs),'k','capsize',0);
    
    bar(3-0.75,mean(selection.targ.(fldnames{k}){s}(:,2)),'facecolor',colr(s,:));errorbar(3-0.75,mean(selection.targ.(fldnames{k}){s}(:,2)),std(selection.targ.(fldnames{k}){s}(:,2))./sqrt(Nsubs),'k','capsize',0);
    bar(3,mean(selection.stop.(fldnames{k}){s}(:,2)),'facecolor',colr(s,:)*0.5);errorbar(3,mean(selection.stop.(fldnames{k}){s}(:,2)),std(selection.stop.(fldnames{k}){s}(:,2))./sqrt(Nsubs),'k','capsize',0);
    bar(3+0.75,mean(selection.blv.(fldnames{k}){s}(:,2)),'facecolor',[.5 .5 .5]);errorbar(3+0.75,mean(selection.blv.(fldnames{k}){s}(:,2)),std(selection.blv.(fldnames{k}){s}(:,2))./sqrt(Nsubs),'k','capsize',0);
    xticks(3*(0:ngroups-1)); xticklabels({'0-35%','35-70%'}); xlabel('distance traveled'); ylabel('regr. weights'); title(fldnames{k});

%     
%     subplot(3,Nstim,s+(k-1)*Nstim); 
%     shadedErrorBar(1:ngroups,mean(multi.targ.(fldnames{k}){s}),std(multi.targ.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); hold on;
%     shadedErrorBar(1:ngroups,mean(multi.stop.(fldnames{k}){s}),std(multi.stop.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5}); 
%     shadedErrorBar(1:ngroups,mean(beta.blv.(fldnames{k}){s}),std(beta.blv.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',[.5 .5 .5]}); 
%      hline(0,'k'); title(fldnames{k}); xlabel('distance %'); ylabel('regr. weights');
end

end
