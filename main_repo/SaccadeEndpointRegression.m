function [rho,pval,beta,multi] = SaccadeEndpointRegression(saccade)

%% Regress endpoint of saccade against target/stop positions before saccade and plot weights over distance %
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'hor','ver'};

figure('name','Saccade endpoint vs Target/Stop location: Over distance %','numbertitle','off');
for k = 1:numel(fldnames)
for s = 1:Nstim
    
    trl_prog = saccade{1}.eyepos.postsac_vs_distperc.groups.val;
    % from target
    rho.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.targ.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.targ.(fldnames{k}).pval,saccade(:,s),'un',0));
    beta.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.targ.(fldnames{k}).beta,saccade(:,s),'un',0));
    multi.targ.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.targ.(fldnames{k}).multi,saccade(:,s),'un',0));
    
    % from stop position
    rho.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.stop.(fldnames{k}).rho,saccade(:,s),'un',0));
    pval.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.stop.(fldnames{k}).pval,saccade(:,s),'un',0));
    beta.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.stop.(fldnames{k}).beta,saccade(:,s),'un',0));
    multi.stop.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.eyepos.postsac_vs_distperc.stop.(fldnames{k}).multi,saccade(:,s),'un',0));

    subplot(2,Nstim,s+(k-1)*Nstim); 
    shadedErrorBar(trl_prog,median(multi.targ.(fldnames{k}){s}),std(multi.targ.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); hold on;
    shadedErrorBar(trl_prog,median(multi.stop.(fldnames{k}){s}),std(multi.stop.(fldnames{k}){s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5}); 
    axis([0 1 -1 2]); hline(0,'k'); title(fldnames{k}); xlabel('distance %'); ylabel('regr. weights');


end
end
