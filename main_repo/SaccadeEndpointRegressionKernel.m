function [rho,pval,beta,multi] = SaccadeEndpointRegressionKernel(saccade)

%% Regression KERNEL of saccade endpoint vs target/stop positions early and late in trial
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'hor','ver'};

hor = figure('name','HOR: Endpoint Regression kernel (target vs stop position) over distance %','numbertitle','off');
ver = figure('name','VER: Endpoint Regression kernel (target vs stop position) over distance %','numbertitle','off');
for s = 1:Nstim
for k = 1:numel(fldnames)
    
    ts = saccade{1}.misc.ts;
    % from target
    rho.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.targ.(fldnames{k}).rho,saccade(:,s),'un',0);
    rho.targ.(fldnames{k}){s} = cat(3,rho.targ.(fldnames{k}){s}{:});
    pval.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.targ.(fldnames{k}).pval,saccade(:,s),'un',0);
    pval.targ.(fldnames{k}){s} = cat(3,pval.targ.(fldnames{k}){s}{:});
    beta.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.targ.(fldnames{k}).beta,saccade(:,s),'un',0);
    beta.targ.(fldnames{k}){s} = cat(3,beta.targ.(fldnames{k}){s}{:});
    multi.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.targ.(fldnames{k}).multi,saccade(:,s),'un',0);
    multi.targ.(fldnames{k}){s} = cat(3,multi.targ.(fldnames{k}){s}{:});
    
    % from stop position
    rho.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.stop.(fldnames{k}).rho,saccade(:,s),'un',0);
    rho.stop.(fldnames{k}){s} = cat(3,rho.stop.(fldnames{k}){s}{:});
    pval.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.stop.(fldnames{k}).pval,saccade(:,s),'un',0);
    pval.stop.(fldnames{k}){s} = cat(3,pval.stop.(fldnames{k}){s}{:});
    beta.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.stop.(fldnames{k}).beta,saccade(:,s),'un',0);
    beta.stop.(fldnames{k}){s} = cat(3,beta.stop.(fldnames{k}){s}{:});
    multi.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.kernel_vs_distperc.sacend.stop.(fldnames{k}).multi,saccade(:,s),'un',0);
    multi.stop.(fldnames{k}){s} = cat(3,multi.stop.(fldnames{k}){s}{:});
    
end
ngroups = size(rho.targ.ver{s},1);
for j = 1:ngroups
    figure(hor);
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts,mean(multi.targ.hor{s}(j,:,:),3),std(multi.targ.hor{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 1.5]); hline(0,'k'); vline(0,'k'); 
    title('HOR'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend('target position');
    subplot(2,Nstim,s+Nstim); hold on;
    shadedErrorBar(ts,mean(multi.stop.hor{s}(j,:,:),3),std(multi.stop.hor{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 1.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend('stop position');
    
    figure(ver);
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts,mean(multi.targ.ver{s}(j,:,:),3),std(multi.targ.ver{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 1.5]); hline(0,'k'); vline(0,'k'); 
    title('VER'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend('target position');
    subplot(2,Nstim,s+Nstim); hold on;
    shadedErrorBar(ts,mean(multi.stop.ver{s}(j,:,:),3),std(multi.stop.ver{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 1.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend('stop position');

end
end
