function [rho,pval,beta,multi] = SaccadeAmplitudeKernel(saccade,type_regr,quantity,split_grp)

%% Regression KERNEL of saccade amplitude vs tracking errors early and late in trial
% type_regr: 'linear', 'robust', or 'ridge'
% quantity: specify quantity to plot, between 'rho', 'beta', 'multi' for
%           correlation, simple linear regression, or multiple linear rergession, respectively.
% split_grp: split in groups of early and late trial

if ~any(strcmpi(type_regr,{'linear','robust','ridge'}))
    error('Specify type of regression in 2nd output ("linear","robust","ridge")');
end

Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'hor','ver'};

if split_grp
hor = figure('name',['HOR: ' type_regr ' Regression kernel (TTE / SPTE / BTE) over distance %'],'numbertitle','off');
ver = figure('name',['VER: ' type_regr ' Regression kernel (TTE / SPTE / BTE) over distance %'],'numbertitle','off');
leg_inp_tte = {'early TTE','late TTE'};
leg_inp_spte = {'early SPTE','late SPTE'};
leg_inp_bte = {'early BTE','late BTE'};

var2plt = 'kernel_vs_distperc';
else
hor = figure('name',[type_regr ' Regression kernel (TTE / SPTE / BTE)'],'numbertitle','off');
leg_inp = {'TTE','SPTE','BTE'}; 

var2plt = 'kernel';
end

for s = 1:Nstim
for k = 1:numel(fldnames)
    
    ts = saccade{1}.misc.ts;
    % from target
    rho.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.targ.(fldnames{k}).rho,saccade(:,s),'un',0);
    rho.targ.(fldnames{k}){s} = cat(3,rho.targ.(fldnames{k}){s}{:});
    pval.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.targ.(fldnames{k}).pval,saccade(:,s),'un',0);
    pval.targ.(fldnames{k}){s} = cat(3,pval.targ.(fldnames{k}){s}{:});
    beta.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.targ.(fldnames{k}).beta.(type_regr),saccade(:,s),'un',0);
    beta.targ.(fldnames{k}){s} = cat(3,beta.targ.(fldnames{k}){s}{:});
    multi.targ.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.targ.(fldnames{k}).multi2.(type_regr),saccade(:,s),'un',0);
    multi.targ.(fldnames{k}){s} = cat(3,multi.targ.(fldnames{k}){s}{:});
    
    % from stop position
    rho.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.stop.(fldnames{k}).rho,saccade(:,s),'un',0);
    rho.stop.(fldnames{k}){s} = cat(3,rho.stop.(fldnames{k}){s}{:});
    pval.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.stop.(fldnames{k}).pval,saccade(:,s),'un',0);
    pval.stop.(fldnames{k}){s} = cat(3,pval.stop.(fldnames{k}){s}{:});
    beta.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.stop.(fldnames{k}).beta.(type_regr),saccade(:,s),'un',0);
    beta.stop.(fldnames{k}){s} = cat(3,beta.stop.(fldnames{k}){s}{:});
    multi.stop.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.stop.(fldnames{k}).multi2.(type_regr),saccade(:,s),'un',0);
    multi.stop.(fldnames{k}){s} = cat(3,multi.stop.(fldnames{k}){s}{:});
    
    % from belief position
    rho.blv.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.blv.(fldnames{k}).rho,saccade(:,s),'un',0);
    rho.blv.(fldnames{k}){s} = cat(3,rho.blv.(fldnames{k}){s}{:});
    pval.blv.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.blv.(fldnames{k}).pval,saccade(:,s),'un',0);
    pval.blv.(fldnames{k}){s} = cat(3,pval.blv.(fldnames{k}){s}{:});
    beta.blv.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.blv.(fldnames{k}).beta.(type_regr),saccade(:,s),'un',0);
    beta.blv.(fldnames{k}){s} = cat(3,beta.blv.(fldnames{k}){s}{:});
    multi.blv.(fldnames{k}){s} = cellfun(@(x) x.eyepos.(var2plt).sacamp.blv.(fldnames{k}).multi2.(type_regr),saccade(:,s),'un',0);
    multi.blv.(fldnames{k}){s} = cat(3,multi.blv.(fldnames{k}){s}{:});

end

ngroups = size(rho.targ.ver{s},1);

if strcmpi(quantity,'rho');         temp = rho; 
elseif strcmpi(quantity,'beta');      temp = beta; 
elseif strcmpi(quantity,'multi');       temp = multi; %temp.blv = beta.blv;
end

if split_grp
    
for j = 1:ngroups
    figure(hor);
    subplot(3,Nstim,s); hold on;
    shadedErrorBar(ts,mean(temp.targ.hor{s}(j,:,:),3),std(temp.targ.hor{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); 
    title('HOR'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_tte{:});
    subplot(3,Nstim,s+Nstim); hold on;
    shadedErrorBar(ts,mean(temp.stop.hor{s}(j,:,:),3),std(temp.stop.hor{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_spte{:});
    subplot(3,Nstim,s+2*Nstim); hold on;
    shadedErrorBar(ts,mean(temp.blv.hor{s}(j,:,:),3),std(temp.blv.hor{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',[.5 .5 .5]*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_bte{:});
    
    figure(ver);
    subplot(3,Nstim,s); hold on;
    shadedErrorBar(ts,mean(temp.targ.ver{s}(j,:,:),3),std(temp.targ.ver{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); 
    title('VER'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_tte{:});
    subplot(3,Nstim,s+Nstim); hold on;
    shadedErrorBar(ts,mean(temp.stop.ver{s}(j,:,:),3),std(temp.stop.ver{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_spte{:});
    subplot(3,Nstim,s+2*Nstim); hold on;
    shadedErrorBar(ts,mean(temp.blv.ver{s}(j,:,:),3),std(temp.blv.ver{s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',[.5 .5 .5]*(1/j)});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp_bte{:});

end

else
    figure(hor);
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts,mean(temp.targ.hor{s},3),std(temp.targ.hor{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)});
    shadedErrorBar(ts,mean(temp.stop.hor{s},3),std(temp.stop.hor{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*0.5});
    shadedErrorBar(ts,mean(temp.blv.hor{s},3),std(temp.blv.hor{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',[.5 .5 .5]});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k');
    title('HOR'); xlabel('time since saccade onset [s]'); ylabel('regr. weights'); 
    subplot(2,Nstim,s+Nstim); hold on;
    shadedErrorBar(ts,mean(temp.targ.ver{s},3),std(temp.targ.ver{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)});
    shadedErrorBar(ts,mean(temp.stop.ver{s},3),std(temp.stop.ver{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*0.5});
    shadedErrorBar(ts,mean(temp.blv.ver{s},3),std(temp.blv.ver{s},[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',[.5 .5 .5]});
    axis([-0.8 0.8 -0.5 0.5]); hline(0,'k'); vline(0,'k');
    title('VER'); xlabel('time since saccade onset [s]'); ylabel('regr. weights');
    xlabel('time since saccade onset [s]'); ylabel('regr. weights'); legend(leg_inp{:});
    
    
end
end
suptitle(quantity);
