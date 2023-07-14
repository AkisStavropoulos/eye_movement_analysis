function [mu,sd] = ErrorsRatioKernel(saccade)

%% Regression KERNEL of saccade amplitude vs tracking errors early and late in trial
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
fldnames = {'hor','ver','all'};

figure('name','Errors Ratio kernel over distance %','numbertitle','off');

for s = 1:Nstim
    for k = 1:numel(fldnames)
        
        ts = saccade{1}.misc.ts;
        % errors ratio kernel
        mu.(fldnames{k}){s} = cellfun(@(x) x.error.ratio_vs_distperc.(fldnames{k}).kernel.mu,saccade(:,s),'un',0);
        mu.(fldnames{k}){s} = cat(3,mu.(fldnames{k}){s}{:});
        sd.(fldnames{k}){s} = cellfun(@(x) x.error.ratio_vs_distperc.(fldnames{k}).kernel.sd,saccade(:,s),'un',0);
        sd.(fldnames{k}){s} = cat(3,sd.(fldnames{k}){s}{:});
            
        ngroups = size(mu.(fldnames{k}){s},1);
        for j = 1:ngroups
            subplot(3,Nstim,s+(k-1)*Nstim); hold on;
            shadedErrorBar(ts,mean(mu.(fldnames{k}){s}(j,:,:),3),std(mu.(fldnames{k}){s}(j,:,:),[],3)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)*(1/j)});
            xlim([-0.8 0.8 ]); hline(0,'k'); vline(0,'k');
            title((fldnames{k})); xlabel('time since saccade onset [s]'); ylabel('errors ratio (TTE/SPTE)');
            
        end
    end
    
end

% Ungrouped
if 0 
figure;
for s = 1:Nstim
    for k = 1:numel(fldnames)
        
        ts = saccade{1}.misc.ts;
        % errors ratio kernel
        mu.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.error.ratio.(fldnames{k}).kernel.mu,saccade(:,s),'un',0)');
        sd.(fldnames{k}){s} = cell2mat(cellfun(@(x) x.error.ratio.(fldnames{k}).kernel.sd,saccade(:,s),'un',0)');
        
        subplot(3,Nstim,s+(k-1)*Nstim); hold on;
        shadedErrorBar(ts,mean(mu.(fldnames{k}){s},2),std(mu.(fldnames{k}){s},[],2)./sqrt(Nsubs),'lineprops',{'linewidth',2,'color',colr(s,:)});
        xlim([-0.8 0.8 ]); hline(0,'k'); vline(0,'k');
        title((fldnames{k})); xlabel('time since saccade onset [s]'); ylabel('errors ratio (TTE/SPTE)');
        
    end
end    
end


