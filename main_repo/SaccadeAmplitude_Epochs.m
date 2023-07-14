function SaccadeAmplitude_Epochs(saccade)

%% Saccade amplitude for different trial epochs (Target ON, Steering, End of Trial)
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

figure('name','saccade amplitude for trial epochs','numbertitle','off','position',[65 720 840 260]);
for s = 1:Nstim
    subplot(1,Nstim,s); hold on;
    % target on
    cprob = cell2mat(cellfun(@(x) x.epochs.distperc.targeton.amplitude.cdf,saccade(:,s),'un',0));
    nx = unique(cell2mat(cellfun(@(x) x.epochs.distperc.targeton.amplitude.nx,saccade(:,s),'un',0)));
    shadedErrorBar(nx,mean(cprob),std(cprob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0});
    % steering
    cprob = cell2mat(cellfun(@(x) x.epochs.distperc.steering.amplitude.cdf,saccade(:,s),'un',0));
    nx = unique(cell2mat(cellfun(@(x) x.epochs.distperc.steering.amplitude.nx,saccade(:,s),'un',0)));
    shadedErrorBar(nx,mean(cprob),std(cprob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    % end of trial
    cprob = cell2mat(cellfun(@(x) x.epochs.distperc.trialend.amplitude.cdf,saccade(:,s),'un',0));
    nx = unique(cell2mat(cellfun(@(x) x.epochs.distperc.trialend.amplitude.nx,saccade(:,s),'un',0)));
    shadedErrorBar(nx,mean(cprob),std(cprob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});

    title('distance %'); ylabel('saccade probability'); xlabel('saccade amplitude [deg]');
    legend('target ON','steering','trial end (>70%)','location','southeast'); xlim([0 40]);
end
