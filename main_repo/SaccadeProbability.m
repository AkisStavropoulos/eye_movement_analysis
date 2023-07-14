function SaccadeProbability(saccade)

%% Saccade Probability
Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

figure('name','saccade probability over distance %','numbertitle','off','position',[65 720 840 260]);
for s = 1:Nstim
    % time
    subplot(1,3,1); hold on;
    prob = cell2mat(cellfun(@(x) x.probability.time.startaligned.val,saccade(:,s),'un',0));
    ts = unique(cell2mat(cellfun(@(x) x.probability.time.startaligned.t,saccade(:,s),'un',0)));
    shadedErrorBar(ts,mean(prob),std(prob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); ylim([0 0.3]);
    xlabel('time [s]'); ylabel('saccade probability');
    % distance to end
    subplot(1,3,2); hold on;
    prob = cell2mat(cellfun(@(x) x.probability.dist2end.startaligned.val,saccade(:,s),'un',0));
    ts = unique(cell2mat(cellfun(@(x) x.probability.dist2end.startaligned.t,saccade(:,s),'un',0)));
    shadedErrorBar(ts,mean(prob),std(prob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); ylim([0 0.3]);
    xlabel('distance to end [cm]'); ylabel('saccade probability');
    % distance %
    subplot(1,3,3); hold on;
    prob = cell2mat(cellfun(@(x) x.probability.distperc.startaligned.val,saccade(:,s),'un',0));
    ts = unique(cell2mat(cellfun(@(x) x.probability.distperc.startaligned.t,saccade(:,s),'un',0)));
    shadedErrorBar(ts,mean(prob),std(prob)./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}); ylim([0 0.3]);
    xlabel('distance %'); ylabel('saccade probability');
    
end
