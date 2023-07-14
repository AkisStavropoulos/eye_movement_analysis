function [rho,pval] = SteeringVsTrackingError(tracking,refplane)

%% Correlation between tracking error and Steering error
 
if strcmp(tracking{1}.misc.sampleflag,'time')
    xaxislim = [0 10];  xaxislabel = 'time [s]';
elseif strcmp(tracking{1}.misc.sampleflag,'distperc')
    xaxislim = [0 1];  xaxislabel = 'distance %';
elseif strcmp(tracking{1}.misc.sampleflag,'distance2end')
    xaxislim = [-600 0];  xaxislabel = 'distance to end [cm]';
end

% extract minimum length of data series
if strcmp(tracking{1}.misc.sampleflag,'time')
    minlength = 8*60;
else
    minlength = min(cellfun(@(x) length(x.misc.ts), tracking(:)));
end

% calculate correlation
rho = cellfun(@(x) x.eyepos.errcorr.targ.(refplane).r(1:minlength)',tracking,'uniformoutput',false);
pval = cellfun(@(x) x.eyepos.errcorr.targ.(refplane).p(1:minlength)',tracking,'uniformoutput',false);

Nsubs = size(rho,1);
Nstim = size(rho,2);

rho_plot = mat2cell([rho{:}], minlength, Nsubs*(ones(Nstim,1)));
pval_plot = mat2cell([pval{:}], minlength, Nsubs*(ones(Nstim,1)));
ts = tracking{1}.misc.ts(1:minlength);

% Plot
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

figure('name',['Correlation between TTE and SE(' refplane ')'],'numbertitle','off','position',[0 685 960 310]);
for s = 1:Nstim
    subplot(1,Nstim,s); hold on;
    plot(ts,rho_plot{s},'color',colr(s,:));
    shadedErrorBar(ts,nanmean(rho_plot{s},2),nanstd(rho_plot{s},[],2)./sqrt(Nsubs),'lineprops',{'color',colr(s,:),'linewidth',3});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim -1 1]); hline(0,'k--');
end
