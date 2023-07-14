function [rho,ts] = EyePositionCorr(tracking)

%% Position gain between eye and target/stop positions
 
Nsubs = size(tracking,1);
Nstim = size(tracking,2);

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

ts = tracking{1}.misc.ts(1:minlength);

for i = 1:Nsubs
    for s = 1:Nstim
        % horizontal
        rho.hor.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.hor_mean.rho.startaligned(1:minlength);
        rho.hor.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.hor_mean.rho.startaligned(1:minlength);
        % vertical
        rho.ver.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.ver_mean.rho.startaligned(1:minlength);
        rho.ver.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.ver_mean.rho.startaligned(1:minlength);
    end
end

%% All subjects
colr = brewermap(Nstim,'Dark2');
figure('name','Correlation coefficient','numbertitle','off');
for s = 1:Nstim
    % horizontal
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts(1:minlength),mean(rho.hor.targ{s}),std(rho.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(rho.hor.stop{s}),std(rho.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    legend('target position','stop position');
    title('horizontal');
    % vertical
    subplot(2,Nstim,s+Nstim); hold on; 
    shadedErrorBar(ts(1:minlength),mean(rho.ver.targ{s}),std(rho.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(rho.ver.stop{s}),std(rho.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    legend('target position','stop position');
    title('vertical');
end

%% Subjects separately
if 0
hor = figure('name','Horizontal Correlations','numbertitle','off'); ver = figure('name','Vertical Correlations','numbertitle','off');
for s = 1:Nstim
    % horizontal
    figure(hor);
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),rho.hor.targ{s},'color',colr(s,:));
    shadedErrorBar(ts(1:minlength),mean(rho.hor.targ{s}),std(rho.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    title('target position');
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),rho.hor.stop{s},'color',colr(s,:)*0.5);
    shadedErrorBar(ts(1:minlength),mean(rho.hor.stop{s}),std(rho.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    title('stop position');
    % vertical
    figure(ver);
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),rho.ver.targ{s},'color',colr(s,:));
    shadedErrorBar(ts(1:minlength),mean(rho.ver.targ{s}),std(rho.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    title('target position');
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),rho.ver.stop{s},'color',colr(s,:)*0.5);
    shadedErrorBar(ts(1:minlength),mean(rho.ver.stop{s}),std(rho.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('correlation coef.'); axis([xaxislim 0 1]);
    title('stop position');
end
end
