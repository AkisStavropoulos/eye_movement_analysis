function [multi,ts] = EyePositionMultiRegr(tracking)

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
        multi.hor.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.hor_mean.multi.startaligned(1:minlength);
        multi.hor.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.hor_mean.multi.startaligned(1:minlength);
        % vertical
        multi.ver.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.ver_mean.multi.startaligned(1:minlength);
        multi.ver.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.ver_mean.multi.startaligned(1:minlength);
    end
end

%% All subjects
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

figure('name','Multiple Regression','numbertitle','off');
for s = 1:Nstim
    % horizontal
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts(1:minlength),mean(multi.hor.targ{s}),std(multi.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(multi.hor.stop{s}),std(multi.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('regression coef.'); axis([xaxislim -0.5 1.5]); hline(1,'k--');
    legend('target position','stop position');
    title('horizontal');
    % vertical
    subplot(2,Nstim,s+Nstim); hold on; 
    shadedErrorBar(ts(1:minlength),mean(multi.ver.targ{s}),std(multi.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(multi.ver.stop{s}),std(multi.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('regression coef.'); axis([xaxislim -0.5 1.5]); hline(1,'k--');
    legend('target position','stop position');
    title('vertical');
end

%% Subjects separately
if 0
fig = figure('name','Horizontal Multiple Regression','numbertitle','off'); 
for s = 1:Nstim
    % horizontal
    figure(fig);
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),multi.hor.targ{s},'color',colr(s,:),'handlevisibility','off');
    shadedErrorBar(ts(1:minlength),mean(multi.hor.targ{s}),std(multi.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:),'linewidth',3});
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),multi.hor.stop{s},'color',colr(s,:)*0.5,'handlevisibility','off');
    shadedErrorBar(ts(1:minlength),mean(multi.hor.stop{s}),std(multi.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5,'linewidth',3});
    xlabel(xaxislabel); ylabel('regression coef.'); axis([xaxislim -0.5 2]); hline(1,'k--'); 
    legend('target position','stop position');
    % vertical
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),multi.ver.targ{s},'color',colr(s,:),'handlevisibility','off');
    shadedErrorBar(ts(1:minlength),mean(multi.ver.targ{s}),std(multi.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:),'linewidth',3});
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),multi.ver.stop{s},'color',colr(s,:)*0.5,'handlevisibility','off');
    shadedErrorBar(ts(1:minlength),mean(multi.ver.stop{s}),std(multi.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5,'linewidth',3});
    xlabel(xaxislabel); ylabel('regression coef.'); axis([xaxislim -0.5 2]); hline(1,'k--'); 
    legend('target position','stop position');
end
end

if 0
fig = figure('name','Multiple Regression','numbertitle','off');
for i = 1:Nsubs
    for s = 1:Nstim
    % horizontal
    figure(fig);
    subplot(2,Nstim,s); plot(ts(1:minlength),multi.hor.targ{s}(i,:),'color',colr(s,:)); hold on; 
    subplot(2,Nstim,s); plot(ts(1:minlength),multi.hor.stop{s}(i,:),'color',colr(s,:)*0.5); 
    xlabel(xaxislabel); ylabel('regression coef. (HOR)'); axis([xaxislim 0 2]); hline(1,'k--'); hold off;
    legend('target position','stop position');
    % vertical
    subplot(2,Nstim,s+Nstim); plot(ts(1:minlength),multi.ver.targ{s}(i,:),'color',colr(s,:)); hold on; 
    subplot(2,Nstim,s+Nstim); plot(ts(1:minlength),multi.ver.stop{s}(i,:),'color',colr(s,:)*0.5); 
    xlabel(xaxislabel); ylabel('regression coef. (VER)'); axis([xaxislim 0 2]); hline(1,'k--'); hold off;
    legend('target position','stop position');
    end
end
end

