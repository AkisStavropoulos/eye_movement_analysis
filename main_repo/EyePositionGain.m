function [g,ts,xaxis] = EyePositionGain(tracking)

%% Position gain between eye and target/stop positions

Nsubs = size(tracking,1);
Nstim = size(tracking,2);

if strcmp(tracking{1}.misc.sampleflag,'time')
    xaxis.lim = [0 10];  xaxis.label = 'time [s]';
elseif strcmp(tracking{1}.misc.sampleflag,'distperc')
    xaxis.lim = [0 1];  xaxis.label = 'distance %';
elseif strcmp(tracking{1}.misc.sampleflag,'distance2end')
    xaxis.lim = [-600 0];  xaxis.label = 'distance to end [cm]';
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
        g.hor.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.hor_mean.beta.startaligned(1:minlength);
        g.hor.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.hor_mean.beta.startaligned(1:minlength);
        % vertical
        g.ver.targ{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.targ.ver_mean.beta.startaligned(1:minlength);
        g.ver.stop{s}(i,:) = tracking{i,s}.eyepos.pred_vs_true.stop.ver_mean.beta.startaligned(1:minlength);
    end
end

%% All subjects
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure('name','Position Gain','numbertitle','off');
for s = 1:Nstim
    % horizontal
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts(1:minlength),mean(g.hor.targ{s}),std(g.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(g.hor.stop{s}),std(g.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    legend('target position','stop position');
    title('horizontal');
    % vertical
    subplot(2,Nstim,s+Nstim); hold on; 
    shadedErrorBar(ts(1:minlength),mean(g.ver.targ{s}),std(g.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts(1:minlength),mean(g.ver.stop{s}),std(g.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    legend('target position','stop position');
    title('vertical');
end

%% Subjects separately
if 0
hor = figure('name','Horizontal Gains','numbertitle','off'); ver = figure('name','Vertical Gains','numbertitle','off');
for s = 1:Nstim
    % horizontal
    figure(hor);
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),g.hor.targ{s},'color',colr(s,:));
    shadedErrorBar(ts(1:minlength),mean(g.hor.targ{s}),std(g.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    title('target position');
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),g.hor.stop{s},'color',colr(s,:)*0.5);
    shadedErrorBar(ts(1:minlength),mean(g.hor.stop{s}),std(g.hor.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    title('stop position');
    % vertical
    figure(ver);
    subplot(2,Nstim,s); hold on; plot(ts(1:minlength),g.ver.targ{s},'color',colr(s,:));
    shadedErrorBar(ts(1:minlength),mean(g.ver.targ{s}),std(g.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    title('target position');
    subplot(2,Nstim,s+Nstim); hold on; plot(ts(1:minlength),g.ver.stop{s},'color',colr(s,:)*0.5);
    shadedErrorBar(ts(1:minlength),mean(g.ver.stop{s}),std(g.ver.stop{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    title('stop position');
end
end

%%
%%
%% Check gain

if 0
    
% vertical
figure;
for i = 1:Nsubs
    for s = 1:Nstim
        for t = 1:minlength
        eye = tracking{i,s}.eyepos.screen.ver_mean(:,t);
        targ = tracking{i,s}.tarpos.screen.ver_mean(:,t);
        stop = tracking{i,s}.stopos.screen.ver_mean(:,t);
        b_targ = regress(eye(:),targ(:));
        b_stop = regress(eye(:),stop(:));

        subplot(2,Nstim,s); hold off;
        plot(targ,eye,'.','color',colr(s,:)); hold on;
        plot(-90:10,-90:10,'k--'); plot(-90:10,(-90:10)*g.ver.targ{s}(i,t),'r');plot(-90:10,(-90:10)*b_targ,'r--');
        xlabel('target position [deg]'); ylabel('eye position [deg]'); axis([-90 10 -90 10]); vline(0,'k-'); hline(0,'k-');
        title(['t = ' num2str(t/60,'%.2f') ' s']);
        
        subplot(2,Nstim,s+Nstim); hold off;
        plot(stop,eye,'.','color',colr(s,:)*0.5); hold on;
        plot(-90:10,-90:10,'k--'); plot(-90:10,(-90:10)*g.ver.stop{s}(i,t),'r');plot(-90:10,(-90:10)*b_stop,'r--');
        xlabel('stop position [deg]'); ylabel('eye position [deg]'); axis([-90 10 -90 10]); vline(0,'k-'); hline(0,'k-');
        
        pause(0.05)
        end
    end
end
% horizontal
figure;
for i = 1:Nsubs
    for s = 1:Nstim
        for t = 1:minlength
        eye = tracking{i,s}.eyepos.screen.hor_mean(:,t);
        targ = tracking{i,s}.tarpos.screen.hor_mean(:,t);
        stop = tracking{i,s}.stopos.screen.hor_mean(:,t);
        b_targ = regress(eye(:),targ(:));
        b_stop = regress(eye(:),stop(:));

        subplot(2,Nstim,s); hold off;
        plot(targ,eye,'.','color',colr(s,:)); hold on;
        plot(-40:40,-40:40,'k--'); plot(-40:40,(-40:40)*g.hor.targ{s}(i,t),'r');plot(-40:40,(-40:40)*b_targ,'r--');
        xlabel('target position [deg]'); ylabel('eye position [deg]'); axis([-40 40 -40 40]); vline(0,'k-'); hline(0,'k-');
        title(['t = ' num2str(t/60,'%.2f') ' s']);
        
        subplot(2,Nstim,s+Nstim); hold off;
        plot(stop,eye,'.','color',colr(s,:)*0.5); hold on;
        plot(-40:40,-40:40,'k--'); plot(-40:40,(-40:40)*g.hor.stop{s}(i,t),'r');plot(-40:40,(-40:40)*b_stop,'r--');
        xlabel('stop position [deg]'); ylabel('eye position [deg]'); axis([-40 40 -40 40]); vline(0,'k-'); hline(0,'k-');
        
        pause(0.05)
        end
    end
end

end