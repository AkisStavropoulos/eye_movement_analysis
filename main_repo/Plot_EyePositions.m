function Plot_EyePositions(tracking_regular,subjectnames)

%% eye position as a function of time, distance %, or absolute distance

Nsubs = size(tracking_regular,1);
Nstim = size(tracking_regular,2);

colr = brewermap(Nstim,'Dark2');

for i = 1:Nsubs
    ver = figure('name',[subjectnames{i} ': vertical eye positions'],'numbertitle','off','position',[0 500 900 500]);
    hor = figure('name',[subjectnames{i} ': horizontal eye positions'],'numbertitle','off','position',[0 5 900 500]);
    for s = 1:Nstim
        % time
        ts = tracking_regular{i,s}.misc.ts;
        % distance to end
        maxD = max(tracking_regular{i,s}.misc.subdist,[],2);
        dist = tracking_regular{i,s}.misc.subdist;
        dist2end = dist - maxD;
        % distance percentage
        distperc = mat2cell(dist,ones(size(dist,1),1),numel(ts));
        distperc = cell2mat(cellfun(@(x,m) x./m, distperc, num2cell(maxD),'uniformoutput',false));

        % vertical
        figure(ver); 
        subplot(4,Nstim,s); hold on;
        plot(ts,tracking_regular{i,s}.eyepos.screen.ver_mean,'color',colr(s,:),'linewidth',0.1); xlabel('time [s]');ylabel('vertical eye [deg]'); axis([0 12 -60 20]);
        
        subplot(4,Nstim,s+Nstim); hold on;
        plot(dist',tracking_regular{i,s}.eyepos.screen.ver_mean','color',colr(s,:)*0.5,'linewidth',0.1); xlabel('distance [cm]');ylabel('vertical eye [deg]'); axis([0 600 -60 20]);
        
        subplot(4,Nstim,s+2*Nstim); hold on;
        plot(dist2end',tracking_regular{i,s}.eyepos.screen.ver_mean','color',colr(s,:)*0.5,'linewidth',0.1); xlabel('distance to end [cm]');ylabel('vertical eye [deg]'); axis([-600 0 -60 20]);
        
        subplot(4,Nstim,s+3*Nstim); hold on;
        plot(distperc',tracking_regular{i,s}.eyepos.screen.ver_mean','color',colr(s,:)*0,'linewidth',0.1); xlabel('distance %');ylabel('vertical eye [deg]'); axis([0 1 -60 20]);
        % horizontal
        figure(hor);
        subplot(4,Nstim,s); hold on;
        plot(ts,tracking_regular{i,s}.eyepos.screen.hor_mean,'color',colr(s,:),'linewidth',0.1); xlabel('time [s]');ylabel('horizontal eye [deg]'); axis([0 12 -40 40]);
        
        subplot(4,Nstim,s+Nstim); hold on;
        plot(dist',tracking_regular{i,s}.eyepos.screen.hor_mean','color',colr(s,:)*0.5,'linewidth',0.1); xlabel('distance [cm]');ylabel('horizontal eye [deg]'); axis([0 600 -40 40]);
        
        subplot(4,Nstim,s+2*Nstim); hold on;
        plot(dist2end',tracking_regular{i,s}.eyepos.screen.hor_mean','color',colr(s,:)*0.5,'linewidth',0.1); xlabel('distance to end [cm]');ylabel('horizontal eye [deg]'); axis([-600 0 -40 40]);
        
        subplot(4,Nstim,s+3*Nstim); hold on;
        plot(distperc',tracking_regular{i,s}.eyepos.screen.hor_mean','color',colr(s,:)*0,'linewidth',0.1); xlabel('distance %');ylabel('horizontal eye [deg]'); axis([0 1 -40 40]);
    end
end
