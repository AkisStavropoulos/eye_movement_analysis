function [varexp_start,varexp_stop] = TargetVsStopTTI(tracking,prs)

%% Target-tracking index comparison (target vs stop position)
if strcmp(tracking{1}.misc.sampleflag,'time')
    T_start = 6; T_stop = 4; T_break = 2;  dt = diff(tracking{1}.misc.ts(1:2)); % dt = 1/60; % NOT SMR SAMPLING RATE, it's is downsampled to 60Hz
elseif strcmp(tracking{1}.misc.sampleflag,'distperc')
    T_start = 1; T_stop = 1; T_break = 0; dt = diff(tracking{1}.misc.ts(1:2)); % dt = dd_distperc; % for distance % sampling
elseif strcmp(tracking{1}.misc.sampleflag,'distance2end')
    T_start = 300; T_stop = 300; T_break = 20;  dt = diff(tracking{1}.misc.ts(1:2)); % dt = dd_dist2end; % for distance % sampling
end
nt_start = round(T_start/dt);
nt_stop = round(T_stop/dt);
nt_break = round(T_break/dt);

Nsubs = size(tracking,1);
Nstim = size(tracking,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
xaxislabel = tracking{1}.misc.sampleflag;

maxlength = min(cellfun(@(x) numel(x.eyepos.pred_vs_true.targ.var_explained.mu.startaligned), tracking),[],'all');

varexp_start = []; varexp_stop = []; sem_start = []; sem_stop = [];
for i = 1:Nsubs
    for s = 1:Nstim
        % target position
        varexp_start.targ{i,s} = tracking{i,s}.eyepos.pred_vs_true.targ.var_explained.mu.startaligned(1:maxlength);
        varexp_stop.targ{i,s} = tracking{i,s}.eyepos.pred_vs_true.targ.var_explained.mu.stopaligned(1:maxlength);
        varexp_start.targ{i,s}(varexp_start.targ{i,s} <= 0) = 0;
        varexp_stop.targ{i,s}(varexp_stop.targ{i,s} <= 0) = 0;
        startLNG(i,s) = length(varexp_start.targ{i,s});
        stopLNG(i,s) = length(varexp_stop.targ{i,s});
        trackdur.targ(i,s) = dt*find(varexp_start.targ{i,s} == 0,1);
        % stop position
        varexp_start.stop{i,s} = tracking{i,s}.eyepos.pred_vs_true.stop.var_explained.mu.startaligned(1:maxlength);
        varexp_stop.stop{i,s} = tracking{i,s}.eyepos.pred_vs_true.stop.var_explained.mu.stopaligned(1:maxlength);
        varexp_start.stop{i,s}(varexp_start.stop{i,s} <= 0) = 0;
        varexp_stop.stop{i,s}(varexp_stop.stop{i,s} <= 0) = 0;
        startLNG(i,s) = length(varexp_start.stop{i,s});
        stopLNG(i,s) = length(varexp_stop.stop{i,s});
        trackdur.stop(i,s) = dt*find(varexp_start.stop{i,s} == 0,1);

        if prs.boots
        sem_start.targ{i,s} = tracking{i,s}.eyepos.pred_vs_true.targ.var_explained.sem.startaligned(1:maxlength); sem_start.targ{i,s}(varexp_start.targ{i,s} <= 0) = 0;
        sem_stop.targ{i,s} = tracking{i,s}.eyepos.pred_vs_true.targ.var_explained.sem.stopaligned(1:maxlength);   sem_stop.targ{i,s}(varexp_stop.targ{i,s} <= 0) = 0;
        sem_start.stop{i,s} = tracking{i,s}.eyepos.pred_vs_true.stop.var_explained.sem.startaligned(1:maxlength); sem_start.stop{i,s}(varexp_start.stop{i,s} <= 0) = 0;
        sem_stop.stop{i,s} = tracking{i,s}.eyepos.pred_vs_true.stop.var_explained.sem.stopaligned(1:maxlength);   sem_stop.stop{i,s}(varexp_stop.stop{i,s} <= 0) = 0;
        end
    end
end

if 0
% subject separate
for i = 1:Nsubs
    figure('name',subject(i).name,'numbertitle','off','position',[0 700 950 250]);hold on;
    for s = 1:Nstim
        subplot(1,Nstim,s); hold on;
        % target position
        if prs.boots
        shadedErrorBar(dt:dt:nt_start*dt, sqrt(varexp_start.targ{i,s}(1:nt_start)),sem_start.targ{i,s}(1:nt_start),'lineprops',{'Color',colr(s,:)});
        shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.targ{i,s}(1:nt_stop))),flip(sem_stop.targ{i,s}(1:nt_stop)),'lineprops',{'Color',colr(s,:),'HandleVisibility','off'});
        else
        plot(dt:dt:nt_start*dt, sqrt(varexp_start.targ{i,s}(1:nt_start)),'Color',colr(s,:));
        plot(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.targ{i,s}(1:nt_stop))),'Color',colr(s,:),'HandleVisibility','off');
        end
        % stop position
        if prs.boots
        shadedErrorBar(dt:dt:nt_start*dt, sqrt(varexp_start.stop{i,s}(1:nt_start)),sem_start.stop{i,s}(1:nt_start),'lineprops',{'Color',colr(s,:)*0.5});
        shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.stop{i,s}(1:nt_stop))),flip(sem_stop.stop{i,s}(1:nt_stop)),'lineprops',{'Color',colr(s,:)*0.5,'HandleVisibility','off'});
        else
        plot(dt:dt:nt_start*dt, sqrt(varexp_start.stop{i,s}(1:nt_start)),'Color',colr(s,:)*0.5);
        plot(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.stop{i,s}(1:nt_stop))),'Color',colr(s,:)*0.5,'HandleVisibility','off');
        end
        xlabel(xaxislabel); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
        legend('target position','stop position');
    end
end
end

% all subjects
figure('name','All Subjects: TTI','numbertitle','off','position',[0 700 950 250]);hold on;
for s = 1:Nstim
    MINindxSTART = min(startLNG(:,s));
    MINindxSTOP = min(stopLNG(:,s));
    
    % target position
    varexp_start.all.targ{s} = nanmean(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_start.targ(:,s),'uniformoutput',false)'),2);
    sem_start.all.targ{s} = nanstd(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_start.targ(:,s),'uniformoutput',false)'),[],2)./sqrt(Nsubs);
    varexp_stop.all.targ{s} = nanmean(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_stop.targ(:,s),'uniformoutput',false)'),2);
    sem_stop.all.targ{s} = nanstd(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_stop.targ(:,s),'uniformoutput',false)'),[],2)./sqrt(Nsubs);
    % stop position
    varexp_start.all.stop{s} = nanmean(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_start.stop(:,s),'uniformoutput',false)'),2);
    sem_start.all.stop{s} = nanstd(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_start.stop(:,s),'uniformoutput',false)'),[],2)./sqrt(Nsubs);
    varexp_stop.all.stop{s} = nanmean(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_stop.stop(:,s),'uniformoutput',false)'),2);
    sem_stop.all.stop{s} = nanstd(cell2mat(cellfun(@(x) x(1:MINindxSTART), varexp_stop.stop(:,s),'uniformoutput',false)'),[],2)./sqrt(Nsubs);
    
    subplot(1,Nstim,s); hold on;
    % target position
    shadedErrorBar(dt:dt:nt_start*dt, sqrt(varexp_start.all.targ{s}(1:nt_start)),sem_start.all.targ{s}(1:nt_start),'lineprops',{'Color',colr(s,:)});
    shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.all.targ{s}(1:nt_stop))),flip(sem_stop.all.targ{s}(1:nt_stop)),'lineprops',{'Color',colr(s,:),'HandleVisibility','off'});
    % stop position
    shadedErrorBar(dt:dt:nt_start*dt, sqrt(varexp_start.all.stop{s}(1:nt_start)),sem_start.all.stop{s}(1:nt_start),'lineprops',{'Color',colr(s,:)*0.5});
    shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(sqrt(varexp_stop.all.stop{s}(1:nt_stop))),flip(sem_stop.all.stop{s}(1:nt_stop)),'lineprops',{'Color',colr(s,:)*0.5,'HandleVisibility','off'});
    xlabel(xaxislabel); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
    legend('target position','stop position');
end

% difference
figure('name','All Subjects: TTI difference','numbertitle','off','position',[0 700 950 250]);hold on;
for s = 1:Nstim
    MINindxSTART = min(startLNG(:,s));
    MINindxSTOP = min(stopLNG(:,s));
    
    diff_start{s} = sqrt([varexp_start.targ{:,s}]) - sqrt([varexp_start.stop{:,s}]);
    diff_stop{s} = sqrt([varexp_stop.targ{:,s}]) - sqrt([varexp_stop.stop{:,s}]);
    diff_start_mu{s} = nanmean(diff_start{s},2);    diff_start_se{s} = nanstd(diff_start{s},[],2)./sqrt(Nsubs);
    diff_stop_mu{s} = nanmean(diff_stop{s},2);    diff_stop_se{s} = nanstd(diff_stop{s},[],2)./sqrt(Nsubs);

    subplot(1,Nstim,s); hold on;
    shadedErrorBar(dt:dt:nt_start*dt, diff_start_mu{s}(1:nt_start),diff_start_se{s}(1:nt_start),'lineprops',{'Color',colr(s,:)});
    shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt), flip(diff_stop_mu{s}(1:nt_stop)),flip(diff_stop_se{s}(1:nt_stop)),'lineprops',{'Color',colr(s,:),'HandleVisibility','off'});
    xlabel(xaxislabel); ylabel('TTI difference'); axis([0 T_start+T_stop+T_break -0.5 0.5]); hline(0,'k');
    legend('target - stop position TTI');
end

