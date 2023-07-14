
clc;
clear
%% extract session

subject.trials = []; subject.behv_stats = []; subject.monkeyname = [];

monkeyname{1} = 'Quigley';
sessionfold{1} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{1} = {'Apr 03 2018','Apr 27 2018'};

monkeyname{2} = 'Bruno';
sessionfold{2} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{2} = {'Aug 03 2017','Aug 05 2017'};


monkeyname{3} = 'Schro';
sessionfold{3} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{3} = {'Feb 08 2018','Nov 29 2017'};

monkeyname{4} = 'Marco';
sessionfold{4} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{4} = {'Feb 11 2021','Feb 22 2021'};


monkeyname{5} = 'Viktor';
sessionfold{5} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{5} = {'Feb 03 2021','Feb 05 2021'};

monkeyname{6} = 'Jimmy';
sessionfold{6} = 'Z:\Data\Monkey2_newzdrive';
sessiondate{6} = {'Mar 10 2021','May 27 2021','Jul 30 2021'};


for i = 1:numel(monkeyname)

trials = []; behv_stats = [];
for j = 1:numel(sessiondate{i})
[trials_tmp,behv_stats_tmp] = extractMonkeySession(monkeyname{i},sessiondate{i}{j},sessionfold{i});
trials = [trials trials_tmp];
end

subject(i).trials = trials;
subject(i).monkeyname = monkeyname{i};

end

Nmax = 2000;

%% Get average timeseries

validtrls = [];
for i = 1:numel(subject)

trials = subject(i).trials;

% trial groups
FFon = arrayfun(@(x) logical(x.logical.firefly_fullON), trials);
try stimtrls = logical(arrayfun(@(x) x.logical.microstim, trials)); catch; stimtrls = false(1,numel(trials)); end
try rewarded = logical(arrayfun(@(x) x.logical.reward, trials)); catch; rewarded = ~isnan(arrayfun(@(x) x.events.t_rew, trials)); end
replay = logical(arrayfun(@(x) x.logical.replay, trials));
indxOn = FFon(:) & ~stimtrls(:) & rewarded(:) & ~replay(:);
indxOff = ~FFon(:) & ~stimtrls(:) & rewarded(:) & ~replay(:);
subject(i).trlgroups.indxOff = indxOff;
subject(i).trlgroups.indxOn = indxOn;

% extract trials
ts = arrayfun(@(x) x.continuous.ts, trials, 'un',0);
eye_hor = arrayfun(@(x) sign(nanmean(x.continuous.xfp))*nanmean([x.continuous.yle x.continuous.yre],2), trials, 'un',0);
eye_ver = arrayfun(@(x) nanmean([x.continuous.zle x.continuous.zre],2), trials, 'un',0);
lin_vel = arrayfun(@(x) x.continuous.v, trials,'un',0);
ang_vel = arrayfun(@(x) sign(nanmean(x.continuous.xfp))*x.continuous.w, trials,'un',0);

% save trials with equal lengths
dt = diff(trials(1).continuous.ts(1:2));
T = -0.2;
t_beg_i = arrayfun(@(x) find(x.continuous.ts >= T,1),trials);
ts = T:dt:T+Nmax*dt; ts(end) = [];

[eye_hor,eye_ver,lin_vel,ang_vel] = equalizeTrialLengths(t_beg_i,Nmax,eye_hor,eye_ver,lin_vel,ang_vel);

validtrls.ts(i,:) = ts;
validtrls.eye_hor.mu(i,:) = nanmean(eye_hor(indxOff,:));
validtrls.eye_ver.mu(i,:) = nanmean(eye_ver(indxOff,:));
validtrls.lin_vel.mu(i,:) = nanmean(lin_vel(indxOff,:));
validtrls.ang_vel.mu(i,:) = nanmean(ang_vel(indxOff,:));

validtrls.eye_hor.all{i} = eye_hor(indxOff,:);
validtrls.eye_ver.all{i} = eye_ver(indxOff,:);
validtrls.lin_vel.all{i} = lin_vel(indxOff,:);
validtrls.ang_vel.all{i} = ang_vel(indxOff,:);

% task params
t_stop = arrayfun(@(x) x.events.t_stop, trials);
t_move = arrayfun(@(x) x.events.t_move, trials); 
t_move_post = t_move; t_move_post(t_move_post < 0) = 0;
start_i = arrayfun(@(x) find(x.continuous.ts >=0,1), trials);
stop_i = arrayfun(@(x,t) find(x.continuous.ts >= t,1), trials, t_stop);

xmp0 = columnize(arrayfun(@(x) x.continuous.xmp(find(~isnan(x.continuous.xmp),1)), trials));
ymp0 = columnize(arrayfun(@(x) x.continuous.ymp(find(~isnan(x.continuous.ymp),1)), trials));
xmp1 = columnize(arrayfun(@(x,i) x.continuous.xmp(i), trials,stop_i));
ymp1 = columnize(arrayfun(@(x,i) x.continuous.ymp(i), trials,stop_i));
[r_sub,th_sub] = cart2polarY(xmp1-xmp0,ymp1-ymp0);

try
xfp = arrayfun(@(x) x.prs.xfp, trials); xfp = xfp(:) - median(xmp0);
yfp = arrayfun(@(x) x.prs.yfp, trials); yfp = yfp(:) - median(ymp0);
catch
xfp = arrayfun(@(x) nanmean(x.continuous.xfp), trials); xfp = xfp(:) - median(xmp0);
yfp = arrayfun(@(x) nanmean(x.continuous.yfp), trials); yfp = yfp(:) - median(ymp0);
end
[r_tar,th_tar] = cart2polarY(xfp,yfp);


trlDur = t_stop;
try taskParams(i).vmax = unique(arrayfun(@(x) x.prs.v_max, trials)); catch; taskParams(i).vmax = unique(arrayfun(@(x) x.prs.vmax, trials)); end
try taskParams(i).wmax = unique(arrayfun(@(x) x.prs.w_max, trials)); catch; taskParams(i).wmax = unique(arrayfun(@(x) x.prs.wmax, trials)); end

taskParams(i).minDist = round(min(r_tar));
taskParams(i).maxDist = round(max(r_tar));
taskParams(i).meanDist = round(nanmedian(r_tar));
taskParams(i).minAng = round(min(th_tar));
taskParams(i).maxAng = round(max(th_tar));
taskParams(i).meanTrlDur = nanmean(trlDur(indxOff) - t_move_post(indxOff));
taskParams(i).sdTrlDur = nanstd(trlDur(indxOff) - t_move_post(indxOff));
taskParams(i).meanMoveOnset = nanmedian(t_move(indxOff));
taskParams(i).TrlDur_2m = nanmedian(trlDur(indxOff & (r_tar(:) > 180 & r_tar(:) <220)));
taskParams(i).biasLin = regress(r_sub(indxOff),r_tar(indxOff));
taskParams(i).biasAng = regress(th_sub(indxOff),th_tar(indxOff));


end

%% Plot all trials
Nsubs = numel(subject);
ts = validtrls.ts(1,:);
Nskip = 20;
Ntrls_max = min(cellfun(@(x) size(x,1), validtrls.eye_hor.all));

% plot
figure;

for s = 1:Nsubs
    
    
    % Horizontal
    subplot(Nsubs,2,2*s-1); hold on; axis([0 3 -30 30]);
    ylabel('Horizontal [deg]'); xlabel('time [s]'); title('horizontal eye');
    for i = 1:Nskip:Ntrls_max
        plot(ts,validtrls.eye_hor.all{s}(i,:),'color',[0 0 0 0.3]);
    end
    legend(subject(s).monkeyname,'location','southeast');

    % Vertical
    subplot(Nsubs,2,2*s); hold on; axis([0 3 -40 10]);
    ylabel('Horizontal [deg]'); xlabel('time [s]'); title('vertical eye');
    for i = 1:Nskip:Ntrls_max
        plot(ts,validtrls.eye_ver.all{s}(i,:),'color',[0 0 0 0.3]);
    end
    
end

%% Plot averages

ts = validtrls.ts(1,:);

figure;

subplot(2,2,1); hold on; axis([-0.2 2.5 -10 20]);
xlabel('time (s)'); ylabel('horizontal eye (deg)'); title('horizontal eye');

subplot(2,2,2); hold on; axis([-0.2 2.5 -30 10]);
xlabel('time (s)'); ylabel('vertical eye (deg)');  title('vertical eye');

subplot(2,2,3); hold on; axis([-0.2 2.5 -10 200]);
xlabel('time (s)'); ylabel('linear velocity (cm/s)');  title('linear velocity');

subplot(2,2,4); hold on; axis([-0.2 2.5 -10 50]);
xlabel('time (s)'); ylabel('angular velocity (deg/s)'); title('angular velocity');

for i = 1:numel(subject)
    legsub{i} = subject(i).monkeyname;
    
    if any(strcmpi(subject(i).monkeyname,{'Jimmy','Viktor','Marco'}))

    subplot(2,2,1); hold on; 
    plot(ts,validtrls.eye_hor.mu(i,:),'linewidth',2); 
        
    subplot(2,2,2); hold on;
    plot(ts,validtrls.eye_ver.mu(i,:),'linewidth',2); 
    
    subplot(2,2,3); hold on; 
    plot(ts,validtrls.lin_vel.mu(i,:),'linewidth',2); 
    
    subplot(2,2,4); hold on; 
    plot(ts,validtrls.ang_vel.mu(i,:),'linewidth',2); 
    
    else % Baylor monkeys

    subplot(2,2,1); hold on;
    plot(ts,validtrls.eye_hor.mu(i,:),'linewidth',2,'color',[0 0 0 0.5]); 
    
    subplot(2,2,2); hold on;
    plot(ts,validtrls.eye_ver.mu(i,:),'linewidth',2,'color',[0 0 0 0.5]); 
    
    subplot(2,2,3); hold on;
    plot(ts,validtrls.lin_vel.mu(i,:),'linewidth',2,'color',[0 0 0 0.5]); 

    subplot(2,2,4); hold on;
    plot(ts,validtrls.ang_vel.mu(i,:),'linewidth',2,'color',[0 0 0 0.5]); 
    end

end
subplot(2,2,1); legend(legsub);
subplot(2,2,4); legend(legsub);

%% Saccade frequency

saccadeRate = [];
for s = 1:Nsubs

Ntrls = sum(subject(s).trlgroups.indxOff);
trlindx = subject(s).trlgroups.indxOff;
t_sac = cell2mat(arrayfun(@(x) x.events.t_sac(:), subject(s).trials(trlindx),'un',0)');

% start aligned
tmin = -0.5;
tmax = 2.1;
wn = 0.1;
stp = 0.01;

[y,x] = movWnHist(tmin,tmax,stp,wn,t_sac);

y = y/wn/Ntrls;
saccadeRate.val(s,:) = y(2:end-1);
saccadeRate.t = x(2:end-1);
end




figure; hold on; axis([-0.2 1.2 0 10]);
xlabel('time (s)'); ylabel('Saccade frequency (Hz)'); title('Saccade frequency');

for i = 1:numel(subject)
    legsub{i} = subject(i).monkeyname;
    
    if any(strcmpi(subject(i).monkeyname,{'Jimmy','Viktor','Marco'}))
    
    plot(saccadeRate.t,saccadeRate.val(i,:),'linewidth',2); 
    
    else % Baylor monkeys

    plot(saccadeRate.t,saccadeRate.val(i,:),'linewidth',2,'color',[0 0 0 0.5]); 

    end

end

legend(legsub);

%% Saccade amplitude

if 0
figure;

for s = 1:Nsubs
    
    flag = 0;
    trlindx = find(subject(s).trlgroups.indxOff);

    for i = 1:numel(trlindx)
    
    subplot(2,1,1); hold on;
    plot(ts,validtrls.eye_hor.all{s}(i,:)); vline(subject(s).trials(trlindx(i)).events.t_sac,'r');
    xlabel('time (s)'); ylabel('horizontal eye (deg)'); title([subject(s).monkeyname ', trial ' num2str(trlindx(i))])

    subplot(2,1,2); hold on;
    plot(ts,validtrls.eye_ver.all{s}(i,:)); vline(subject(s).trials(trlindx(i)).events.t_sac,'r');
    xlabel('time (s)'); ylabel('vertical eye (deg)');

    keyboard; clf;

    if flag; break; end

    end

end
end


presac = round(0.05/dt);
postsac = round(0.15/dt);

epoch_steer = [0.5 1.5];
epoch_targON = [0 0.3];

saccadeMagnitude = [];
for s = 1:Nsubs
    
    trlindx = find(subject(s).trlgroups.indxOff);
    Ntrls = numel(trlindx);
    
    tsac = []; sacmag = [];
    for i = 1:Ntrls
        tsac_new = subject(s).trials(trlindx(i)).events.t_sac;
        Nt = numel(subject(s).trials(trlindx(i)).continuous.ts);

        if ~isempty(tsac_new)
            
            hor = nanmean([subject(s).trials(trlindx(i)).continuous.yle subject(s).trials(trlindx(i)).continuous.yre],2);
            ver = nanmean([subject(s).trials(trlindx(i)).continuous.zle subject(s).trials(trlindx(i)).continuous.zre],2);
            tsac_i = arrayfun(@(t) find(subject(s).trials(trlindx(i)).continuous.ts >= t,1), tsac_new);
            rmindx = tsac_i <= presac | tsac_i >= Nt-postsac-presac;
            tsac_i(rmindx) = [];
            tsac_new(rmindx) = [];

            sacmag_new = arrayfun(@(t) sqrt( [nanmean(hor(t-presac:t)) -  nanmean(hor(t+postsac:t+postsac+presac))].^2 + ...
                [nanmean(ver(t-presac:t)) -  nanmean(ver(t+postsac:t+postsac+presac))].^2), tsac_i);
            
            tsac = [tsac  ; tsac_new(:)];
            sacmag = [sacmag ; sacmag_new(:)];
        end
    end
    
    % select trial epoch
    tsac_targON = tsac(tsac>epoch_targON(1) & tsac<epoch_targON(2));
    sacmag_targON = sacmag(tsac>epoch_targON(1) & tsac<epoch_targON(2));
    
    tsac_steer = tsac(tsac>epoch_steer(1) & tsac<epoch_steer(2));
    sacmag_steer = sacmag(tsac>epoch_steer(1) & tsac<epoch_steer(2));
    
    % saccade magnitude cdf (steering)
    nx = -1:41;
    [ny,nx] =  hist(sacmag_steer,nx);
    ny = cumsum(ny,'omitnan')/nansum(ny);
    saccadeMagnitude.steering.val(s,:) = ny(2:end-1);
    saccadeMagnitude.steering.x = nx(2:end-1);
    
    % saccade magnitude cdf (target ON)
    nx = -1:41;
    [ny,nx] =  hist(sacmag_targON,nx);
    ny = cumsum(ny,'omitnan')/nansum(ny);
    saccadeMagnitude.targON.val(s,:) = ny(2:end-1);
    saccadeMagnitude.targON.x = nx(2:end-1);
    
end


figure; 

subplot(1,2,1); hold on; 
xlabel('time (s)'); ylabel('Saccade magnitude (deg)'); title('Target ON');

subplot(1,2,2); hold on; 
xlabel('time (s)'); ylabel('Saccade magnitude (deg)'); title('Steering');

for i = 1:numel(subject)
    legsub{i} = subject(i).monkeyname;
    
    if any(strcmpi(subject(i).monkeyname,{'Jimmy','Viktor','Marco'}))

    subplot(1,2,1); hold on;
    plot(saccadeMagnitude.targON.x,saccadeMagnitude.targON.val(i,:),'linewidth',2);
    subplot(1,2,2); hold on;
    plot(saccadeMagnitude.steering.x,saccadeMagnitude.steering.val(i,:),'linewidth',2); 
    
    else % Baylor monkeys
    subplot(1,2,1); hold on;
    plot(saccadeMagnitude.targON.x,saccadeMagnitude.targON.val(i,:),'linewidth',2,'color',[0 0 0 0.5]); 
    subplot(1,2,2); hold on;
    plot(saccadeMagnitude.steering.x,saccadeMagnitude.steering.val(i,:),'linewidth',2,'color',[0 0 0 0.5]); 
    end

end

legend(legsub,'location','southeast');

%% Average velocity


figure; hold on;

for s = 1:Nsubs
    
    linVelmu(s) = nanmean(nanmean(validtrls.lin_vel.all{s},2));
    linVelsd(s) = nanstd(nanmean(validtrls.lin_vel.all{s},2));
    angVelmu(s) = nanmean(nanmean(validtrls.ang_vel.all{s},2));
    angVelsd(s) = nanstd(nanmean(validtrls.ang_vel.all{s},2));

    if any(strcmpi(subject(s).monkeyname,{'Jimmy','Viktor','Marco'}))

    bar(s,linVelmu(s)); errorbar(s,linVelmu(s),linVelsd(s),'k','capsize',0); 

    else % Baylor monkeys
    
    bar(s,linVelmu(s),'facecolor',[.5 .5 .5]); errorbar(s,linVelmu(s),linVelsd(s),'k','capsize',0); 

    end

end

xticks(1:Nsubs); xticklabels(legsub); ylabel('Linear Velocity (cm/s)'); title('Average velocity');

%% Average trial duration
vel_thresh = 10;

figure; hold on;

for s = 1:Nsubs
    
    Ntrls = size(validtrls.lin_vel.all{s},1);
    Nt = size(validtrls.lin_vel.all{s},2);
    ts = validtrls.ts(1,:);
    startindx = find(ts>=0,1);

    moveindx = arrayfun(@(i) find(validtrls.lin_vel.all{s}(i,:) > vel_thresh, 1), 1:Ntrls); 
    moveindx = moveindx-startindx; moveindx(moveindx<=0)=1;
    stopindx = arrayfun(@(i) Nt-find(flip(validtrls.lin_vel.all{s}(i,:)) > vel_thresh, 1), 1:Ntrls);
    
    trlength_mu(s) = nanmean(dt*(stopindx-moveindx));
    trlength_sd(s) = nanstd(dt*(stopindx-moveindx));

%     stopindx = arrayfun(@(i) find(~isnan(validtrls.lin_vel.all{s}(i,:)),1,'last'), 1:Ntrls);
%     trlength_mu(s) = nanmean(ts(lenindx));
%     trlength_sd(s) = nanstd(ts(lenindx));

    if any(strcmpi(subject(s).monkeyname,{'Jimmy','Viktor','Marco'}))

    bar(s,trlength_mu(s)); errorbar(s,trlength_mu(s),trlength_sd(s),'k','capsize',0); 

    else % Baylor monkeys
    
    bar(s,trlength_mu(s),'facecolor',[.5 .5 .5]); errorbar(s,trlength_mu(s),trlength_sd(s),'k','capsize',0); 

    end

end

xticks(1:Nsubs); xticklabels(legsub); ylabel('Trial duration (s)'); title('Average Trial duration');


figure; hold on;

for s = 1:Nsubs
    
    if any(strcmpi(subject(s).monkeyname,{'Jimmy','Viktor','Marco'}))

    bar(s,taskParams(s).meanTrlDur); errorbar(s,taskParams(s).meanTrlDur,taskParams(s).sdTrlDur,'k','capsize',0); 

    else % Baylor monkeys
    
    bar(s,taskParams(s).meanTrlDur,'facecolor',[.5 .5 .5]); errorbar(s,taskParams(s).meanTrlDur,taskParams(s).sdTrlDur,'k','capsize',0); 

    end

end

xticks(1:Nsubs); xticklabels(legsub); ylabel('Trial duration (s)'); title('Average Trial duration');

%% Bias


figure; hold on;

for s = 1:Nsubs
    
    if any(strcmpi(subject(s).monkeyname,{'Jimmy','Viktor','Marco'}))

    bar(s,taskParams(s).biasLin); 

    else % Baylor monkeys
    
    bar(s,taskParams(s).biasLin,'facecolor',[.5 .5 .5]); 

    end

end

xticks(1:Nsubs); xticklabels(legsub); ylabel('Linear bias'); title('Linear bias'); ylim([0.6 1]);























%% Quigley

monkeyname = 'Quigley';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Apr 03 2018','Apr 27 2018'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 10;


for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

%% Ody
% 
% monkeyname = 'Ody';
% sessionfold = 'Z:\Data\Monkey2_newzdrive';
% sessiondate = {'Aug 21 2019','Oct 17 2019'};
% 
% % define filter
% sig = 8; %filter width
% sz = 2*sig; %filter size
% t2 = linspace(-sz/2, sz/2, sz);
% h = exp(-t2.^2/(2*sig^2));
% h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
% 
% Nskip = 20;
% 
% for i = 1:numel(sessiondate)
% [taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
% taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
% end

%% Bruno

monkeyname = 'Bruno';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Aug 03 2017','Aug 05 2017'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 20;

for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

%% Schro

monkeyname = 'Schro';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Feb 08 2018','Nov 29 2017'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 5;

for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

%% Marco

monkeyname = 'Marco';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Feb 11 2021','Feb 22 2021'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 20;


for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

%% Viktor

monkeyname = 'Viktor';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Feb 03 2021','Feb 05 2021'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 20;


for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

%% Jimmy - Spike 2

monkeyname = 'Jimmy';
sessionfold = 'Z:\Data\Monkey2_newzdrive';
sessiondate = {'Mar 10 2021','May 27 2021','Jul 30 2021'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 20;

for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end


%% Jimmy - Unity (Stimulation)

monkeyname = 'Jimmy';
sessionfold = 'W:\Monkeys';
sessiondate = {'Feb 28 2023'};%,'Mar 21 2023','Mar 22 2023'};

% define filter
sig = 8; %filter width
sz = 2*sig; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = 1;%h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

Nskip = 5;

for i = 1:numel(sessiondate)
[taskParams.(monkeyname){i},eyeTrace.(monkeyname){i},vel.(monkeyname){i},ts] = plotSessionEyeMovements(monkeyname,sessiondate{i},sessionfold,Nskip,h);
taskParams.(monkeyname){i}.sessiondate = sessiondate{i};
end

