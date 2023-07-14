function [trials,taskParams,eyeTrace,vel,ts] = plotSessionEyeMovements(monkeyname,sessiondate,sessionfold,Nskip,eyePosFilter)

% choose data folder

if strcmpi(sessionfold(1),'Z')
    recSpec = 'U-probe';
    if strcmpi(monkeyname,'Quigley')
        recSpec = 'Sim_recordings';
elseif strcmpi(monkeyname,'Ody') || strcmpi(monkeyname,'Bruno') 
        recSpec = 'Utah Array';
elseif strcmpi(monkeyname,'Schro') && strcmpi(sessiondate,'Feb 08 2018')
        recSpec = 'Utah Array';
    end
elseif strcmpi(sessionfold(1),'W')
    recSpec = 'Stimulation';
end

datafold = fullfile(sessionfold,monkeyname,recSpec,sessiondate,'neural data\Pre-processing X E');

% load data
fname = dir(fullfile(datafold,'m*s*.mat'));
cd(datafold) 
load(fname(1).name)

trials = trials_behv;

fprintf('Monkey: %s\nSession: %s \nNtrials = %d\n',monkeyname,fname(1).name,numel(trials));

% select trialgroup
% stop_err = cellfun(@(x,y) sqrt(x(end).^2 + y(end).^2) ,behv_stats.pos_rel.x_targ,behv_stats.pos_rel.y_targ)';
% r_tar = behv_stats.pos_final.r_targ';
% validtrls = stop_err < 0.5*r_tar;

FFon = logical(arrayfun(@(x) x.logical.firefly_fullON, trials));
try stimtrls = logical(arrayfun(@(x) x.logical.microstim, trials)); catch; stimtrls = false(1,numel(trials)); end
try rewarded = logical(arrayfun(@(x) x.logical.reward, trials)); catch; rewarded = ~isnan(arrayfun(@(x) x.events.t_rew, trials)); end
replay = logical(arrayfun(@(x) x.logical.replay, trials));
indxOn = FFon(:) & ~stimtrls(:) & rewarded(:) & ~replay(:);
indxOff = ~FFon(:) & ~stimtrls(:) & rewarded(:) & ~replay(:);

% extract trials
Nmax = 2000;
ts_off = arrayfun(@(x) x.continuous.ts, trials(indxOff), 'un',0);
eye_hor_off = arrayfun(@(x) conv(nanmean([x.continuous.yle x.continuous.yre],2),eyePosFilter,'same'), trials(indxOff), 'un',0);
eye_ver_off = arrayfun(@(x) conv(nanmean([x.continuous.zle x.continuous.zre],2),eyePosFilter,'same'), trials(indxOff), 'un',0);

ts_on = arrayfun(@(x) x.continuous.ts, trials(indxOn), 'un',0);
eye_hor_on = arrayfun(@(x) conv(nanmean([x.continuous.yle  x.continuous.yre],2),eyePosFilter,'same'), trials(indxOn), 'un',0);
eye_ver_on = arrayfun(@(x) conv(nanmean([x.continuous.zle  x.continuous.zre],2),eyePosFilter,'same'), trials(indxOn), 'un',0);

% save trials with equal lengths
dt = diff(trials(1).continuous.ts(1:2));
T = -0.2;
t_beg_i = arrayfun(@(x) find(x.continuous.ts >= T,1),trials);
ts = T:dt:T+Nmax*dt; ts(end) = [];

[eyeTrace.targetOFF.hor,eyeTrace.targetOFF.ver] = equalizeTrialLengths ...
                    (t_beg_i(indxOff),Nmax,eye_hor_off,eye_ver_off);

[eyeTrace.targetON.hor,eyeTrace.targetON.ver] = equalizeTrialLengths ...
                    (t_beg_i(indxOn),Nmax,eye_hor_on,eye_ver_on);

% basic stuff about task
xmp0 = arrayfun(@(x) x.continuous.xmp(find(~isnan(x.continuous.xmp),1)), trials);
ymp0 = arrayfun(@(x) x.continuous.ymp(find(~isnan(x.continuous.ymp),1)), trials);
try
xfp = arrayfun(@(x) x.prs.xfp, trials); xfp = xfp(:) - median(xmp0);
yfp = arrayfun(@(x) x.prs.yfp, trials); yfp = yfp(:) - median(ymp0);
catch
xfp = arrayfun(@(x) nanmean(x.continuous.xfp), trials); xfp = xfp(:) - median(xmp0);
yfp = arrayfun(@(x) nanmean(x.continuous.yfp), trials); yfp = yfp(:) - median(ymp0);
end
[r_tar,th_tar] = cart2polarY(xfp,yfp);

trlDur = arrayfun(@(x) x.events.t_end, trials(indxOff));
try taskParams.vmax = unique(arrayfun(@(x) x.prs.v_max, trials)); catch; taskParams.vmax = unique(arrayfun(@(x) x.prs.vmax, trials)); end
try taskParams.wmax = unique(arrayfun(@(x) x.prs.w_max, trials)); catch; taskParams.wmax = unique(arrayfun(@(x) x.prs.wmax, trials)); end

taskParams.minDist = round(min(r_tar));
taskParams.maxDist = round(max(r_tar));
taskParams.minAng = round(min(th_tar));
taskParams.maxAng = round(max(th_tar));
taskParams.meanTrlDur = nanmedian(trlDur);


fprintf('Target distances:     %d to %d cm\n',taskParams.minDist,taskParams.maxDist);
fprintf('Target angles:        %d to %d deg\n',taskParams.minAng,taskParams.maxAng);
fprintf('Mean Trial Duration:  %.2f s \n',taskParams.meanTrlDur);
fprintf('Max Linear Velocity:  %d cm/s \n',taskParams.vmax);
fprintf('Max Angular Velocity: %d deg/s \n\n',taskParams.wmax);

% extract velocity
vel_tmp = arrayfun(@(x) x.continuous.v, trials,'un',0);
vel_on = vel_tmp(indxOn);
vel_off = vel_tmp(indxOff);
[vel.targetOFF] = equalizeTrialLengths(t_beg_i(indxOff),Nmax,vel_tmp(indxOff));
[vel.targetON] = equalizeTrialLengths(t_beg_i(indxOn),Nmax,vel_tmp(indxOn));


% plot
figure;

% ffOFF
subplot(3,2,1); hold on; axis([0 3 -30 30]);
ylabel('Horizontal [deg]'); xlabel('time [s]'); title('FF off');
for i = 1:Nskip:sum(indxOff)
plot(ts_off{i},eye_hor_off{i},'color',[0 0 0 0.3]);
end
subplot(3,2,3); hold on; axis([0 3 -40 10]);
ylabel('Vertical [deg]'); xlabel('time [s]');
for i = 1:Nskip:sum(indxOff)
plot(ts_off{i},eye_ver_off{i},'color',[0 0 0 0.3]);
end
subplot(3,2,5); hold on; axis([0 3 -5 205]);
ylabel('Linear velocity [cm/s]'); xlabel('time [s]'); 
for i = 1:round(Nskip/10):sum(indxOn)
plot(ts_off{i},vel_off{i},'color',[0 0 0 0.3]);
end

% ffON
subplot(3,2,2); hold on; axis([0 3 -30 30]);
ylabel('Horizontal [deg]'); xlabel('time [s]'); title('FF on');
for i = 1:round(Nskip/10):sum(indxOn)
plot(ts_on{i},eye_hor_on{i},'color',[0 0 0 0.3]);
end
subplot(3,2,4); hold on; axis([0 3 -40 10]);
ylabel('Vertical [deg]'); xlabel('time [s]');
for i = 1:round(Nskip/10):sum(indxOn)
plot(ts_on{i},eye_ver_on{i},'color',[0 0 0 0.3]);
end
subplot(3,2,6); hold on; axis([0 3 -5 205]);
ylabel('Linear velocity [cm/s]'); xlabel('time [s]');
for i = 1:round(Nskip/10):sum(indxOn)
plot(ts_on{i},vel_on{i},'color',[0 0 0 0.3]);
end


sgtitle([monkeyname ' - ' sessiondate])

