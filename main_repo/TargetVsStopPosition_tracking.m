function tracking = TargetVsStopPosition_tracking(trials,rmsac,prs)

%% disregard trials without recording
dt = 1/60; % NOT SMR SAMPLING RATE, it is downsampled to 60 Hz
NORECindx = [];
for j = 1:length(trials)
    if trials(j).prs.NOREC
        NORECindx = [NORECindx j];
    end
end
trials(NORECindx) = [];

rmindx = [];
for j = 1:length(trials)
    EYE = [trials(j).continuous.zle(:) trials(j).continuous.zre(:)...
           trials(j).continuous.yle(:) trials(j).continuous.yre(:)];
   if any(std(EYE,[],1) <= 1) 
       rmindx = [rmindx j];
   end
   blinkindx = (trials(j).continuous.blink);
   nanindx = (isnan(sum(EYE,2)));
   if ((sum(blinkindx) + sum(nanindx)) >= 0.60*length(trials(j).continuous.ts))
       rmindx = [rmindx j];
   end
end
trials(rmindx) = [];

%% Set-up data
Ntrials = length(trials);
tau = zeros(Ntrials,1); w = cell(1,Ntrials); v = cell(1,Ntrials); x_monk = cell(1,Ntrials); y_monk = cell(1,Ntrials); phi = cell(1,Ntrials); x_fly = cell(1,Ntrials); y_fly = cell(1,Ntrials); x_stop = cell(1,Ntrials); y_stop = cell(1,Ntrials); d_sub = cell(1,Ntrials); zle = cell(1,Ntrials); zre = cell(1,Ntrials); yle = cell(1,Ntrials); yre = cell(1,Ntrials); t_sac = cell(1,Ntrials); ts = cell(1,Ntrials); t_stop = zeros(Ntrials,1); endofsteering = zeros(Ntrials,1); behverrors = zeros(Ntrials,1);
for j = 1:Ntrials
    tau(j) = trials(j).prs.tau;
    w{j} = trials(j).continuous.w(:);       v{j} = trials(j).continuous.v(:);   
    [x_monk{j}, y_monk{j}, ~, ~, ~, phi{j}] = gen_traj(w{j}, v{j}, trials(j).continuous.xmp(1), trials(j).continuous.ymp(1),dt);
    x_monk{j} = x_monk{j}(1:end-1);         y_monk{j} = y_monk{j}(1:end-1);
    phi{j} = phi{j}(1:end-1);
    
    x_fly{j} = trials(j).prs.fireflyposx - x_monk{j};     y_fly{j} = trials(j).prs.fireflyposy - y_monk{j};
    x_stop{j} = x_monk{j}(end) - x_monk{j};               y_stop{j} = y_monk{j}(end) - y_monk{j};

    R = @(phi) [cosd(phi)   sind(phi) ; -sind(phi)  cosd(phi)]; 
    XY_fly = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi{j}, x_fly{j}, y_fly{j},'uniformoutput',false));
    XY_stop = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi{j}, x_stop{j}, y_stop{j},'uniformoutput',false));
    x_fly{j} = XY_fly(:,1);    y_fly{j} = XY_fly(:,2);
    x_stop{j} = XY_stop(:,1);   y_stop{j} = XY_stop(:,2);

    d_sub{j} = cumsum(trials(j).continuous.v)*dt;
    zle{j} = trials(j).continuous.zle;    %zle{j}(trials(j).continuous.blink) = nan;
    yle{j} = trials(j).continuous.yle;    %yle{j}(trials(j).continuous.blink) = nan;
    zre{j} = trials(j).continuous.zre;    %zre{j}(trials(j).continuous.blink) = nan;
    yre{j} = trials(j).continuous.yre;    %yre{j}(trials(j).continuous.blink) = nan;
    t_sac{j} = trials(j).events.t_sac - trials(j).events.t_beg;
    t_stop(j) = trials(j).events.t_end - trials(j).events.t_beg;
    ts{j} = trials(j).continuous.ts;
    if sum(size(trials(j).mc.JS_X_Raw)) <= 2
        endofsteering(j) = nan;
    else
        endofsteering(j) = find(trials(j).mc.JS_X_Raw,1,'last');
    end
    behverrors(j) = sqrt(x_fly{j}(end).^2 + y_fly{j}(end).^2);
end
%% prs
TauAdjust = prs.TauAdjust; % adjust frame cutoffs to varying duration of stopping time caused by varying taus
Ntaus = prs.Ntaus;
PeakVel_cutoff = prs.PeakVel_cutoff;
Dist_cutoff = prs.DistPerc;
DistPerc = prs.DistPerc;
Dist2End = prs.Dist2End;
boots = prs.boots; % bootstrap 1 or not 0
ngroups = prs.ngroups; % accuracy magnitude groups
delta = prs.interoculardist/2;
zt = -prs.height;
saccade_duration = prs.saccade_duration(2); % min saccade duration
position_thresh = prs.position_thresh;
velthresh = prs.velthresh;
fixation_completed = prs.fixation_start;
fly_ONduration = prs.fly_ONduration;
Nboots = prs.bootstrap_trl;
factor_downsample = 1; %prs.factor_downsample; % downsampling factor for storage
ntrls = length(x_fly);
%% Extract data
for i=1:ntrls
    sacstart = []; sacend = []; sacampli = [];
    t_sac2 = t_sac{i};
    
    % ignore saccades right before the end of the trial
    if any(t_sac2 + saccade_duration >= ts{i}(end))
        ii = find(t_sac2 + saccade_duration/2 >= ts{i}(end));
        t_sac2(ii) = [];
    end
    % identify time of target fixation
    sac_indx = t_sac2>0 & t_sac2<fly_ONduration;
    if any(sac_indx)
        t_sacs = t_sac2(sac_indx);
        for j=1:length(t_sacs)
            sacstart(j) = find(ts{i}>(t_sacs(j)), 1);
            sacend(j) = find(ts{i}>(t_sacs(j) + saccade_duration), 1);
            sacampli(j) = nanmean([sum(abs(zle{i}(sacstart(j)) - zle{i}(sacend(j)))^2 + abs(yle{i}(sacstart(j)) - yle{i}(sacend(j)))^2) ...
                sum(abs(zre{i}(sacstart(j)) - zre{i}(sacend(j)))^2 + abs(yre{i}(sacstart(j)) - yre{i}(sacend(j)))^2)]);
        end
        t_fix(i) = t_sacs(sacampli == max(sacampli)) + saccade_duration;%  disp(['i = ' num2str(i) ', j = ' num2str(j)])
    else 
        t_fix(i) = 0 + fixation_completed; % if no saccade detected, assume monkey was already fixating on target
    end 
    
    % Save frame of target disappearance
    target_off(i) = find(ts{i} >= fly_ONduration, 1);

    % select data between target OFF and end of movement
    pretrial = 0; posttrial = 0;
    timeindx = find(ts{i}>=(fly_ONduration-pretrial) & ts{i}<(t_stop(i)+posttrial));
    
    % subject position
    d_sub{i} = d_sub{i}(timeindx);
    phi{i} = phi{i}(timeindx);
    w{i} = w{i}(timeindx);
    xmt{i} = x_monk{i}(timeindx) - x_monk{i}(find(~isnan(x_monk{i}),1)); 
    ymt{i} = y_monk{i}(timeindx) - y_monk{i}(find(~isnan(y_monk{i}),1));
    xmt{i}(isnan(xmt{i})) = xmt{i}(find(~isnan(xmt{i}),1));
    ymt{i}(isnan(ymt{i})) = ymt{i}(find(~isnan(ymt{i}),1));
    subj_r{i} = sqrt(xmt{i}.^2 + ymt{i}.^2);
    subj_th{i} = atan2d(xmt{i},ymt{i}); % [r_sub2eye{i},th_sub2eye{i}] = cart2polarY(xmt{i} - xmt{i}(end), ymt{i} - ymt{i}(end));
    
    % target position
    xt{i} = x_fly{i}(timeindx);
    yt{i} = y_fly{i}(timeindx);
    xt{i}(isnan(xt{i})) = xt{i}(find(~isnan(xt{i}),1));
    yt{i}(isnan(yt{i})) = yt{i}(find(~isnan(yt{i}),1));
    targ_r{i} = sqrt(xt{i}.^2 + yt{i}.^2);
    targ_th{i} = atan2d(xt{i},yt{i});

    % stopping position
    xs{i} = x_stop{i}(timeindx);
    ys{i} = y_stop{i}(timeindx);
    xs{i}(isnan(xs{i})) = xs{i}(find(~isnan(xs{i}),1));
    ys{i}(isnan(ys{i})) = ys{i}(find(~isnan(ys{i}),1));
    stop_r{i} = sqrt(xs{i}.^2 + ys{i}.^2);
    stop_th{i} = atan2d(xs{i},ys{i});

    % eye position
    yle{i} = yle{i}(timeindx); 
    yre{i} = yre{i}(timeindx);
    zle{i} = zle{i}(timeindx); 
    zre{i} = zre{i}(timeindx);
    
    % eye projection on arena
    [eye_r{i},eye_th{i},x_mean_obs{i},y_mean_obs{i}] = eye2world(zle{i},zre{i},yle{i},yre{i},zt);   
    [eye_r{i},eye_th{i},x_mean_obs{i},y_mean_obs{i}] = subject2arena_polar(eye_r{i},eye_th{i},subj_r{i},subj_th{i},phi{i}); 
        
    % prediction for eye position based on target location (if the monkey really followed the target)    
    y = yt{i}; y(y < 0) = nan;
    [yle_targ{i}, zle_targ{i}, yre_targ{i}, zre_targ{i}] = world2eye(xt{i},y,zt,delta);
%     [~,nanindx] = ReplaceWithNans([abs(zle_pred{i}) abs(zre_pred{i}) abs(yle_pred{i}) abs(yre_pred{i})], position_thresh, 0, 'normal');
%     yle_pred{i}(nanindx) = nan;    yre_pred{i}(nanindx) = nan;
%     zle_pred{i}(nanindx) = nan;    zre_pred{i}(nanindx) = nan;
    ver_mean_targ{i} = nanmean([zle_targ{i} , zre_targ{i}],2); 
    hor_mean_targ{i} = nanmean([yle_targ{i} , yre_targ{i}],2); 
    ver_diff_targ{i} = 0.5*(zle_targ{i} - zre_targ{i}); 
    hor_diff_targ{i} = 0.5*(yle_targ{i} - yre_targ{i}); 
        
    % prediction for eye position based on stopping location
    y = ys{i}; y(y < 0) = nan;
    [yle_stop{i}, zle_stop{i}, yre_stop{i}, zre_stop{i}] = world2eye(xs{i},y,zt,delta);
    ver_mean_stop{i} = nanmean([zle_stop{i} , zre_stop{i}],2); 
    hor_mean_stop{i} = nanmean([yle_stop{i} , yre_stop{i}],2); 
    
    % actual eye position
    ver_mean{i} = nanmean([zle{i} , zre{i}],2); % mean vertical eye position (of the two eyes)
    hor_mean{i} = nanmean([yle{i} , yre{i}],2); % mean horizontal eye position
    ver_diff{i} = 0.5*(zle{i} - zre{i}); % 0.5*difference between vertical eye positions (of the two eyes)
    hor_diff{i} = 0.5*(yle{i} - yre{i}); % 0.5*difference between horizontal eye positions
        
    if rmsac
    % remove saccades and concatenate eye position
    startpos.ver = ver_mean{i}(1);
    startpos.hor = hor_mean{i}(1);
    eyevel.ver = diff(ver_mean{i})/dt;
    eyevel.hor = diff(hor_mean{i})/dt;
    eyevel.mag = sqrt(eyevel.ver.^2 + eyevel.hor.^2);

    nanindx = unique([find(abs(eyevel.ver) > velthresh) ; find(abs(eyevel.hor) > velthresh) ; find(abs(eyevel.mag) > velthresh)]);
    eyevel.ver(nanindx) = nan;
    eyevel.hor(nanindx) = nan;
    pre_ver_mean = ver_mean{i};
    pre_hor_mean = hor_mean{i};
    ver_mean{i} = [0 ; startpos.ver + cumsum(eyevel.ver,'omitnan')*dt];
    hor_mean{i} = [0 ; startpos.hor + cumsum(eyevel.hor,'omitnan')*dt];
    end

    % align eye to target position at target offset
    if 0
    ver_mean{i} = ver_mean{i} + (ver_mean_stop{i}(1) - ver_mean{i}(1));
    hor_mean{i} = hor_mean{i} + (hor_mean_stop{i}(1) - hor_mean{i}(1));    
    end
        
    % saccade direction
    sacxy{i} = []; sacdir{i} = []; sacxy_targ{i} = []; sacdir_targ{i} = []; sacxy_stop{i} = []; sacdir_stop{i} = []; sac_time{i} = []; sac_trl{i} = [];
    if ~isempty(timeindx)
        for j=1:length(t_sac2)
            sacstartindx = find(ts{i}>(t_sac2(j) - saccade_duration/2), 1) - timeindx(1) - 2;
            sacendindx = find(ts{i}>(t_sac2(j) + saccade_duration/2), 1) - timeindx(1) + 2;
            if sacstartindx>0 && sacendindx<=(numel(timeindx)-30) % only consider saccades 500ms (= 30 samples) before stopping
                sacxy{i} = [sacxy{i} [ver_mean{i}(sacendindx) - ver_mean{i}(sacstartindx); 
                                      hor_mean{i}(sacendindx) - hor_mean{i}(sacstartindx)]];
                sacdir{i} = [sacdir{i} atan2d(sacxy{i}(1,end),sacxy{i}(2,end))];
                % prediction relative to target position
                sacxy_targ{i} = [sacxy_targ{i} [ver_mean_targ{i}(sacstartindx) - ver_mean{i}(sacstartindx); 
                                                hor_mean_targ{i}(sacstartindx) - hor_mean{i}(sacstartindx)]];
                sacdir_targ{i} = [sacdir_targ{i} atan2d(sacxy_targ{i}(1,end),sacxy_targ{i}(2,end))];
                % prediction relative to stopping position
                sacxy_stop{i} = [sacxy_stop{i} [ver_mean_targ{i}(sacstartindx) - ver_mean{i}(sacstartindx);
                                                hor_mean_targ{i}(sacstartindx) - hor_mean{i}(sacstartindx)]];
                sacdir_stop{i} = [sacdir_stop{i} atan2d(sacxy_stop{i}(1,end),sacxy_stop{i}(2,end))];
                % saccade timing
                sac_time{i} = [sac_time{i} sacstartindx*dt]; % time since target off
                sac_trl{i} = [sac_trl{i} i];
            end
        end
    end
end

%% save saccadic eye movements
tracking.saccade.true.val = cell2mat(sacxy);
tracking.saccade.true.dir = cell2mat(sacdir);
tracking.saccade.pred_targ.val = cell2mat(sacxy_targ);
tracking.saccade.pred_targ.dir = cell2mat(sacdir_targ);
tracking.saccade.pred_stop.val = cell2mat(sacxy_stop);
tracking.saccade.pred_stop.dir = cell2mat(sacdir_stop);
tracking.saccade.time = cell2mat(sac_time);
tracking.saccade.trl = cell2mat(sac_trl);

%% Target and Stop position tracking errors
Nt = max(cellfun(@(x) length(x),xt)); % max number of timepoints

for i=1:ntrls
    nt = length(ver_mean{i}); 
    t_cut = nt; %     distind = find(d_sub{i} >= DistPerc*max(d_sub{i}),1); t_cut = distind; 

    tar_sign(i) = sign(hor_mean_targ{i}(1));
    % Target-Tracking Error
    TTE_ver(i,:) = [ver_mean_targ{i}(1:t_cut) - ver_mean{i}(1:t_cut) ;  nan(Nt - t_cut,1)];
    TTE_hor(i,:) = [hor_mean_targ{i}(1:t_cut) - hor_mean{i}(1:t_cut) ;  nan(Nt - t_cut,1)];
    TTE_all(i,:) = sqrt(TTE_ver(i,:).^2 + TTE_hor(i,:).^2);
    
    % Stop-Position-Tracking Error
    SPTE_ver(i,:) = [ver_mean_stop{i}(1:t_cut) - ver_mean{i}(1:t_cut) ;  nan(Nt - t_cut,1)];
    SPTE_hor(i,:) = [hor_mean_stop{i}(1:t_cut) - hor_mean{i}(1:t_cut) ;  nan(Nt - t_cut,1)];
    SPTE_all(i,:) = sqrt(SPTE_ver(i,:).^2 + SPTE_hor(i,:).^2);
    
    % Trial Errors (in screen coordinates)
    trlerrors.screen.ver(i,:) = [ver_mean_targ{i}(1:t_cut) - ver_mean_stop{i}(1:t_cut) ; nan(Nt - t_cut,1)];
    trlerrors.screen.hor(i,:) = [hor_mean_targ{i}(1:t_cut) - hor_mean_stop{i}(1:t_cut) ; nan(Nt - t_cut,1)];
    trlerrors.screen.all(i,:) = sqrt(trlerrors.screen.ver(i,:).^2 + trlerrors.screen.hor(i,:).^2);
    trlerrors.arena.all(i,:) = behverrors(i)*ones(1,Nt);
    
end

ts = [1:Nt]*dt;
% correlation of tracking error and steering error in screen coordinates
[tracking.eyepos.errcorr.targ.screen.r,tracking.eyepos.errcorr.targ.screen.p] = ...
       arrayfun(@(t) nancorr(trlerrors.screen.all(:,t),TTE_all(:,t)), 1:Nt);
[tracking.eyepos.errcorr.stop.screen.r,tracking.eyepos.errcorr.stop.screen.p] = ...
       arrayfun(@(t) nancorr(trlerrors.screen.all(:,t),SPTE_all(:,t)), 1:Nt);
% correlation of tracking error and steering error in arena coordinates
[tracking.eyepos.errcorr.targ.arena.r,tracking.eyepos.errcorr.targ.arena.p] = ...
       arrayfun(@(t) nancorr(trlerrors.arena.all(:,t),TTE_all(:,t)), 1:Nt);
[tracking.eyepos.errcorr.stop.arena.r,tracking.eyepos.errcorr.stop.arena.p] = ...
       arrayfun(@(t) nancorr(trlerrors.arena.all(:,t),SPTE_all(:,t)), 1:Nt);

tracking.eyepos.error.targ.hor = TTE_hor;
tracking.eyepos.error.targ.ver = TTE_ver;
tracking.eyepos.error.targ.all = TTE_all;
tracking.eyepos.error.stop.hor = SPTE_hor;
tracking.eyepos.error.stop.ver = SPTE_ver;
tracking.eyepos.error.stop.all = SPTE_all;

tracking.eyepos.error.behv_arena.all = trlerrors.arena.all;
tracking.eyepos.error.behv_screen.ver = trlerrors.screen.ver;
tracking.eyepos.error.behv_screen.hor = trlerrors.screen.hor;    

%% save true and predicted eye positions
Nt = max(cellfun(@(x) length(x),xt)); % max number of timepoints

tracking.tarpos.arena.r = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],targ_r,'uniformoutput',false)),factor_downsample)';
tracking.tarpos.arena.theta = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],targ_th,'uniformoutput',false)),factor_downsample)';
tracking.tarpos.arena.x = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],xt,'uniformoutput',false)),factor_downsample)';
tracking.tarpos.arena.y = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],yt,'uniformoutput',false)),factor_downsample)';
tracking.tarpos.screen.ver_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_targ,'uniformoutput',false)),factor_downsample)';
tracking.tarpos.screen.hor_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_targ,'uniformoutput',false)),factor_downsample)';
% tracking.tarpos.screen.ver_diff = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff_targ,'uniformoutput',false)),factor_downsample)';
% tracking.tarpos.screen.hor_diff = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff_targ,'uniformoutput',false)),factor_downsample)';

tracking.stopos.arena.r = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],stop_r,'uniformoutput',false)),factor_downsample)';
tracking.stopos.arena.theta = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],stop_th,'uniformoutput',false)),factor_downsample)';
tracking.stopos.arena.x = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],xs,'uniformoutput',false)),factor_downsample)';
tracking.stopos.arena.y = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ys,'uniformoutput',false)),factor_downsample)';
tracking.stopos.screen.ver_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_stop,'uniformoutput',false)),factor_downsample)';
tracking.stopos.screen.hor_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_stop,'uniformoutput',false)),factor_downsample)';

tracking.eyepos.arena.r = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],eye_r,'uniformoutput',false)),factor_downsample)';
tracking.eyepos.arena.theta = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],eye_th,'uniformoutput',false)),factor_downsample)';
tracking.eyepos.arena.x = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],x_mean_obs,'uniformoutput',false)),factor_downsample)';
tracking.eyepos.arena.y = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],y_mean_obs,'uniformoutput',false)),factor_downsample)';
tracking.eyepos.screen.ver_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean,'uniformoutput',false)),factor_downsample)';
tracking.eyepos.screen.hor_mean = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean,'uniformoutput',false)),factor_downsample)';
% tracking.eyepos.screen.ver_diff = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff,'uniformoutput',false)),factor_downsample)';
% tracking.eyepos.screen.hor_diff = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff,'uniformoutput',false)),factor_downsample)';

tracking.misc.ts = ts(:);
tracking.misc.tau.val = tau(:);
tracking.misc.target_off.indx = target_off(:);
tracking.misc.endofsteering.val = endofsteering(:);
tracking.misc.behverrors.all.val = behverrors(:);
tracking.misc.behverrors.r.val = cellfun(@(xt,xs) xt(1)-xs(end), targ_r, subj_r);
tracking.misc.behverrors.th.val = cellfun(@(xt,xs) xt(1)-xs(end), targ_th, subj_th);
tracking.misc.subdist = downsample(cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],d_sub,'uniformoutput',false)),factor_downsample)';
tracking.misc.sampleflag = 'time';
tracking.misc.rmsac = logical(rmsac);

Nt = max(cellfun(@(x) length(x),xt)); % max number of timepoints
%% compare predicted & eye positions from all trials aligned to target fixation
% target position
ver_mean_targ1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_targ,'UniformOutput',false));
hor_mean_targ1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_targ,'UniformOutput',false));
% stoping position
ver_mean_stop1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_stop,'UniformOutput',false));
hor_mean_stop1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_stop,'UniformOutput',false));
% true eye position
ver_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean,'UniformOutput',false));
hor_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean,'UniformOutput',false));

% cross-correlation between observed and predicted eye positions
zeropad_length = 200; % pad these many bins around each trial
maxlag = 400;
% vertical
obs = [zeros(zeropad_length,ntrls); ver_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;

targ = [zeros(zeropad_length,ntrls); ver_mean_targ1; zeros(zeropad_length,ntrls)]; 
targ_shuffled = targ(:,randperm(ntrls));
targ = targ(:); targ(isnan(targ)) = 0; targ_shuffled = targ_shuffled(:); targ_shuffled(isnan(targ_shuffled)) = 0;
[tracking.eyepos.pred_vs_true.targ.ver_mean.xcorr.c, tracking.eyepos.pred_vs_true.targ.ver_mean.xcorr.lag] = xcorr(targ,obs,maxlag,'coeff');
tracking.eyepos.pred_vs_true.targ.ver_mean.xcorr.c_shuffled = xcorr(targ_shuffled,obs,maxlag,'coeff');

stop = [zeros(zeropad_length,ntrls); ver_mean_stop1; zeros(zeropad_length,ntrls)]; 
stop_shuffled = stop(:,randperm(ntrls));
stop = stop(:); stop(isnan(stop)) = 0; stop_shuffled = stop_shuffled(:); stop_shuffled(isnan(stop_shuffled)) = 0;
[tracking.eyepos.pred_vs_true.stop.ver_mean.xcorr.c, tracking.eyepos.pred_vs_true.stop.ver_mean.xcorr.lag] = xcorr(stop,obs,maxlag,'coeff');
tracking.eyepos.pred_vs_true.stop.ver_mean.xcorr.c_shuffled = xcorr(stop_shuffled,obs,maxlag,'coeff');
% horizontal
obs = [zeros(zeropad_length,ntrls); hor_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;

targ = [zeros(zeropad_length,ntrls); hor_mean_targ1; zeros(zeropad_length,ntrls)]; 
targ_shuffled = targ(:,randperm(ntrls));
targ = targ(:); targ(isnan(targ)) = 0; targ_shuffled = targ_shuffled(:); targ_shuffled(isnan(targ_shuffled)) = 0;
[tracking.eyepos.pred_vs_true.targ.hor_mean.xcorr.c, tracking.eyepos.pred_vs_true.targ.hor_mean.xcorr.lag] = xcorr(targ,obs,maxlag,'coeff');
tracking.eyepos.pred_vs_true.targ.hor_mean.xcorr.c_shuffled = xcorr(targ_shuffled,obs,maxlag,'coeff');

stop = [zeros(zeropad_length,ntrls); hor_mean_stop1; zeros(zeropad_length,ntrls)]; 
stop_shuffled = stop(:,randperm(ntrls));
stop = stop(:); stop(isnan(stop)) = 0; stop_shuffled = stop_shuffled(:); stop_shuffled(isnan(stop_shuffled)) = 0;
[tracking.eyepos.pred_vs_true.stop.hor_mean.xcorr.c, tracking.eyepos.pred_vs_true.stop.hor_mean.xcorr.lag] = xcorr(stop,obs,maxlag,'coeff');
tracking.eyepos.pred_vs_true.stop.hor_mean.xcorr.c_shuffled = xcorr(stop_shuffled,obs,maxlag,'coeff');

% timecourse of component-wise corr between predicted & true eye position
[tracking.eyepos.pred_vs_true.targ.ver_mean.rho.startaligned,tracking.eyepos.pred_vs_true.targ.ver_mean.pval.startaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_targ1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.targ.hor_mean.rho.startaligned,tracking.eyepos.pred_vs_true.targ.hor_mean.pval.startaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_targ1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.stop.ver_mean.rho.startaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.pval.startaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_stop1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.stop.hor_mean.rho.startaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.pval.startaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_stop1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% timecourse of component-wise regression between predicted & true eye position
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.ver_mean.beta.startaligned,tracking.eyepos.pred_vs_true.targ.ver_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.hor_mean.beta.startaligned,tracking.eyepos.pred_vs_true.targ.hor_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.stop.ver_mean.beta.startaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.stop.hor_mean.beta.startaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
% timecourse of component-wise multiple regression between target/stop & true eye position
multi = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ1(i,:)' ver_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.ver_mean.multi.startaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.multi.startaligned] = deal(multi(1,:),multi(2,:));
multi = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ1(i,:)' hor_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.hor_mean.multi.startaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.multi.startaligned] = deal(multi(1,:),multi(2,:));

% timecourse of variance explained
targ = permute(cat(3,ver_mean_targ1 , hor_mean_targ1),[3 1 2]);
stop = permute(cat(3,ver_mean_stop1 , hor_mean_stop1),[3 1 2]);
obs = permute(cat(3,ver_mean1 , hor_mean1),[3 1 2]);
if boots
    var_explained_targ = nan(Nboots,Nt);    var_explained_targ_shuffled = nan(Nboots,Nt);
    var_explained_stop = nan(Nboots,Nt);    var_explained_stop_shuffled = nan(Nboots,Nt);
    for i=1:Nboots
        randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
        % target position
        squared_err = sum(nanmean((obs(:,:,randtrls) - targ(:,:,randtrls)).^2,3));
        var_pred = sum(nanvar(targ(:,:,randtrls),[],3));
        var_explained_targ(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - targ(:,:,randtrls2)).^2,3));
        var_pred_shuffled = sum(nanvar(targ(:,:,randtrls2),[],3));
        var_explained_targ_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
        % stopping position
        squared_err = sum(nanmean((obs(:,:,randtrls) - stop(:,:,randtrls)).^2,3));
        var_pred = sum(nanvar(stop(:,:,randtrls),[],3));
        var_explained_stop(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - stop(:,:,randtrls2)).^2,3));
        var_pred_shuffled = sum(nanvar(stop(:,:,randtrls2),[],3));
        var_explained_stop_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
        
    end
    tracking.eyepos.pred_vs_true.targ.var_explained.mu.startaligned = nanmean(var_explained_targ);
    tracking.eyepos.pred_vs_true.targ.var_explained.sem.startaligned = nanstd(var_explained_targ);
    tracking.eyepos.pred_vs_true.targ.var_explained_shuffled.mu.startaligned = nanmean(var_explained_targ_shuffled);
    tracking.eyepos.pred_vs_true.targ.var_explained_shuffled.sem.startaligned = nanstd(var_explained_targ_shuffled);

    tracking.eyepos.pred_vs_true.stop.var_explained.mu.startaligned = nanmean(var_explained_stop);
    tracking.eyepos.pred_vs_true.stop.var_explained.sem.startaligned = nanstd(var_explained_stop);
    tracking.eyepos.pred_vs_true.stop.var_explained_shuffled.mu.startaligned = nanmean(var_explained_stop_shuffled);
    tracking.eyepos.pred_vs_true.stop.var_explained_shuffled.sem.startaligned = nanstd(var_explained_stop_shuffled);

else
    % target position
    squared_err = sum(nanmean((obs - targ).^2,3));
    var_pred = sum(nanvar(targ,[],3));
    var_explained_targ = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking.eyepos.pred_vs_true.targ.var_explained.mu.startaligned = var_explained_targ;
    % stopping position
    squared_err = sum(nanmean((obs - stop).^2,3));
    var_pred = sum(nanvar(stop,[],3));
    var_explained_stop = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking.eyepos.pred_vs_true.stop.var_explained.mu.startaligned = var_explained_stop;
end
%% compare predicted & observed eye positions from all trials aligned to END of movement
% disregard eye movements towards the end of the trial
cutoff = cellfun(@(x) find(x >= DistPerc*max(x),1), d_sub,'uniformoutput',false); % cutoff = cellfun(@(x) find(x >= max(x)-Dist2End,1), d_sub,'uniformoutput',false); 
% target position
ver_mean_targ2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], ver_mean_targ, cutoff, 'UniformOutput',false));
hor_mean_targ2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], hor_mean_targ, cutoff, 'UniformOutput',false));
% stopping position
ver_mean_stop2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], ver_mean_stop, cutoff, 'UniformOutput',false));
hor_mean_stop2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], hor_mean_stop, cutoff, 'UniformOutput',false));
% true eye position
ver_mean2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], ver_mean, cutoff, 'UniformOutput',false));
hor_mean2 = cell2mat(cellfun(@(x,ind) [flipud(x(1:ind)) ; nan(Nt - ind,1)], hor_mean, cutoff, 'UniformOutput',false));

% timecourse of component-wise corr between predicted & true eye position
[tracking.eyepos.pred_vs_true.targ.ver_mean.rho.stopaligned,tracking.eyepos.pred_vs_true.targ.ver_mean.pval.stopaligned] = arrayfun(@(i) corr(ver_mean2(i,:)',ver_mean_targ2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.targ.hor_mean.rho.stopaligned,tracking.eyepos.pred_vs_true.targ.hor_mean.pval.stopaligned] = arrayfun(@(i) corr(hor_mean2(i,:)',hor_mean_targ2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.stop.ver_mean.rho.stopaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.pval.stopaligned] = arrayfun(@(i) corr(ver_mean2(i,:)',ver_mean_stop2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking.eyepos.pred_vs_true.stop.hor_mean.rho.stopaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.pval.stopaligned] = arrayfun(@(i) corr(hor_mean2(i,:)',hor_mean_stop2(i,:)','Type','Spearman','rows','complete'), 1:Nt);

% component-wise regression between predicted & true eye position
beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_targ2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.ver_mean.beta.stopaligned,tracking.eyepos.pred_vs_true.targ.ver_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_targ2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.hor_mean.beta.stopaligned,tracking.eyepos.pred_vs_true.targ.hor_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.stop.ver_mean.beta.stopaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.stop.hor_mean.beta.stopaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));

% timecourse of component-wise multiple regression between target/stop & true eye position
multi = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ2(i,:)' ver_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.ver_mean.multi.stopaligned,tracking.eyepos.pred_vs_true.stop.ver_mean.multi.stopaligned] = deal(multi(1,:),multi(2,:));
multi = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ2(i,:)' hor_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking.eyepos.pred_vs_true.targ.hor_mean.multi.stopaligned,tracking.eyepos.pred_vs_true.stop.hor_mean.multi.stopaligned] = deal(multi(1,:),multi(2,:));

% timecourse of variance explained
targ = permute(cat(3,ver_mean_targ2 , hor_mean_targ2),[3 1 2]);
stop = permute(cat(3,ver_mean_stop2 , hor_mean_stop2),[3 1 2]);
obs = permute(cat(3,ver_mean2 , hor_mean2),[3 1 2]);
if boots
    var_explained_targ = nan(Nboots,Nt); var_explained_targ_shuffled = nan(Nboots,Nt);
    var_explained_stop = nan(Nboots,Nt); var_explained_stop_shuffled = nan(Nboots,Nt);
    for i=1:Nboots
        randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
        % target position
        squared_err = sum(nanmean((obs(:,:,randtrls) - targ(:,:,randtrls)).^2,3)); var_pred = sum(nanvar(targ(:,:,randtrls),[],3));
        var_explained_targ(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - targ(:,:,randtrls2)).^2,3)); var_pred_shuffled = sum(nanvar(targ(:,:,randtrls2),[],3));
        var_explained_targ_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
        % stopping location
        squared_err = sum(nanmean((obs(:,:,randtrls) - stop(:,:,randtrls)).^2,3)); var_pred = sum(nanvar(stop(:,:,randtrls),[],3));
        var_explained_stop(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - stop(:,:,randtrls2)).^2,3)); var_pred_shuffled = sum(nanvar(stop(:,:,randtrls2),[],3));
        var_explained_stop_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
    end
    tracking.eyepos.pred_vs_true.targ.var_explained.mu.stopaligned = nanmean(var_explained_targ);
    tracking.eyepos.pred_vs_true.targ.var_explained.sem.stopaligned = nanstd(var_explained_targ);
    tracking.eyepos.pred_vs_true.targ.var_explained_shuffled.mu.stopaligned = nanmean(var_explained_targ_shuffled);
    tracking.eyepos.pred_vs_true.targ.var_explained_shuffled.sem.stopaligned = nanstd(var_explained_targ_shuffled);

    tracking.eyepos.pred_vs_true.stop.var_explained.mu.stopaligned = nanmean(var_explained_stop);
    tracking.eyepos.pred_vs_true.stop.var_explained.sem.stopaligned = nanstd(var_explained_stop);
    tracking.eyepos.pred_vs_true.stop.var_explained_shuffled.mu.stopaligned = nanmean(var_explained_stop_shuffled);
    tracking.eyepos.pred_vs_true.stop.var_explained_shuffled.sem.stopaligned = nanstd(var_explained_stop_shuffled);
else
    squared_err = sum(nanmean((obs - targ).^2,3));
    var_pred = sum(nanvar(targ,[],3));
    var_explained_targ = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking.eyepos.pred_vs_true.targ.var_explained.mu.stopaligned = var_explained_targ;

    squared_err = sum(nanmean((obs - stop).^2,3));
    var_pred = sum(nanvar(stop,[],3));
    var_explained_stop = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking.eyepos.pred_vs_true.stop.var_explained.mu.stopaligned = var_explained_stop;

end
