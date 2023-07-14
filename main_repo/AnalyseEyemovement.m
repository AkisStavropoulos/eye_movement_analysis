function eye_movement = AnalyseEyemovement(eye_fixation,trials,spatialerr,prs,BELIEF)
% BELIEF: 1 to use subject's believed position
%         0 to use actual position

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
% extract bias
if BELIEF
    b0 = [1.2 1.2];
    [b,fval] = fminsearch(@(b) scaleV(trials,b),b0);
    b_w = b(1);  b_v = b(2);
%     b_x = regress(x_sub(:),x_tar(:));                  b_y = regress(y_sub(:),y_tar(:));
%     a_x = regress(x_sub(:),[x_tar(:) x_tar(:).*tau(:)]);  a_y = regress(y_sub(:),[y_tar(:) y_tar(:).*tau(:)]);
else
    b_w = 1;    b_v = 1;
end

% set up data
for j = 1:length(trials)
    tau(j) = trials(j).prs.tau;
    w{j} = trials(j).continuous.w(:);
    v{j} = trials(j).continuous.v(:);   
%     x_monk{j} = trials(j).continuous.xmp;
%     y_monk{j} = trials(j).continuous.ymp;
    [x_monk{j}, y_monk{j}, ~, ~, ~, phi{j}] = gen_traj(w{j}, v{j}, trials(j).continuous.xmp(1), trials(j).continuous.ymp(1),[]);
    x_monk{j} = x_monk{j}(1:end-1);
    y_monk{j} = y_monk{j}(1:end-1);
    phi{j} = phi{j}(1:end-1);
%     b_x = a_x(1) + a_x(2)*tau(j);
%     b_y = a_y(1) + a_y(2)*tau(j);
    blv_v = b_v*v{j};
    blv_w = b_w*w{j};
    [blv_x, blv_y, ~, ~, ~, blv_phi] = gen_traj(blv_w, blv_v, trials(j).continuous.xmp(1), trials(j).continuous.ymp(1),[]);
    blv_x = blv_x(1:end-1);
    blv_y = blv_y(1:end-1);
    blv_phi = blv_phi(1:end-1);
    
    x_fly{j} = trials(j).prs.fireflyposx - blv_x; 
    y_fly{j} = trials(j).prs.fireflyposy - blv_y;
    x_stop{j} = blv_x(end) - blv_x;
    y_stop{j} = blv_y(end) - blv_y;
    x_monk{j} = blv_x;
    y_monk{j} = blv_y;

    if length(blv_phi) ~= length(x_fly{j}); blv_phi = [blv_phi ; blv_phi(end)*ones(length(length(x_fly{j}) - length(blv_phi)),1)]; end
    phi{j} = cumsum(blv_w)*dt;
    A = []; B = [];
    for t = 1:length(x_fly{j})
    R = [cosd(phi{j}(t))   sind(phi{j}(t)) ;...
         -sind(phi{j}(t))  cosd(phi{j}(t))];
     A(t,:) = [x_fly{j}(t) y_fly{j}(t)]*R;
     B(t,:) = [x_stop{j}(t) y_stop{j}(t)]*R;
    end
    if 0
        plot(x_monk{j},y_monk{j},'k'); hold on; plot(trials(j).prs.fireflyposx,trials(j).prs.fireflyposy,'rx');
        plot(x_fly{j},y_fly{j},'r'); plot(A(:,1),A(:,2),'b'); plot(A(end,1),A(end,2),'bx');
        xlim([-500 500]); ylim([-50 950]); hold off; title(['Trial No. ' num2str(j)])
    end

    x_fly{j} = A(:,1);    y_fly{j} = A(:,2);
    x_stop{j} = B(:,1);   y_stop{j} = B(:,2);

    d_sub{j} = cumsum(trials(j).continuous.v)*(1/60);
    zle{j} = trials(j).continuous.zle;    %zle{j}(trials(j).continuous.blink) = nan;
    yle{j} = trials(j).continuous.yle;    %yle{j}(trials(j).continuous.blink) = nan;
    zre{j} = trials(j).continuous.zre;    %zre{j}(trials(j).continuous.blink) = nan;
    yre{j} = trials(j).continuous.yre;    %yre{j}(trials(j).continuous.blink) = nan;
    t_sac{j} = trials(j).events.t_sac - trials(j).events.t_beg;
    t_stop(j) = trials(j).events.t_end - trials(j).events.t_beg;
    ts{j} = trials(j).continuous.ts;
    endofsteering(j) = find(trials(j).mc.JS_X_Raw,1,'last');
    behverrors(j) = sqrt(x_fly{j}(end).^2 + y_fly{j}(end).^2);
end
% figure; hold on;    plot(x_monk{j},y_monk{j},'k'); plot(x_tar(j),y_tar(j),'kx') 
%                     plot(x_monk{j}./b_x , y_monk{j}./b_y,'b'); plot(x_tar(j),y_tar(j),'bx') 
%                     plot(x_monk{j},y_monk{j},'r'); plot(x_tar(j)*b_x , y_tar(j)*b_y,'rx') 

% sanity check
if 0
for j = 1:length(trials)
    err_x(j) = x_fly{j}(end);
    err_y(j) = y_fly{j}(end);
end
figure;plot(-err_x,-err_y,'.')
vline(0,'w');hline(0,'w');

figure; hold on;  plot(x_tar,x_sub,'.');  plot(x_tar,x_tar*b_x,'r');   plot(-500:500,-500:500,'k--');   title('X');   xlim([-500 500]);  ylim([-500 500]);
figure; hold on;  plot(y_tar,y_sub,'.');  plot(y_tar,y_tar*b_y,'r');   plot(-500:500,-500:500,'k--');   title('Y');   xlim([0 500]);     ylim([0 500]);

figure; scatter3(x_tar,y_tar,x_tar-x_sub); xlabel('x_{tar}'); ylabel('y_{tar}');  zlabel('x_{sub}')
figure; scatter3(x_tar,y_tar,y_tar-y_sub); xlabel('x_{tar}'); ylabel('y_{tar}');  zlabel('y_{sub}')
end
%% prs
TauAdjust = prs.TauAdjust; % adjust frame cutoffs to varying duration of stopping time caused by varying taus
Ntaus = prs.Ntaus;
PeakVel_cutoff = prs.PeakVel_cutoff;
Dist_cutoff = prs.Dist_cutoff;
DistPerc = prs.DistPerc;
boots = prs.boots; % bootstrap 1 or not 0
norm_TTI = prs.norm; % normalize TTI by initial value of TTI
ngroups = prs.ngroups; % accuracy magnitude groups
delta = prs.interoculardist/2;
zt = -prs.height;
saccade_duration = prs.saccade_duration(2); % min saccade duration
position_thresh = prs.position_thresh;
fixation_completed = prs.fixation_start;
fly_ONduration = prs.fly_ONduration;
Nboots = prs.bootstrap_trl;
factor_downsample = 1; %prs.factor_downsample; % downsampling factor for storage
ntrls = length(x_fly);
pretrial = prs.pretrial;
posttrial = prs.posttrial;
if TauAdjust
    cutoffType = 'Tau';
elseif PeakVel_cutoff
    cutoffType = 'PeakVel';
elseif Dist_cutoff
    cutoffType = [num2str(DistPerc*100) '% Distance'];
else
    cutoffType = 'end of trial';
end
eye_movement.prs.CutOffType = cutoffType;
eye_movement.prs.bootstrap = boots;
eye_movement.prs.AccGroups = ngroups;
eye_movement.prs.ActualError = ~BELIEF;
%% sort trials by euclideian error
[~,errorindx] = sort(behverrors);

%% eye position immediately after the first saccade following target onset
for i=1:ntrls
%     figure;plot(ts{i},[zle{i} yle{i} zre{i} yre{i}])
%     vline(t_sac{i});title(['Trial No. ' num2str(i)])
    sacstart = []; sacend = []; sacampli = [];
    t_sac2 = t_sac{i};
    
    % ignore saccades right before the end of the trial
    if any(t_sac2 + saccade_duration/2 >= ts{i}(end))
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
    % remove saccade periods from eye position data
%     sacstart = []; sacend = [];
%     for j=1:length(t_sac2)
%         sacstart(j) = find(ts{i}>(t_sac2(j) - saccade_duration/2), 1);
%         if isempty(find(ts{i}>(t_sac2(j) + saccade_duration/2), 1))
%             sacstart(j) = [];   disp(['Trial ' num2str(i)])
%         else
%             sacend(j) = find(ts{i}>(t_sac2(j) + saccade_duration/2), 1);
%         end
%         yle{i}(sacstart(j):sacend(j)) = nan; % left eye horizontal position
%         yre{i}(sacstart(j):sacend(j)) = nan; % right eye horizontal position
%         zle{i}(sacstart(j):sacend(j)) = nan; % left eye vertical position
%         zre{i}(sacstart(j):sacend(j)) = nan; % right eye vertical
%         position 
%     end
%     t_fix(i) = 0;
    
    % Save frame of target disappearance
    target_off(i) = find(ts{i} >= 1, 1);

    % select data between target OFF and end of movement
    pretrial = 0; posttrial = 0;
    timeindx = find(ts{i}>=(0-pretrial) & ts{i}<(t_stop(i)+posttrial));
    
    
    % subject position
    d_sub{i} = d_sub{i}(timeindx);
    phi{i} = phi{i}(timeindx);
    w{i} = w{i}(timeindx);
    xmt{i} = x_monk{i}(timeindx) - x_monk{i}(find(~isnan(x_monk{i}),1)); 
    ymt{i} = y_monk{i}(timeindx) - y_monk{i}(find(~isnan(y_monk{i}),1));
    xmt{i}(isnan(xmt{i})) = xmt{i}(find(~isnan(xmt{i}),1));
    ymt{i}(isnan(ymt{i})) = ymt{i}(find(~isnan(ymt{i}),1));

    % subject r and theta (in arena coordinates)
    sub_r{i} = sqrt(xmt{i}.^2 + ymt{i}.^2);
    sub_th{i} = atan2d(xmt{i},ymt{i});
    [r_sub2eye{i},th_sub2eye{i}] = cart2polarY(xmt{i} - xmt{i}(end), ymt{i} - ymt{i}(end));

    % variance in subject position
    if ~isempty(spatialerr)
    [~,x_indx] = min(abs(xmt{i}-spatialerr.x),[],2);
    [~,y_indx] = min(abs(ymt{i}-spatialerr.y),[],2);
    var_xmt{i} = (spatialerr.x_err_smooth(sub2ind(size(spatialerr.x_err_smooth), x_indx, y_indx))).^2;
    var_ymt{i} = (spatialerr.y_err_smooth(sub2ind(size(spatialerr.y_err_smooth), x_indx, y_indx))).^2;
    end
    
    % target position
    xt{i} = x_fly{i}(timeindx);
    yt{i} = y_fly{i}(timeindx);
    xt{i}(isnan(xt{i})) = xt{i}(find(~isnan(xt{i}),1));
    yt{i}(isnan(yt{i})) = yt{i}(find(~isnan(yt{i}),1));
    
    % subject position
    xs{i} = x_stop{i}(timeindx);
    ys{i} = y_stop{i}(timeindx);
    xs{i}(isnan(xs{i})) = xs{i}(find(~isnan(xs{i}),1));
    ys{i}(isnan(ys{i})) = ys{i}(find(~isnan(ys{i}),1));

    % target r and theta (in subject coordinates)
    tar_r_sub{i} = sqrt(xt{i}.^2 + yt{i}.^2);
    tar_th_sub{i} = atan2d(xt{i},yt{i});

    % eye position
    yle{i} = yle{i}(timeindx); 
    yre{i} = yre{i}(timeindx);
    zle{i} = zle{i}(timeindx); 
    zre{i} = zre{i}(timeindx);
    
    % eye projection on arena
    if 0
    [x_mean_obs{i},y_mean_obs{i}] = eye2subject(yle{i},zle{i},yre{i},zre{i},delta,zt);
    [x_mean_obs{i},y_mean_obs{i}] = subject2arena(x_mean_obs{i},y_mean_obs{i},xmt{i},ymt{i},phi{i});
    eye_r_arena{i} = sqrt((x_mean_obs{i} - xmt{i}(1)).^2 + (y_mean_obs{i} - ymt{i}(1)).^2);
    eye_th_arena{i} = atan2d(x_mean_obs{i} - xmt{i}(1),y_mean_obs{i} - ymt{i}(1));
    else
    [eye_r{i},eye_th{i},x_mean_obs{i},y_mean_obs{i}] = eyepos2flypos(zle{i},zre{i},yle{i},yre{i},zt);   
    [eye_r_arena{i},eye_th_arena{i},x_mean_obs{i},y_mean_obs{i}] = subject2arena_polar(eye_r{i},eye_th{i},sub_r{i},sub_th{i},phi{i}); 
    end
    
    if 0
        colr = jet(length(x_mean_obs{i}));
        for n = 1:length(x_mean_obs{i})
        plot(x_mean_obs{i}(n),y_mean_obs{i}(n),'.','color',colr(n,:)); hold on;
        end
        hold on; plot(xt{i}(1),yt{i}(1),'rx'); plot(x_mean_obs{i}(1),y_mean_obs{i}(1),'go');
        A = [];
        for t = 1:length(xt{i})
            R = [cosd(phi{i}(t))   -sind(phi{i}(t)) ;...
                sind(phi{i}(t))  cosd(phi{i}(t))];
            A(t,:) = [xt{i}(t) yt{i}(t)]*R + [xmt{i}(t) ymt{i}(t)];
        end
        xtemp = A(:,1);
        ytemp = A(:,2);

        plot(xtemp,ytemp,'ro'); xlim([-500 500]); ylim([0 1000]); hold off; title(['Trial No. ' num2str(i)])

% %         plot(xLE{i},yLE{i},'b'); hold on; plot(xRE{i},yRE{i},'r'); hold off;
%         plot(xt{i}(1),yt{i}(1),'rx'); hold on; plot(x_mean_obs{i}(1),y_mean_obs{i}(1),'bo'); % initial positions
%         plot(xt{i},yt{i},'r-'); plot(x_mean_obs{i},y_mean_obs{i},'b-'); % trajectories
    
        subplot(2,1,1); plot(ts{i}(timeindx),eye_r_arena{i},'k'); hold on; xlim([0 8]);
        subplot(2,1,2); plot(ts{i}(timeindx),eye_th_arena{i},'k'); hold on; xlim([0 8]);
    end
    
    % ground truth prediction for eye position (if the monkey really followed the target)
    y0 = atan2d(0, sqrt(0.^2 + zt^2));
    z0 = atan2d(zt , sqrt(0.^2 + 0.^2));
    
    y = yt{i}; y(y < 0) = nan;
    yle_pred{i} = atan2d(xt{i} + delta, sqrt(y.^2 + zt^2));
    yre_pred{i} = atan2d(xt{i} - delta, sqrt(y.^2 + zt^2));
    zle_pred{i} = atan2d(zt , sqrt(y.^2 + (xt{i} + delta).^2));
    zre_pred{i} = atan2d(zt , sqrt(y.^2 + (xt{i} - delta).^2));
    [~,nanindx] = ReplaceWithNans([abs(zle_pred{i}) abs(zre_pred{i}) abs(yle_pred{i}) abs(yre_pred{i})], position_thresh, 0, 'normal');
    yle_pred{i}(nanindx) = nan;    yre_pred{i}(nanindx) = nan;
    zle_pred{i}(nanindx) = nan;    zre_pred{i}(nanindx) = nan;

    ver_mean_pred{i} = nanmean([zle_pred{i} , zre_pred{i}],2); % mean vertical eye position (of the two eyes)    
    hor_mean_pred{i} = nanmean([yle_pred{i} , yre_pred{i}],2); % mean horizontal eye position
    ver_diff_pred{i} = 0.5*(zle_pred{i} - zre_pred{i}); % 0.5*difference between vertical eye positions (of the two eyes)
    hor_diff_pred{i} = 0.5*(yle_pred{i} - yre_pred{i}); % 0.5*difference between horizontal eye positions
    % actual eye position
    ver_mean{i} = nanmean([zle{i} , zre{i}],2); % mean vertical eye position (of the two eyes)
    hor_mean{i} = nanmean([yle{i} , yre{i}],2); % mean horizontal eye position
    ver_diff{i} = 0.5*(zle{i} - zre{i}); % 0.5*difference between vertical eye positions (of the two eyes)
    hor_diff{i} = 0.5*(yle{i} - yre{i}); % 0.5*difference between horizontal eye positions
    if 0 
        figure;plot(hor_mean_pred{i},ver_mean_pred{i})
        hold on;plot(hor_mean{i},ver_mean{i}); xlim([-50 50]); ylim([-90 10]);
    end
        
    % subject trajectory in eye coordinates (fixed on starting position)
    [YE_sub{i},ZE_sub{i}] = subject2eye(xmt{i},ymt{i},zt,0);

    % subject's stopping position in subject's FOV
    [hor_sub{i},ver_sub{i}] = subject2eye(xs{i},ys{i},zt,0);

    % eye position in eye coordinates fixed on starting position
    YE_eye{i} = eye_th_arena{i};
    ZE_eye{i} = atan2d(zt,eye_r_arena{i});
    if 0
    hold off; plot(hor_mean_pred{i}(1),ver_mean_pred{i}(1),'rx'); hold on;
    plot(YE_eye{i},ZE_eye{i},'b'); plot(YE_eye{i}(1),ZE_eye{i}(1),'mo','markerfacecolor','m'); 
    plot(YE_eye{i}(wdrawindx(i)),ZE_eye{i}(wdrawindx(i)),'ro','markerfacecolor','r');
    plot(YE_sub{i},ZE_sub{i},'g'); plot(YE_sub{i}(end),ZE_sub{i}(end),'ko','markerfacecolor','k');
    xlim([-30 30]); ylim([-90 0]);
    end
    
    if 0
    % eye r and theta relative to starting eye position (eye coordinates)
    eye_r_screen{i} = sqrt((ver_mean{i} - ver_mean{i}(1)).^2 + (hor_mean{i} - hor_mean{i}(1)).^2);
    eye_th_screen{i} = atan2d(hor_mean{i}(1) - hor_mean{i},ver_mean{i}(1) - ver_mean{i});

    % target r and theta relative to starting target position (eye coordinates)
    tar_r_screen{i} = sqrt((ver_mean_pred{i} - ver_mean_pred{i}(1)).^2 + (hor_mean_pred{i} - hor_mean_pred{i}(1)).^2);
    tar_th_screen{i} = atan2d(hor_mean_pred{i}(1) - hor_mean_pred{i},ver_mean_pred{i}(1) - ver_mean_pred{i});
    else
    % r and theta relative to transformed (0,0) position in eye coordinates
    eye_r_screen{i} = sqrt((ver_mean{i} - z0).^2 + (hor_mean{i} - y0).^2);
    eye_th_screen{i} = atan2d(hor_mean{i} - y0,ver_mean{i} - z0);
    tar_r_screen{i} = sqrt((ver_mean_pred{i} - z0).^2 + (hor_mean_pred{i} - y0).^2);
    tar_th_screen{i} = atan2d(hor_mean_pred{i} - y0, ver_mean_pred{i} - z0);
    end
    
    % reverse transform
    if 0
        [x_tar,y_tar] = eye2subject(yle_pred{i},zle_pred{i},yre_pred{i},zre_pred{i},delta,zt);
        hold on; plot(x_tar(1),y_tar(1),'kx'); plot(x_tar,y_tar,'k'); hold off;
    end
    
    % sanity check
    if 0
    subplot(1,2,1);hold off;    plot(x_monk{i},y_monk{i})
    hold on;                    plot(trials(i).prs.fireflyposx,trials(i).prs.fireflyposy,'rx');

    subplot(1,2,2);plot(hor_mean_pred{i},ver_mean_pred{i})
    hold on;plot(hor_mean_pred{i}(1),ver_mean_pred{i}(1),'ro');
    ind = find(~isnan(hor_mean{i}),1);
    plot(hor_mean{i}(ind),ver_mean{i}(ind),'k+');
    plot(hor_mean{i}(1:50),ver_mean{i}(1:50),'ko');
    hold off; xlim([-30 30]);
    end
    
    % saccade direction
    sacxy{i} = []; sacxy_pred{i} = []; sacdir{i} = []; sacdir_pred{i} = []; sac_time{i} = []; sac_trl{i} = []; errxy{i} = []; errdir{i} = [];
    if ~isempty(timeindx)
        for j=1:length(t_sac2)
            sacstartindx = find(ts{i}>(t_sac2(j) - saccade_duration/2), 1) - timeindx(1) - 2;
            sacendindx = find(ts{i}>(t_sac2(j) + saccade_duration/2), 1) - timeindx(1) + 2;
            if sacstartindx>0 && sacendindx<=(numel(timeindx)-50) % only consider saccades 300ms (= 50 samples) before stopping
                sacxy{i} = [sacxy{i} [ver_mean{i}(sacendindx) - ver_mean{i}(sacstartindx); 
                    hor_mean{i}(sacendindx) - hor_mean{i}(sacstartindx)]];
                sacdir{i} = [sacdir{i} atan2d(sacxy{i}(1,end),sacxy{i}(2,end))];
                sacxy_pred{i} = [sacxy_pred{i} [ver_mean_pred{i}(sacstartindx) - ver_mean{i}(sacstartindx); 
                    hor_mean_pred{i}(sacstartindx) - hor_mean{i}(sacstartindx)]];
                sacdir_pred{i} = [sacdir_pred{i} atan2d(sacxy_pred{i}(1,end),sacxy_pred{i}(2,end))];
                sac_time{i} = [sac_time{i} (sacstartindx+sacendindx)/2*dt]; % time since target fixation
                sac_trl{i} = [sac_trl{i} i];
            end
        end
    end
end

%% save saccadic eye movements
eye_movement.saccade.true.val = cell2mat(sacxy);
eye_movement.saccade.true.dir = cell2mat(sacdir);
eye_movement.saccade.pred.val = cell2mat(sacxy_pred);
eye_movement.saccade.pred.dir = cell2mat(sacdir_pred);
eye_movement.saccade.time = cell2mat(sac_time);
eye_movement.saccade.trl = cell2mat(sac_trl);

%% correlation between behv error and eye-movement prediction error
Nt = max(cellfun(@(x) length(x),tar_r_sub)); % max number of timepoints
eye_mean_err = nan(ntrls, 1);
Wn = 2;

for i=1:ntrls
    nt = length(ver_mean{i}); %     if t_cut > nt;  t_cut = nt; end
    
    distind = find(d_sub{i} >= DistPerc*max(d_sub{i}),1);
    t_cut = nt - distind;     %        t_cut = 6; t_cut = round(t_cut/dt);

    %     eye_mean_err(i) = nanmean(sqrt((ver_mean{i}(1:nt-t_cut) -
    %     ver_mean_pred{i}(1:nt-t_cut)).^2 + (hor_mean{i}(1:nt-t_cut) - hor_mean_pred{i}(1:nt-t_cut)).^2));
    eye_mean_err(i) = sqrt(nanmean((ver_mean_pred{i}(1:nt-t_cut) - ver_mean{i}(1:nt-t_cut)).^2 + (hor_mean_pred{i}(1:nt-t_cut) - hor_mean{i}(1:nt-t_cut)).^2));
    
    TT_err_temp = sqrt((ver_mean{i}(1:nt-t_cut) - ver_mean_pred{i}(1:nt-t_cut)).^2 + (hor_mean{i}(1:nt-t_cut) - hor_mean_pred{i}(1:nt-t_cut)).^2);
    TT_err(i,:) = [movmean(TT_err_temp,Wn,'omitnan') ; nan(Nt - length(TT_err_temp),1)];
    
    TT_err_ver(i,:) = [ver_mean_pred{i}(1:nt-t_cut) - ver_mean{i}(1:nt-t_cut) ; nan(Nt - length(TT_err_temp),1)];
    sign_x = sign(hor_mean_pred{i}(1));
    TT_err_hor(i,:) = [(hor_mean_pred{i}(1:nt-t_cut) - hor_mean{i}(1:nt-t_cut)) ; nan(Nt - length(TT_err_temp),1)];
    
    TT_err_rate{i} = [0  diff(movmean(TT_err(i,:),Wn,'omitnan'))/dt];
                          
    TT_err_mean_rate(i) = nanmean(TT_err_rate{i});
    
    % trial errors (in screen coordinates)
    if prs.alldeg
        tempbehv_ver = ver_mean_pred{i}(1:nt-t_cut) - ver_sub{i}(1:nt-t_cut);
        tempbehv_hor = hor_mean_pred{i}(1:nt-t_cut) - hor_sub{i}(1:nt-t_cut);
        tempbehv_all = sqrt(tempbehv_ver.^2 + tempbehv_hor.^2);
        
        trlerrors.ver(i,:) = [tempbehv_ver ; nan(Nt - length(tempbehv_ver),1)];
        trlerrors.hor(i,:) = [tempbehv_hor ; nan(Nt - length(tempbehv_hor),1)];
        trlerrors.all(i,:) = [tempbehv_all ; nan(Nt - length(tempbehv_all),1)];
    else
        trlerrors.all(i,:) = behverrors(i)*ones(1,Nt);
    end

%
%     plot(movmean(TT_err{i},50,'omitnan')); 
%     plot(diff(movmean(TT_err{i},50,'omitnan'))/dt);
%     plot(nanmean(diff(movmean(TT_err{i},50,'omitnan'))/dt),'rx');
%     plot(eye_mean_err(i),'ro')
end
TT_err_mean = nanmean(TT_err);
TT_err_std = nanstd(TT_err);
ts = [1:Nt]*dt;

[eye_movement.eyepos.behvcorr.r,eye_movement.eyepos.behvcorr.p] = nancorr(behverrors(:), eye_mean_err(:));
eye_movement.eyepos.behvcorr.eye_err = eye_mean_err;
% eye_movement.eyepos.behvcorr.eye_err_rate = TT_err_mean_rate;
eye_movement.eyepos.behvcorr.eye_err_ts.ver = TT_err_ver';
eye_movement.eyepos.behvcorr.eye_err_ts.hor = TT_err_hor';
eye_movement.eyepos.behvcorr.eye_err_ts.all = TT_err';
eye_movement.eyepos.behvcorr.eye_err_ts.mu = TT_err_mean(:);
eye_movement.eyepos.behvcorr.eye_err_ts.std = TT_err_std(:);
eye_movement.eyepos.behvcorr.ts = ts(:);
eye_movement.eyepos.behvcorr.behv_err.all = trlerrors.all;
if prs.alldeg
eye_movement.eyepos.behvcorr.behv_err.ver = trlerrors.ver;
eye_movement.eyepos.behvcorr.behv_err.hor = trlerrors.hor;    
end
%% Withdrawal of eyes from target
% % choose earliest of 1)max eye distance, 2) max eye error rate
% for i=1:ntrls
%     distind = find(d_sub{i} > .5*max(d_sub{i}),1);
%     ind = find(~isnan(ver_mean{i}),1);
%     d_ver_eye1 = ver_mean{i}(ind) - ver_mean{i};
%     d_ver_eye = movmean(d_ver_eye1,50,'omitnan');
%     d_ver_rate = [0; diff(movmean(d_ver_eye1,60,'omitnan'))];
%     
%     [~,ind1] = max(d_ver_eye(100:end)); ind1 = 100 + ind1 - 1;% maximum vertical distance
%     
%     % eye movement speed
%     ver_speed = [0; diff(movmean(ver_mean{i},40,'omitnan'))]./dt; ver_speed = movmean(ver_speed,10,'omitnan');
%     hor_speed = [0; diff(movmean(hor_mean{i},40,'omitnan'))]./dt; hor_speed = movmean(hor_speed,10,'omitnan');
%     eye_speed = sqrt(ver_speed.^2 + hor_speed.^2);
%     [~,ind2] = findpeaks(eye_speed(100:end),'minpeakprominence',10,'minpeakwidth',15,'minpeakheight',15);
%     if ~isempty(ind2); ind2 = ind2(1) + 100 - 30; else; ind2 = nan; end
%     
%     ind3 = find(eye_speed(120:end) < 1.5); ind3 = 120 + ind3;
%     temp = diff(ind3); temp(temp > 1) = 0; 
%     if length(temp)>20
%         [~,loc,~,~] = findpeaks(temp,'minpeakwidth',20); 
%         if ~isempty(loc); ind3 = ind3(loc(1)); else; ind3 = nan; end
%     else; ind3 = nan; 
%     end
%     
%     indx = min([ind1 ind2 ind3]);
%     indx = ind1;
%     % rate of tracking error increase
% %     H = 3;
% %     W = 20;
% %     [pk,loc,w,pr] = findpeaks(TT_err_rate{i},'minpeakheight',H,'minpeakwidth',W);
% %     [~,ind2] = min((loc - ind1).^2);
% %     nback = find(flip(TT_err_rate{i}(1:loc(ind2))) <= H, 1);
% %     ind2 = loc(ind2) - nback;
% %     if isempty(ind2); ind2 = nan; end
% %     
% %     indx  = nanmedian([ind1 ind2]);
% %     indx = round(indx);
%     wdrawindx(i) = indx;
% %     if indx < 100
% %         keyboard
% %     end
%     if 0 
%         figure;plot(hor_mean_pred{i},ver_mean_pred{i})
%         hold on;plot(hor_mean{i},ver_mean{i},'k'); plot(hor_mean{i}(1),ver_mean{i}(1),'ko'); 
%         plot(hor_mean{i}(indx),ver_mean{i}(indx),'ro')
%         figure;hold on;plot(ver_speed);plot(hor_speed);plot(eye_speed,'k'); hline(0,'k--'); vline(indx);
%         figure;plot(TT_err_rate{i}); title('TT error rate');vline(indx)
%         figure;plot(d_ver_eye); title(['trial ' num2str(i) ' - vertical distance']);vline(indx)
%         hold on;plot(d_ver_rate*10); 
%     end
% end
%% save true and predicted eye positions
for i=1:ntrls
    eye_movement.arena.flypos.r.val{i} = downsample(tar_r_sub{i},factor_downsample);
    eye_movement.arena.flypos.theta.val{i} = downsample(tar_th_sub{i},factor_downsample);
    eye_movement.arena.flypos.x.val{i} = downsample(xt{i},factor_downsample);
    eye_movement.arena.flypos.y.val{i} = downsample(yt{i},factor_downsample);
    eye_movement.arena.eyepos.r.val{i} = downsample(eye_r_arena{i},factor_downsample);
    eye_movement.arena.eyepos.theta.val{i} = downsample(eye_th_arena{i},factor_downsample);
    eye_movement.arena.eyepos.x.val{i} = downsample(x_mean_obs{i},factor_downsample);
    eye_movement.arena.eyepos.y.val{i} = downsample(y_mean_obs{i},factor_downsample);
    eye_movement.arena.subpos.r.val{i} = downsample(sub_r{i},factor_downsample);
    eye_movement.arena.subpos.theta.val{i} = downsample(sub_th{i},factor_downsample);
    eye_movement.arena.subpos.x.val{i} = downsample(x_monk{i},factor_downsample);
    eye_movement.arena.subpos.y.val{i} = downsample(y_monk{i},factor_downsample);
    eye_movement.arena.subpos.w.val{i} = downsample(w{i},factor_downsample);
    eye_movement.arena.sub2eye.r.val{i} = downsample(r_sub2eye{i},factor_downsample);
    eye_movement.arena.sub2eye.theta.val{i} = downsample(th_sub2eye{i},factor_downsample);

%     eye_movement.screen.flypos.ver.val(i) = ver_mean_pred{i}(1);
%     eye_movement.screen.flypos.hor.val(i) = hor_mean_pred{i}(1);
%     eye_movement.screen.subpos.ver.val{i} = downsample(ZE_sub{i},factor_downsample);
%     eye_movement.screen.subpos.hor.val{i} = downsample(YE_sub{i},factor_downsample);
%     eye_movement.screen.eyepos.ver.val{i} = downsample(ZE_eye{i},factor_downsample);
%     eye_movement.screen.eyepos.hor.val{i} = downsample(YE_eye{i},factor_downsample);
    
    eye_movement.subdis.d{i} = [downsample(d_sub{i},factor_downsample) ; nan(Nt - length(downsample(d_sub{i},factor_downsample)),1)];
    eye_movement.eyepos.pred.ver_mean.val{i} = downsample(ver_mean_pred{i},factor_downsample);
    eye_movement.eyepos.pred.hor_mean.val{i} = downsample(hor_mean_pred{i},factor_downsample);
    eye_movement.eyepos.pred.ver_diff.val{i} = downsample(ver_diff_pred{i},factor_downsample);
    eye_movement.eyepos.pred.hor_diff.val{i} = downsample(hor_diff_pred{i},factor_downsample);
    eye_movement.eyepos.pred.r.val{i} = downsample(tar_r_screen{i},factor_downsample);
    eye_movement.eyepos.pred.theta.val{i} = downsample(tar_th_screen{i},factor_downsample);
    
    eye_movement.eyepos.true.ver_mean.val{i} = downsample(ver_mean{i},factor_downsample);
    eye_movement.eyepos.true.hor_mean.val{i} = downsample(hor_mean{i},factor_downsample);
    eye_movement.eyepos.true.ver_diff.val{i} = downsample(ver_diff{i},factor_downsample);
    eye_movement.eyepos.true.hor_diff.val{i} = downsample(hor_diff{i},factor_downsample);
    eye_movement.eyepos.true.r.val{i} = downsample(eye_r_screen{i},factor_downsample);
    eye_movement.eyepos.true.theta.val{i} = downsample(eye_th_screen{i},factor_downsample);

    eye_movement.eyepos.blv.ver_mean.val{i} = downsample(ver_sub{i},factor_downsample);
    eye_movement.eyepos.blv.hor_mean.val{i} = downsample(hor_sub{i},factor_downsample);
    
    eye_movement.eyepos.target_off.indx(i) = target_off(i);
%     eye_movement.eyepos.withdraw.indx(i) = wdrawindx(i);
    eye_movement.tau.val(i) = tau(i);
    eye_movement.endofsteering.val(i) = endofsteering(i);
end

Nt = max(cellfun(@(x) length(x),tar_r_sub)); % max number of timepoints
%% compare predicted & eye positions from all trials aligned to target fixation
% gather predicted eye position
ver_mean_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean_pred,'UniformOutput',false));
hor_mean_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean_pred,'UniformOutput',false));
ver_diff_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff_pred,'UniformOutput',false));
hor_diff_pred1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff_pred,'UniformOutput',false));
% compute variance of predicted eye position
if ~isempty(spatialerr)
ver_var_pred1 = ComputeVarianceEyeVer(xt,yt,zt,var_xmt,var_ymt,'normal');
hor_var_pred1 = ComputeVarianceEyeHor(xt,yt,zt,var_xmt,var_ymt,'normal');
end
% gather true eye position
ver_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_mean,'UniformOutput',false));
hor_mean1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_mean,'UniformOutput',false));
ver_diff1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],ver_diff,'UniformOutput',false));
hor_diff1 = cell2mat(cellfun(@(x) [x(:) ; nan(Nt - length(x),1)],hor_diff,'UniformOutput',false));

% cross-correlation between observed and predicted eye positions
zeropad_length = 200; % pad these many bins around each trial
maxlag = 400;
pred = [zeros(zeropad_length,ntrls); ver_mean_pred1; zeros(zeropad_length,ntrls)]; 
pred_shuffled = pred(:,randperm(ntrls));
pred = pred(:); pred(isnan(pred)) = 0; 
pred_shuffled = pred_shuffled(:); pred_shuffled(isnan(pred_shuffled)) = 0;
obs = [zeros(zeropad_length,ntrls); ver_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;
[eye_movement.eyepos.pred_vs_true.ver_mean.xcorr.c, eye_movement.eyepos.pred_vs_true.ver_mean.xcorr.lag] = xcorr(pred,obs,maxlag,'coeff');
eye_movement.eyepos.pred_vs_true.ver_mean.xcorr.c_shuffled = xcorr(pred_shuffled,obs,maxlag,'coeff');

pred = [zeros(zeropad_length,ntrls); hor_mean_pred1; zeros(zeropad_length,ntrls)]; 
pred_shuffled = pred(:,randperm(ntrls));
pred = pred(:); pred(isnan(pred)) = 0; 
pred_shuffled = pred_shuffled(:); pred_shuffled(isnan(pred_shuffled)) = 0;
obs = [zeros(zeropad_length,ntrls); hor_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;
[eye_movement.eyepos.pred_vs_true.hor_mean.xcorr.c, eye_movement.eyepos.pred_vs_true.hor_mean.xcorr.lag] = xcorr(pred,obs,maxlag,'coeff');
eye_movement.eyepos.pred_vs_true.hor_mean.xcorr.c_shuffled = xcorr(pred_shuffled,obs,maxlag,'coeff');

% timecourse of component-wise corr between predicted & true eye position
% [eye_movement.eyepos.pred_vs_true.ver_mean.rho.startaligned,eye_movement.eyepos.pred_vs_true.ver_mean.pval.startaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.hor_mean.rho.startaligned,eye_movement.eyepos.pred_vs_true.hor_mean.pval.startaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.ver_diff.rho.startaligned,eye_movement.eyepos.pred_vs_true.ver_diff.pval.startaligned] = arrayfun(@(i) corr(ver_diff1(i,:)',ver_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.hor_diff.rho.startaligned,eye_movement.eyepos.pred_vs_true.hor_diff.pval.startaligned] = arrayfun(@(i) corr(hor_diff1(i,:)',hor_diff_pred1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% 
% timecourse of component-wise regression between predicted & true eye position
% beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta10.startaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta00.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta10.startaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta00.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_diff1(i,:)',[ver_diff_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta10.startaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta00.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_diff1(i,:)',[hor_diff_pred1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta10.startaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta00.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta1.startaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta0.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta1.startaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta0.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_diff1(i,:)',[ver_diff_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta1.startaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta0.startaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_diff1(i,:)',[hor_diff_pred1(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta1.startaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta0.startaligned] = deal(beta(1,:),beta(2,:));

% timecourse of cosine similarity between predicted & true eye position
% pred = permute(cat(3,ver_mean_pred1 , hor_mean_pred1),[3 1 2]);
% true = permute(cat(3,ver_mean1 , hor_mean1),[3 1 2]);
% cos_similarity = nan(Nboots,Nt); cos_similarity_shuffled = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1); 
%     randtrls2 = randsample(ntrls,ntrls,1);
%     cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
%     cos_similarity_shuffled(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls2)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cos_similarity.mu.startaligned = mean(cos_similarity);
% eye_movement.eyepos.pred_vs_true.cos_similarity.sem.startaligned = std(cos_similarity);
% eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.mu.startaligned = mean(cos_similarity_shuffled);
% eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.sem.startaligned = std(cos_similarity_shuffled);

% timecourse of cosine similarity between predicted & true eye position
% ngroups = 5;
% ntrls_per_group = (ntrls - mod(ntrls,ngroups))/ngroups;
% errorindx = errorindx(1:ntrls_per_group*ngroups);
% cos_similarity = nan(ngroups,Nt);
% for i=1:ngroups
%     trlgroup = errorindx(ntrls_per_group*(i-1) + 1:ntrls_per_group*i);
%     cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,trlgroup),true(:,i,trlgroup)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cos_similarity_vs_error.mu.startaligned = cos_similarity;

% timecourse of centered cosine similarity between predicted & true eye position -
% pred = permute(cat(3,ver_mean_pred1 , hor_mean_pred1),[3 1 2]); pred = (pred - repmat(nanmean(pred,3),[1 1 ntrls]))./repmat(nanstd(pred,[],3),[1 1 ntrls]);
% true = permute(cat(3,ver_mean1 , hor_mean1),[3 1 2]); true = (true - repmat(nanmean(true,3),[1 1 ntrls]))./repmat(nanstd(true,[],3),[1 1 ntrls]);
% cntr_cos_similarity = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1);
%     cntr_cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.mu.startaligned = mean(cntr_cos_similarity);
% eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.sem.startaligned = std(cntr_cos_similarity);

% timecourse of variance explained
pred = permute(cat(3,ver_mean_pred1 , hor_mean_pred1),[3 1 2]);
obs = permute(cat(3,ver_mean1 , hor_mean1),[3 1 2]);
if boots
    var_explained = nan(Nboots,Nt);
    var_explained_shuffled = nan(Nboots,Nt);
    for i=1:Nboots
        randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
        squared_err = sum(nanmean((obs(:,:,randtrls) - pred(:,:,randtrls)).^2,3));
        var_pred = sum(nanvar(pred(:,:,randtrls),[],3));
        var_explained(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - pred(:,:,randtrls2)).^2,3));
        var_pred_shuffled = sum(nanvar(pred(:,:,randtrls2),[],3));
        var_explained_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
    end
    eye_movement.eyepos.pred_vs_true.var_explained.mu.startaligned = nanmean(var_explained);
    eye_movement.eyepos.pred_vs_true.var_explained.sem.startaligned = nanstd(var_explained);
    eye_movement.eyepos.pred_vs_true.var_explained_shuffled.mu.startaligned = nanmean(var_explained_shuffled);
    eye_movement.eyepos.pred_vs_true.var_explained_shuffled.sem.startaligned = nanstd(var_explained_shuffled);
else
    if norm_TTI
    squared_err = sum((obs - pred).^2);
    var_pred = sum(nanvar(pred,[],3));
    var_explained = squeeze(1 - (squared_err./var_pred)); % try taking sqrt or cosine of var_explained
    init_val = var_explained(1,:);
    sem_var_explained = nanstd(var_explained./init_val,[],2)./sqrt(size(pred,3));
    var_explained = nanmean(var_explained./init_val,2);
    else
    if 0
    squared_err = sum((obs - pred).^2);
    var_pred = sum(nanvar(pred,[],3));
    var_explained = squeeze(1 - (squared_err./var_pred)); % try taking sqrt or cosine of var_explained
    sem_var_explained = nanstd(var_explained,[],2)./sqrt(size(pred,3));
    var_explained = nanmean(var_explained,2);
    
    % Kaushik's way
    else
    squared_err = sum(nanmean((obs - pred).^2,3));
    var_pred = sum(nanvar(pred,[],3));
    var_explained = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    end
    end
    eye_movement.eyepos.pred_vs_true.var_explained.mu.startaligned = var_explained;    
%     eye_movement.eyepos.pred_vs_true.var_explained.sem.startaligned = sem_var_explained;
end

% timecourse of upper bound on variance explained
if ~isempty(spatialerr)
expected_squared_err = (nanmean(ver_var_pred1,2) + nanmean(hor_var_pred1,2)); 
var_pred = sum(nanvar(pred,[],3)); 
var_true = sum(nanvar(obs,[],3));
eye_movement.eyepos.pred_vs_true.var_explained_upperbound.mu.startaligned = 1 - (expected_squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
eye_movement.eyepos.pred_vs_true.expected_squared_err.mu.startaligned = expected_squared_err(:);
eye_movement.eyepos.pred_vs_true.var_pred.mu.startaligned = var_pred(:);
eye_movement.eyepos.pred_vs_true.var_true.mu.startaligned = var_true(:);
end
% timecourse of variance explained for various accuracies
ntrls_per_group = (ntrls - mod(ntrls,ngroups))/ngroups;
errorindx = errorindx(1:ntrls_per_group*ngroups);
var_explained = nan(ngroups,Nt);
for i=1:ngroups
    trlgroup = errorindx(ntrls_per_group*(i-1) + 1:ntrls_per_group*i);
    squared_err = sum(nanmean((obs(:,:,trlgroup) - pred(:,:,trlgroup)).^2,3));
    var_pred = sum(nanvar(pred(:,:,trlgroup),[],3));
    var_explained(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    trl_dur(i,:) = t_stop(trlgroup);
    trl_err(i,:) = behverrors(trlgroup);
end
eye_movement.eyepos.pred_vs_true.var_explained_vs_error.mu.startaligned = var_explained;
eye_movement.eyepos.pred_vs_true.var_explained_vs_error.trldur = trl_dur;
eye_movement.eyepos.pred_vs_true.var_explained_vs_error.trlerr = trl_err;

% timecourse of similarity measure based on Mahalanobis distance
% dist_mahalanobis = nan(Nboots,Nt); dist_mahalanobis_shuffled = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1);
%     [~,indx1] = min(abs(squeeze(true(1,:,randtrls)) - permute(repmat(eye_fixation.eyepos.true.ver_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
%     [~,indx2] = min(abs(squeeze(true(2,:,randtrls)) - permute(repmat(eye_fixation.eyepos.true.hor_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
%     ver_mahalanobis = nanmean(((squeeze(true(1,:,randtrls)) - squeeze(pred(1,:,randtrls))).^2)./(eye_fixation.eyepos.true.ver_mean.sig(indx1).^2),2);
%     hor_mahalanobis = nanmean(((squeeze(true(2,:,randtrls)) - squeeze(pred(2,:,randtrls))).^2)./(eye_fixation.eyepos.true.hor_mean.sig(indx2).^2),2);
%     randtrls2 = randsample(ntrls,ntrls,1);
%     ver_mahalanobis_shuffled = nanmean(((squeeze(true(1,:,randtrls)) - squeeze(pred(1,:,randtrls2))).^2)./(eye_fixation.eyepos.true.ver_mean.sig(indx1).^2),2);
%     hor_mahalanobis_shuffled = nanmean(((squeeze(true(2,:,randtrls)) - squeeze(pred(2,:,randtrls2))).^2)./(eye_fixation.eyepos.true.hor_mean.sig(indx2).^2),2);
%     dist_mahalanobis(i,:) = sqrt(ver_mahalanobis + hor_mahalanobis);
%     dist_mahalanobis_shuffled(i,:) = sqrt(ver_mahalanobis_shuffled + hor_mahalanobis_shuffled);
% end
% dist_mahalanobis = nanmean(dist_mahalanobis); d0 = min(dist_mahalanobis(1:50)); % reused for stopaligned analysis
% eye_movement.eyepos.pred_vs_true.mahalanobis_distance.mu.startaligned = dist_mahalanobis;
% eye_movement.eyepos.pred_vs_true.mahalanobis_similarity.mu.startaligned = 2./(1 + dist_mahalanobis/d0);
% dist_mahalanobis_shuffled = nanmean(dist_mahalanobis_shuffled);
% eye_movement.eyepos.pred_vs_true.mahalanobis_distance_shuffled.mu.startaligned = dist_mahalanobis_shuffled;
% eye_movement.eyepos.pred_vs_true.mahalanobis_similarity_shuffled.mu.startaligned = 2./(1 + dist_mahalanobis_shuffled/d0);

% timecourse of upper bound on similarity measure based on Mahalanobis distance
% [~,indx1] = min(abs(squeeze(true(1,:,:)) - permute(repmat(eye_fixation.eyepos.true.ver_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
% ver_mahalanobis_expected = nanmean(ver_var_pred1./(eye_fixation.eyepos.true.hor_mean.sig(indx1).^2),2);
% [~,indx2] = min(abs(squeeze(true(2,:,:)) - permute(repmat(eye_fixation.eyepos.true.hor_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
% hor_mahalanobis_expected = nanmean(hor_var_pred1./(eye_fixation.eyepos.true.ver_mean.sig(indx2).^2),2);
% dist_mahalanobis_expected = sqrt(ver_mahalanobis_expected + hor_mahalanobis_expected);
% eye_movement.eyepos.pred_vs_true.mahalanobis_upperbound.mu.startaligned = 2./(1 + dist_mahalanobis_expected/d0);

%% compare predicted & observed eye positions from all trials aligned to END of movement
% predicted eye position

ver_mean_pred2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],ver_mean_pred,'UniformOutput',false));
hor_mean_pred2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],hor_mean_pred,'UniformOutput',false));
ver_diff_pred2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],ver_diff_pred,'UniformOutput',false));
hor_diff_pred2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],hor_diff_pred,'UniformOutput',false));
% compute variance of predicted eye position
if ~isempty(spatialerr)
ver_var_pred2 = ComputeVarianceEyeVer(xt,yt,zt,var_xmt,var_ymt,'reverse');
hor_var_pred2 = ComputeVarianceEyeHor(xt,yt,zt,var_xmt,var_ymt,'reverse');
end
% true eye position
ver_mean2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],ver_mean,'UniformOutput',false));
hor_mean2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],hor_mean,'UniformOutput',false));
ver_diff2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],ver_diff,'UniformOutput',false));
hor_diff2 = cell2mat(cellfun(@(x) [flipud(x(:)) ; nan(Nt - length(x),1)],hor_diff,'UniformOutput',false));

% Adjust for tau-induced variable length of trial (high tau, higher difficulty to stop)
for j = 1:length(tau)
    if TauAdjust
        cutoff = round(Ntaus*tau(j)/dt);
    elseif PeakVel_cutoff
        [~,peakind] = max(diff(d_sub{j})/dt);
        cutoff = length(d_sub{j}) - peakind;
    elseif Dist_cutoff
        distind = find(d_sub{j} >= DistPerc*max(d_sub{j}),1);
        cutoff = length(d_sub{j}) - distind;
    else
        cutoff = 4/dt;
    end
    if cutoff < 1
        cutoff = 1;
    end
    ver_mean_pred2(:,j) = [ver_mean_pred2(cutoff:end,j) ; nan(cutoff-1,1)];
    hor_mean_pred2(:,j) = [hor_mean_pred2(cutoff:end,j) ; nan(cutoff-1,1)];
    ver_diff_pred2(:,j) = [ver_diff_pred2(cutoff:end,j) ; nan(cutoff-1,1)];
    hor_diff_pred2(:,j) = [hor_diff_pred2(cutoff:end,j) ; nan(cutoff-1,1)];
    
    ver_mean2(:,j) = [ver_mean2(cutoff:end,j) ; nan(cutoff-1,1)];
    hor_mean2(:,j) = [hor_mean2(cutoff:end,j) ; nan(cutoff-1,1)];
    ver_diff2(:,j) = [ver_diff2(cutoff:end,j) ; nan(cutoff-1,1)];
    hor_diff2(:,j) = [hor_diff2(cutoff:end,j) ; nan(cutoff-1,1)];
end

% timecourse of component-wise corr between predicted & true eye position
% [eye_movement.eyepos.pred_vs_true.ver_mean.rho.stopaligned,eye_movement.eyepos.pred_vs_true.ver_mean.pval.stopaligned] = arrayfun(@(i) corr(ver_mean2(i,:)',ver_mean_pred2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.hor_mean.rho.stopaligned,eye_movement.eyepos.pred_vs_true.hor_mean.pval.stopaligned] = arrayfun(@(i) corr(hor_mean2(i,:)',hor_mean_pred2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.ver_diff.rho.stopaligned,eye_movement.eyepos.pred_vs_true.ver_diff.pval.stopaligned] = arrayfun(@(i) corr(ver_diff2(i,:)',ver_diff_pred2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% [eye_movement.eyepos.pred_vs_true.hor_diff.rho.stopaligned,eye_movement.eyepos.pred_vs_true.hor_diff.pval.stopaligned] = arrayfun(@(i) corr(hor_diff2(i,:)',hor_diff_pred2(i,:)','Type','Spearman','rows','complete'), 1:Nt);

% component-wise regression between predicted & true eye position
% beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_pred2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta10.stopaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta00.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_pred2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta10.stopaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta00.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_diff2(i,:)',[ver_diff_pred2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta10.stopaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta00.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_diff2(i,:)',[hor_diff_pred2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta10.stopaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta00.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_pred2(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_mean.beta1.stopaligned,eye_movement.eyepos.pred_vs_true.ver_mean.beta0.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_pred2(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_mean.beta1.stopaligned,eye_movement.eyepos.pred_vs_true.hor_mean.beta0.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(ver_diff2(i,:)',[ver_diff_pred2(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.ver_diff.beta1.stopaligned,eye_movement.eyepos.pred_vs_true.ver_diff.beta0.stopaligned] = deal(beta(1,:),beta(2,:));
% beta = cell2mat(arrayfun(@(i) regress(hor_diff2(i,:)',[hor_diff_pred2(i,:)' ones(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [eye_movement.eyepos.pred_vs_true.hor_diff.beta1.stopaligned,eye_movement.eyepos.pred_vs_true.hor_diff.beta0.stopaligned] = deal(beta(1,:),beta(2,:));

% timecourse of cosine similarity between predicted & true eye position
% pred = permute(cat(3,ver_mean_pred2 , hor_mean_pred2),[3 1 2]);
% true = permute(cat(3,ver_mean2 , hor_mean2),[3 1 2]);
% cos_similarity = nan(Nboots,Nt); cos_similarity_shuffled = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
%     cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
%     cos_similarity_shuffled(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls2)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cos_similarity.mu.stopaligned = mean(cos_similarity);
% eye_movement.eyepos.pred_vs_true.cos_similarity.sem.stopaligned = std(cos_similarity);
% eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.mu.stopaligned = mean(cos_similarity_shuffled);
% eye_movement.eyepos.pred_vs_true.cos_similarity_shuffled.sem.stopaligned = std(cos_similarity_shuffled);

% timecourse of cosine similarity between predicted & true eye position -
% ngroups = 5;
% ntrls_per_group = (ntrls - mod(ntrls,ngroups))/ngroups;
% errorindx = errorindx(1:ntrls_per_group*ngroups);
% cos_similarity = nan(ngroups,Nt);
% for i=1:ngroups
%     trlgroup = errorindx(ntrls_per_group*(i-1) + 1:ntrls_per_group*i);
%     cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,trlgroup),true(:,i,trlgroup)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cos_similarity_vs_error.mu.stopaligned = cos_similarity;

% timecourse of centered cosine similarity between predicted & true eye position -
% pred = permute(cat(3,ver_mean_pred2 , hor_mean_pred2),[3 1 2]); pred = (pred - repmat(nanmean(pred,3),[1 1 ntrls]))./repmat(nanstd(pred,[],3),[1 1 ntrls]);
% true = permute(cat(3,ver_mean2 , hor_mean2),[3 1 2]); true = (true - repmat(nanmean(true,3),[1 1 ntrls]))./repmat(nanstd(true,[],3),[1 1 ntrls]);
% cntr_cos_similarity = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1);
%     cntr_cos_similarity(i,:) = arrayfun(@(i) CosSimilarity(pred(:,i,randtrls),true(:,i,randtrls)), 1:Nt);
% end
% eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.mu.stopaligned = mean(cntr_cos_similarity);
% eye_movement.eyepos.pred_vs_true.cntr_cos_similarity.sem.stopaligned = std(cntr_cos_similarity);

% timecourse of variance explained
pred = permute(cat(3,ver_mean_pred2 , hor_mean_pred2),[3 1 2]);
obs = permute(cat(3,ver_mean2 , hor_mean2),[3 1 2]);
if boots
    var_explained = nan(Nboots,Nt); var_explained_shuffled = nan(Nboots,Nt);
    for i=1:Nboots
        randtrls = randsample(ntrls,ntrls,1); randtrls2 = randsample(ntrls,ntrls,1);
        squared_err = sum(nanmean((obs(:,:,randtrls) - pred(:,:,randtrls)).^2,3)); var_pred = sum(nanvar(pred(:,:,randtrls),[],3));
        var_explained(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
        squared_err_shuffled = sum(nanmean((obs(:,:,randtrls) - pred(:,:,randtrls2)).^2,3)); var_pred_shuffled = sum(nanvar(pred(:,:,randtrls2),[],3));
        var_explained_shuffled(i,:) =  1 - (squared_err_shuffled(:)./var_pred_shuffled(:));
    end
    eye_movement.eyepos.pred_vs_true.var_explained.mu.stopaligned = nanmean(var_explained);
    eye_movement.eyepos.pred_vs_true.var_explained.sem.stopaligned = nanstd(var_explained);
    eye_movement.eyepos.pred_vs_true.var_explained_shuffled.mu.stopaligned = nanmean(var_explained_shuffled);
    eye_movement.eyepos.pred_vs_true.var_explained_shuffled.sem.stopaligned = nanstd(var_explained_shuffled);
else
    if 0
    squared_err = sum((obs - pred).^2);
    var_pred = sum(nanvar(pred,[],3));
    var_explained = squeeze(1 - (squared_err./var_pred)); % try taking sqrt or cosine of var_explained
    sem_var_explained = nanstd(var_explained,[],2)./sqrt(size(pred,3));
    var_explained = nanmean(var_explained,2);
    
    % Kaushik's way
    else
    squared_err = sum(nanmean((obs - pred).^2,3));
    var_pred = sum(nanvar(pred,[],3));
    var_explained = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    end    
    eye_movement.eyepos.pred_vs_true.var_explained.mu.stopaligned = var_explained;
%     eye_movement.eyepos.pred_vs_true.var_explained.sem.stopaligned = sem_var_explained;
end

% timecourse of upper bound on variance explained
if ~isempty(spatialerr)
    expected_squared_err = (nanmean(ver_var_pred2,2) + nanmean(hor_var_pred2,2)); var_pred = sum(nanvar(pred,[],3)); var_true = sum(nanvar(obs,[],3));
    eye_movement.eyepos.pred_vs_true.var_explained_upperbound.mu.stopaligned = 1 - (expected_squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    eye_movement.eyepos.pred_vs_true.expected_squared_err.mu.stopaligned = expected_squared_err(:);
    eye_movement.eyepos.pred_vs_true.var_pred.mu.stopaligned = var_pred(:);
    eye_movement.eyepos.pred_vs_true.var_true.mu.stopaligned = var_true(:);
end
% timecourse of variance explained for various accuracies
ntrls_per_group = (ntrls - mod(ntrls,ngroups))/ngroups;
errorindx = errorindx(1:ntrls_per_group*ngroups);
var_explained = nan(ngroups,Nt);
for i=1:ngroups
    trlgroup = errorindx(ntrls_per_group*(i-1) + 1:ntrls_per_group*i);
    squared_err = sum(nanmean((obs(:,:,trlgroup) - pred(:,:,trlgroup)).^2,3)); var_pred = sum(nanvar(pred(:,:,trlgroup),[],3));
    var_explained(i,:) = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
end
eye_movement.eyepos.pred_vs_true.var_explained_vs_error.mu.stopaligned = var_explained;

% timecourse of similarity measure based on Mahalanobis distance
% dist_mahalanobis = nan(Nboots,Nt); dist_mahalanobis_shuffled = nan(Nboots,Nt);
% for i=1:Nboots
%     randtrls = randsample(ntrls,ntrls,1);
%     [~,indx1] = min(abs(squeeze(true(1,:,randtrls)) - permute(repmat(eye_fixation.eyepos.true.ver_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
%     [~,indx2] = min(abs(squeeze(true(2,:,randtrls)) - permute(repmat(eye_fixation.eyepos.true.hor_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
%     ver_mahalanobis = nanmean(((squeeze(true(1,:,randtrls)) - squeeze(pred(1,:,randtrls))).^2)./(eye_fixation.eyepos.true.ver_mean.sig(indx1).^2),2);
%     hor_mahalanobis = nanmean(((squeeze(true(2,:,randtrls)) - squeeze(pred(2,:,randtrls))).^2)./(eye_fixation.eyepos.true.hor_mean.sig(indx2).^2),2);
%     randtrls2 = randsample(ntrls,ntrls,1);
%     ver_mahalanobis_shuffled = nanmean(((squeeze(true(1,:,randtrls)) - squeeze(pred(1,:,randtrls2))).^2)./(eye_fixation.eyepos.true.ver_mean.sig(indx1).^2),2);
%     hor_mahalanobis_shuffled = nanmean(((squeeze(true(2,:,randtrls)) - squeeze(pred(2,:,randtrls2))).^2)./(eye_fixation.eyepos.true.hor_mean.sig(indx2).^2),2);
%     dist_mahalanobis(i,:) = sqrt(ver_mahalanobis + hor_mahalanobis);
%     dist_mahalanobis_shuffled(i,:) = sqrt(ver_mahalanobis_shuffled + hor_mahalanobis_shuffled);
% end
% dist_mahalanobis = nanmean(dist_mahalanobis);
% eye_movement.eyepos.pred_vs_true.mahalanobis_distance.mu.stopaligned = dist_mahalanobis;
% eye_movement.eyepos.pred_vs_true.mahalanobis_similarity.mu.stopaligned = 2./(1 + dist_mahalanobis/d0);
% dist_mahalanobis_shuffled = nanmean(dist_mahalanobis_shuffled);
% eye_movement.eyepos.pred_vs_true.mahalanobis_distance_shuffled.mu.stopaligned = dist_mahalanobis_shuffled;
% eye_movement.eyepos.pred_vs_true.mahalanobis_similarity_shuffled.mu.stopaligned = 2./(1 + dist_mahalanobis_shuffled/d0);

% timecourse of upper bound on similarity measure based on Mahalanobis distance
% [~,indx1] = min(abs(squeeze(true(1,:,:)) - permute(repmat(eye_fixation.eyepos.true.ver_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
% ver_mahalanobis_expected = nanmean(ver_var_pred2./(eye_fixation.eyepos.true.hor_mean.sig(indx1).^2),2);
% [~,indx2] = min(abs(squeeze(true(2,:,:)) - permute(repmat(eye_fixation.eyepos.true.hor_mean.mu(:),[1 Nt ntrls]),[2 3 1])),[],3);
% hor_mahalanobis_expected = nanmean(hor_var_pred2./(eye_fixation.eyepos.true.ver_mean.sig(indx2).^2),2);
% dist_mahalanobis_expected = sqrt(ver_mahalanobis_expected + hor_mahalanobis_expected);
% eye_movement.eyepos.pred_vs_true.mahalanobis_upperbound.mu.stopaligned = 2./(1 + dist_mahalanobis_expected/d0);
% 
