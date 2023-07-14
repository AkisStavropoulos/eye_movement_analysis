function tracking_dist2end = EyeDataOverDistanceEndAligned(tracking,dd,prs,align_flag)
%% Calculate same measures of eye movements as a function of distance traveled (end aligned)
boots = prs.boots;

% dd = 5; % spatial resolution of re-sampling

ntrls = numel(tracking.misc.behverrors.all.val);
Dmax = 1500; % max(tracking.misc.subdist,[],'all'); % maximum traveled distance across trials
ds = -flip(0:dd:Dmax);
Nt = numel(ds);

sacdist = [];
for j = 1:ntrls
    
    dist = tracking.misc.subdist(j,:);
    travelD = max(dist);
    dist2end = dist - travelD;
    
    ds_trl = ceil(dist2end):dd:0;
    
    ind = arrayfun(@(x) find(dist2end <= x,1,'last'), ds_trl);
    
    ver_mean(j,:) = [nan(1,Nt - numel(ind))  tracking.eyepos.screen.ver_mean(j,ind)];
    hor_mean(j,:) = [nan(1,Nt - numel(ind))  tracking.eyepos.screen.hor_mean(j,ind)];
    ver_mean_targ(j,:) = [nan(1,Nt - numel(ind))  tracking.tarpos.screen.ver_mean(j,ind)];
    hor_mean_targ(j,:) = [nan(1,Nt - numel(ind))  tracking.tarpos.screen.hor_mean(j,ind)];
    ver_mean_stop(j,:) = [nan(1,Nt - numel(ind))  tracking.stopos.screen.ver_mean(j,ind)];
    hor_mean_stop(j,:) = [nan(1,Nt - numel(ind))  tracking.stopos.screen.hor_mean(j,ind)];
    
    % align eye to either target or stopping position at target offset
    indx1 = find(~isnan(ver_mean(j,:)),1);
    if strcmp(align_flag,'align2targ')
    ver_mean(j,:) = ver_mean(j,:) + (ver_mean_targ(j,indx1) - ver_mean(j,indx1));
%     hor_mean(j,:) = hor_mean(j,:) + (hor_mean_targ(j,indx1) - hor_mean(j,indx1));    
    elseif strcmp(align_flag,'align2stop')
    ver_mean(j,:) = ver_mean(j,:) + (ver_mean_stop(j,indx1) - ver_mean(j,indx1));
%     hor_mean(j,:) = hor_mean(j,:) + (hor_mean_stop(j,indx1) - hor_mean(j,indx1));  
    elseif strcmp(align_flag,'alignall')
        ver_mean(j,:) = ver_mean(j,:) - ver_mean(j,1);
        hor_mean(j,:) = hor_mean(j,:) - hor_mean(j,1);
        ver_mean_targ(j,:) = ver_mean_targ(j,:) - ver_mean_targ(j,1);
        hor_mean_targ(j,:) = hor_mean_targ(j,:) - hor_mean_targ(j,1);
        ver_mean_stop(j,:) = ver_mean_stop(j,:) - ver_mean_stop(j,1);
        hor_mean_stop(j,:) = hor_mean_stop(j,:) - hor_mean_stop(j,1);
    end

    % get distance percentage of saccade occurence
    sactrl = tracking.saccade.trl == j;
    if sum(sactrl) > 0
        sactime = tracking.saccade.time(sactrl);
        sacindx = arrayfun(@(x) find(tracking.misc.ts <= x,1,'last'), sactime);
        sacdist = [sacdist dist2end(sacindx)];
    end

end

%% Get errors
% Target-Tracking Error
TTE_ver = ver_mean_targ - ver_mean ;
TTE_hor = hor_mean_targ - hor_mean ;
TTE_all = sqrt(TTE_ver.^2 + TTE_hor.^2);

% Stop-Position-Tracking Error
SPTE_ver = ver_mean_stop - ver_mean ;
SPTE_hor = hor_mean_stop - hor_mean ;
SPTE_all = sqrt(SPTE_ver.^2 + SPTE_hor.^2);

trlerrors.screen.ver = ver_mean_targ - ver_mean_stop ;
trlerrors.screen.hor = hor_mean_targ - hor_mean_stop ;
trlerrors.screen.all = sqrt(trlerrors.screen.ver.^2 + trlerrors.screen.hor.^2);
trlerrors.arena.all = tracking.misc.behverrors.all.val.*ones(ntrls,Nt);

% steering error in screen coordinates
[tracking_dist2end.eyepos.errcorr.targ.screen.r,tracking_dist2end.eyepos.errcorr.targ.screen.p] = ...
       arrayfun(@(t) nancorr(trlerrors.screen.all(:,t),TTE_all(:,t)), 1:Nt);
[tracking_dist2end.eyepos.errcorr.stop.screen.r,tracking_dist2end.eyepos.errcorr.stop.screen.p] = ...
       arrayfun(@(t) nancorr(trlerrors.screen.all(:,t),SPTE_all(:,t)), 1:Nt);
% steering error in arena coordinates
[tracking_dist2end.eyepos.errcorr.targ.arena.r,tracking_dist2end.eyepos.errcorr.targ.arena.p] = ...
       arrayfun(@(t) nancorr(trlerrors.arena.all(:,t),TTE_all(:,t)), 1:Nt);
[tracking_dist2end.eyepos.errcorr.stop.arena.r,tracking_dist2end.eyepos.errcorr.stop.arena.p] = ...
       arrayfun(@(t) nancorr(trlerrors.arena.all(:,t),SPTE_all(:,t)), 1:Nt);
    
tracking_dist2end.eyepos.error.targ.hor = TTE_hor;
tracking_dist2end.eyepos.error.targ.ver = TTE_ver;
tracking_dist2end.eyepos.error.targ.all = TTE_all;
tracking_dist2end.eyepos.error.stop.hor = SPTE_hor;
tracking_dist2end.eyepos.error.stop.ver = SPTE_ver;
tracking_dist2end.eyepos.error.stop.all = SPTE_all;

tracking_dist2end.eyepos.error.behv_arena.all = trlerrors.arena.all;
tracking_dist2end.eyepos.error.behv_screen.ver = trlerrors.screen.ver;
tracking_dist2end.eyepos.error.behv_screen.hor = trlerrors.screen.hor;    

%% Save new positions
tracking_dist2end.tarpos.screen.ver_mean = ver_mean_targ;
tracking_dist2end.tarpos.screen.hor_mean = hor_mean_targ;
tracking_dist2end.stopos.screen.ver_mean = ver_mean_stop;
tracking_dist2end.stopos.screen.hor_mean = hor_mean_stop;
tracking_dist2end.eyepos.screen.ver_mean = ver_mean;
tracking_dist2end.eyepos.screen.hor_mean = hor_mean;

tracking_dist2end.misc.ts = ds(:);
tracking_dist2end.misc.tau.val = tracking.misc.tau.val;
tracking_dist2end.misc.behverrors.all.val = tracking.misc.behverrors.all.val;
tracking_dist2end.misc.behverrors.r.val = tracking.misc.behverrors.r.val;
tracking_dist2end.misc.behverrors.th.val = tracking.misc.behverrors.th.val;
tracking_dist2end.misc.sampleflag = 'distance2end';
tracking_dist2end.misc.rmsac = tracking.misc.rmsac;
%% save saccadic eye movements
tracking_dist2end.saccade.true.val = tracking.saccade.true.val;
tracking_dist2end.saccade.true.dir = tracking.saccade.true.dir;
tracking_dist2end.saccade.pred_targ.val = tracking.saccade.pred_targ.val;
tracking_dist2end.saccade.pred_targ.dir = tracking.saccade.pred_targ.dir;
tracking_dist2end.saccade.pred_stop.val = tracking.saccade.pred_stop.val;
tracking_dist2end.saccade.pred_stop.dir = tracking.saccade.pred_stop.dir;
tracking_dist2end.saccade.trl = tracking.saccade.trl;

tracking_dist2end.saccade.time  = sacdist;

%% compare predicted & eye positions from all trials aligned to target fixation
% target position
ver_mean_targ1 = (ver_mean_targ)';
hor_mean_targ1 = (hor_mean_targ)';
% stopping position
ver_mean_stop1 = (ver_mean_stop)';
hor_mean_stop1 = (hor_mean_stop)';
% true eye position
ver_mean1 = (ver_mean)';
hor_mean1 = (hor_mean)';
% cross-correlation between observed and predicted eye positions
zeropad_length = 200; % pad these many bins around each trial
maxlag = 400;
% vertical
obs = [zeros(zeropad_length,ntrls); ver_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;

targ = [zeros(zeropad_length,ntrls); ver_mean_targ1; zeros(zeropad_length,ntrls)]; 
targ_shuffled = targ(:,randperm(ntrls));
targ = targ(:); targ(isnan(targ)) = 0; targ_shuffled = targ_shuffled(:); targ_shuffled(isnan(targ_shuffled)) = 0;
[tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.xcorr.c, tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.xcorr.lag] = xcorr(targ,obs,maxlag,'coeff');
tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.xcorr.c_shuffled = xcorr(targ_shuffled,obs,maxlag,'coeff');

stop = [zeros(zeropad_length,ntrls); ver_mean_stop1; zeros(zeropad_length,ntrls)]; 
stop_shuffled = stop(:,randperm(ntrls));
stop = stop(:); stop(isnan(stop)) = 0; stop_shuffled = stop_shuffled(:); stop_shuffled(isnan(stop_shuffled)) = 0;
[tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.xcorr.c, tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.xcorr.lag] = xcorr(stop,obs,maxlag,'coeff');
tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.xcorr.c_shuffled = xcorr(stop_shuffled,obs,maxlag,'coeff');
% horizontal
obs = [zeros(zeropad_length,ntrls); hor_mean1; zeros(zeropad_length,ntrls)];  
obs = obs(:); obs(isnan(obs)) = 0;

targ = [zeros(zeropad_length,ntrls); hor_mean_targ1; zeros(zeropad_length,ntrls)]; 
targ_shuffled = targ(:,randperm(ntrls));
targ = targ(:); targ(isnan(targ)) = 0; targ_shuffled = targ_shuffled(:); targ_shuffled(isnan(targ_shuffled)) = 0;
[tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.xcorr.c, tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.xcorr.lag] = xcorr(targ,obs,maxlag,'coeff');
tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.xcorr.c_shuffled = xcorr(targ_shuffled,obs,maxlag,'coeff');

stop = [zeros(zeropad_length,ntrls); hor_mean_stop1; zeros(zeropad_length,ntrls)]; 
stop_shuffled = stop(:,randperm(ntrls));
stop = stop(:); stop(isnan(stop)) = 0; stop_shuffled = stop_shuffled(:); stop_shuffled(isnan(stop_shuffled)) = 0;
[tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.xcorr.c, tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.xcorr.lag] = xcorr(stop,obs,maxlag,'coeff');
tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.xcorr.c_shuffled = xcorr(stop_shuffled,obs,maxlag,'coeff');

% timecourse of component-wise corr between predicted & true eye position
[tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.rho.startaligned,tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.pval.startaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_targ1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.rho.startaligned,tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.pval.startaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_targ1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.rho.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.pval.startaligned] = arrayfun(@(i) corr(ver_mean1(i,:)',ver_mean_stop1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.rho.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.pval.startaligned] = arrayfun(@(i) corr(hor_mean1(i,:)',hor_mean_stop1(i,:)','Type','Spearman','rows','complete'), 1:Nt);
% timecourse of component-wise regression between predicted & true eye position
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.beta.startaligned,tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.beta.startaligned,tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.beta.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.beta.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.int.startaligned] = deal(beta(1,:),beta(2,:));
% timecourse of component-wise multiple regression between target/stop & true eye position
multi = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ1(i,:)' ver_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.multi.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.multi.startaligned] = deal(multi(1,:),multi(2,:));
multi = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ1(i,:)' hor_mean_stop1(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.multi.startaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.multi.startaligned] = deal(multi(1,:),multi(2,:));

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
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.mu.startaligned = nanmean(var_explained_targ);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.sem.startaligned = nanstd(var_explained_targ);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained_shuffled.mu.startaligned = nanmean(var_explained_targ_shuffled);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained_shuffled.sem.startaligned = nanstd(var_explained_targ_shuffled);

    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.mu.startaligned = nanmean(var_explained_stop);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.sem.startaligned = nanstd(var_explained_stop);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained_shuffled.mu.startaligned = nanmean(var_explained_stop_shuffled);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained_shuffled.sem.startaligned = nanstd(var_explained_stop_shuffled);

else
    % target position
    squared_err = sum(nanmean((obs - targ).^2,3));
    var_pred = sum(nanvar(targ,[],3));
    var_explained_targ = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.mu.startaligned = var_explained_targ;
    % stopping position
    squared_err = sum(nanmean((obs - stop).^2,3));
    var_pred = sum(nanvar(stop,[],3));
    var_explained_stop = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.mu.startaligned = var_explained_stop;
end
%% compare predicted & observed eye positions from all trials aligned to END of movement
% target position
ver_mean_targ2 = fliplr(ver_mean_targ)';
hor_mean_targ2 = fliplr(hor_mean_targ)';
% stopping position
ver_mean_stop2 = fliplr(ver_mean_stop)';
hor_mean_stop2 = fliplr(hor_mean_stop)';
% true eye position
ver_mean2 = fliplr(ver_mean)';
hor_mean2 = fliplr(hor_mean)';

% timecourse of component-wise corr between predicted & true eye position
[tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.rho.stopaligned,tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.pval.stopaligned] = arrayfun(@(i) corr(ver_mean2(i,:)',ver_mean_targ2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.rho.stopaligned,tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.pval.stopaligned] = arrayfun(@(i) corr(hor_mean2(i,:)',hor_mean_targ2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.rho.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.pval.stopaligned] = arrayfun(@(i) corr(ver_mean2(i,:)',ver_mean_stop2(i,:)','Type','Spearman','rows','complete'), 1:Nt);
[tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.rho.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.pval.stopaligned] = arrayfun(@(i) corr(hor_mean2(i,:)',hor_mean_stop2(i,:)','Type','Spearman','rows','complete'), 1:Nt);

% component-wise regression between predicted & true eye position
beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_targ2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.beta.stopaligned,tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_targ2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.beta.stopaligned,tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(ver_mean2(i,:)',[ver_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.beta.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));
beta = cell2mat(arrayfun(@(i) regress(hor_mean2(i,:)',[hor_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.beta.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.int.stopaligned] = deal(beta(1,:),beta(2,:));

% timecourse of component-wise multiple regression between target/stop & true eye position
multi = cell2mat(arrayfun(@(i) regress(ver_mean1(i,:)',[ver_mean_targ2(i,:)' ver_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.ver_mean.multi.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.ver_mean.multi.stopaligned] = deal(multi(1,:),multi(2,:));
multi = cell2mat(arrayfun(@(i) regress(hor_mean1(i,:)',[hor_mean_targ2(i,:)' hor_mean_stop2(i,:)' zeros(ntrls,1)]), 1:Nt, 'UniformOutput', false)); [tracking_dist2end.eyepos.pred_vs_true.targ.hor_mean.multi.stopaligned,tracking_dist2end.eyepos.pred_vs_true.stop.hor_mean.multi.stopaligned] = deal(multi(1,:),multi(2,:));

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
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.mu.stopaligned = nanmean(var_explained_targ);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.sem.stopaligned = nanstd(var_explained_targ);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained_shuffled.mu.stopaligned = nanmean(var_explained_targ_shuffled);
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained_shuffled.sem.stopaligned = nanstd(var_explained_targ_shuffled);

    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.mu.stopaligned = nanmean(var_explained_stop);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.sem.stopaligned = nanstd(var_explained_stop);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained_shuffled.mu.stopaligned = nanmean(var_explained_stop_shuffled);
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained_shuffled.sem.stopaligned = nanstd(var_explained_stop_shuffled);
else
    squared_err = sum(nanmean((obs - targ).^2,3));
    var_pred = sum(nanvar(targ,[],3));
    var_explained_targ = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking_dist2end.eyepos.pred_vs_true.targ.var_explained.mu.stopaligned = var_explained_targ;

    squared_err = sum(nanmean((obs - stop).^2,3));
    var_pred = sum(nanvar(stop,[],3));
    var_explained_stop = 1 - (squared_err(:)./var_pred(:)); % try taking sqrt or cosine of var_explained
    tracking_dist2end.eyepos.pred_vs_true.stop.var_explained.mu.stopaligned = var_explained_stop;

end

