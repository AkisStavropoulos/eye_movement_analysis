%%
%%
%% EYE MOVEMENTS

%%
%%
%% Load data
inp = input('Clear everything? 1 yes, 0 no: ');
if inp; clear; end
data_folder = 'C:\Users\ges6\Documents\MATLAB\Data\Single firefly with motion cuing\Subjects\New\';
cd(data_folder);
folders = dir(data_folder);   folders = folders([folders.isdir]);

indx = cell2mat(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders,'uniformoutput',false));
subject_name = {folders(indx).name};

% Load data
savenow = 0;
if savenow
    [subject,button_push] = SaveStats(subject,data_folder,savenow);
else
    fldcont = dir(data_folder);
    dataindx = arrayfun(@(x) contains(x.name,'EyeAnalysisData.mat'), fldcont);
    if any(dataindx)
        fixeye = 0;
        tic;
        disp('.... Loading EyeAnalysisData.mat')
        load(fldcont(dataindx).name);
        toc;
    else
    tic;
    [subject_backup,~] = LoadData(data_folder,subject_name);
    fixeye = 1;
    toc;
    
    default_prs;
    % name subjects
    subj_names = dir(data_folder);
    subj_names = subj_names(~ismember({subj_names.name},{'.','..', '.DS_Store'}));
    dir_indx = arrayfun(@(x) x.isdir==1,subj_names);
    subj_names = subj_names(dir_indx);
    for k = 1:numel(subject_backup); subject_backup(k).name = subj_names(k).name; end

    end
end


% Melissa Fix
subject = subject_backup;
melissa_indx = arrayfun(@(x) strcmp(x.name,'Melissa'),subject);

if any(melissa_indx)
    Ntrials = numel(subject(melissa_indx).trials);
    mcnames = fieldnames(subject(melissa_indx).trials(1).mc);
    
    for i = 1:Ntrials
        if isempty(subject(melissa_indx).trials(i).mc.JS_X_Raw)
            trlength = numel(subject(melissa_indx).trials(i).continuous.ts);
            
            for k = 1:numel(mcnames)
                subject(melissa_indx).trials(i).mc.(mcnames{k}) = nan(trlength,1);
            end
        end
    end
    disp(['....Fixed Melissa MC data'])
end

% remove subjects
subject = subject_backup;
if strcmp(subject_backup(5).name,'Ding'); subject_backup(5) = []; end
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
% keepindx = [1 3 4 5 8 9 10 11 13];
keepindx = [1 2 7 9 10 12 13 14]; % the correct one
subject = subject_backup(keepindx);
subject_name = subject_name(keepindx);

%% Fix vertical eye position slip

if fixeye
thresh = 20;
subject = subject_backup;
Nsubs = numel(subject);

figure;
for i = 1:Nsubs
    clf;
    init_eye_hor = abs(arrayfun(@(x) 0.5*(x.continuous.yle(60)+x.continuous.yre(60)), subject(i).trials));
    init_eye_ver = arrayfun(@(x) 0.5*(x.continuous.zle(60)+x.continuous.zre(60)), subject(i).trials);
            
    % compute actual target position and add to correction
    delta = prs.interoculardist/2;
    zt = -prs.height;
    xt = arrayfun(@(x) (x.prs.fireflyposx-x.continuous.xmp(60)), subject(i).trials);
    yt = arrayfun(@(x) (x.prs.fireflyposy-x.continuous.ymp(60)), subject(i).trials);
    
    [yle_targ, zle_targ, yre_targ, zre_targ] = world2eye(xt,yt,zt,delta);
    init_tar_ver = nanmean([zle_targ' , zre_targ'],2); 
    init_tar_hor = nanmean([yle_targ' , yre_targ'],2); 
    tar_ver_mean = mean(init_tar_ver);
    tar_ver_std = std(init_tar_ver);

    
    eye_ver = arrayfun(@(x) 0.5*(x.continuous.zle+x.continuous.zre), subject(i).trials,'un',0);
    
    subplot(2,1,1); hold on;
    plot(init_eye_ver); 
    plot(init_tar_ver); 
    ylim([-40 30]);
    title(subject(i).name);
    
    subplot(2,1,2); hold on;
    cellfun(@plot,eye_ver); hline(0,'k'); ylim([-80 50])
    
    % fix
    Ntrials = numel(init_eye_ver);
    zeroindx = init_eye_ver==0;
    offset_ver = init_eye_ver - tar_ver_mean;
    final_eye_ver = init_eye_ver - offset_ver + randn(1,Ntrials)*tar_ver_std;
    final_eye_ver(zeroindx) = 0; 
    offset2 = final_eye_ver - init_eye_ver;
    ver_targ_dev = init_tar_ver - tar_ver_mean;
    init_eye_ver = init_eye_ver + offset2;
    
    for j = 1:Ntrials
        if ~any(j==find(zeroindx))
            subject(i).trials(j).continuous.zle = subject(i).trials(j).continuous.zle + offset2(j) + ver_targ_dev(j);
            subject(i).trials(j).continuous.zre = subject(i).trials(j).continuous.zre + offset2(j) + ver_targ_dev(j);
        end
    end
    
    init_eye_ver = arrayfun(@(x) 0.5*(x.continuous.zle(60)+x.continuous.zre(60)), subject(i).trials);
    subplot(2,1,1); hold on;
    plot(init_eye_ver)
end

% change subject_backup too
for i = 1:Nsubs
    ind = arrayfun(@(x) strcmp(x.name,subject(i).name),subject_backup);
    if any(ind)
        subject_backup(ind) = subject(i);
    end
end

save('EyeAnalysisData','subject_backup','-v7.3');
disp('.........saved corrected vertical eye positions')
end

%%
%% BEHAVIORAL PERFORMANCE FIGURES
%%
%% Bias Detection and Multiplicative Model Fit (+ compare 1st-button-push bias with end-of-trial bias)
keepindx = [1 2 7 9 10 12 13 14];
subject = subject_backup(keepindx);

polyorder = 1;  plt = 0;    intercept = 0;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params); 
[bias.r,bias.th] = ErrScatterFit(subject,params,polyorder,intercept,plt);
Nstim = size(poolindx,2);

% Plot angular vs radial bias
count = 1;
% plot biases
if ~(size(bias.r,2)>3); colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
else
colr_temp1 = brewermap(9,'*BrBG'); 
colr_temp2 = brewermap(9,'*YlOrBr');
colr_temp3 = brewermap(7,'*RdPu');
colr = [colr_temp1(1:3,:) ; colr_temp2(1:3,:) ; colr_temp3(1:3,:)];
end
figure; hold on;    axis equal;
xlim([0 1]);ylim([0 1]);    hline(1,'k--'); vline(1,'k--');
plot(0:2,0:2,'k','HandleVisibility','off'); xlim([0 1.5]);  ylim([0 1.5]);

conf_int = .68;
for s = 1:size(bias.r,2)
    [Xellipse,Yellipse,x0,y0] = error_ellipse(bias.r(:,s),bias.th(:,s),conf_int,0);
    
    plot(bias.r(:,s),bias.th(:,s),'x','Color',colr(s,:),'MarkerSize',5);
    plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
    plot(x0, y0,'d','Color',colr(s,:),'MarkerFaceColor',colr(s,:),'MarkerSize',15,'HandleVisibility','off');
    if mod(s,3)==0; legend(legend_input{1+3*(count-1):s},'Location','northwest'); count=count+1;  end
end
if ~(size(bias.r,2)>3); [~,legend_input] = get_poolindx(subject,[]); legend(legend_input{:},'Location','northwest'); end  
xlabel('Radial Bias');ylabel('Angular Bias'); 

% F-test
regress_ftest;

%% Radial vs angular bias for different TAU BINS for each modality

params = 'stimtau';
[mu,sem,bias] = bias_taugroups(subject,params,0);

% check significance across subjects
% t-test
% vestibular
[h_r.ves(1),p_r.ves(1)] = ttest(bias.r(:,1),bias.r(:,2));   [h_th.ves(1),p_th.ves(1)] = ttest(bias.th(:,1),bias.th(:,2));
[h_r.ves(2),p_r.ves(2)] = ttest(bias.r(:,2),bias.r(:,3));   [h_th.ves(2),p_th.ves(2)] = ttest(bias.th(:,2),bias.th(:,3));
[h_r.ves(3),p_r.ves(3)] = ttest(bias.r(:,1),bias.r(:,3));   [h_th.ves(3),p_th.ves(3)] = ttest(bias.th(:,1),bias.th(:,3));
% visual
[h_r.vis(1),p_r.vis(1)] = ttest(bias.r(:,4),bias.r(:,5));   [h_th.vis(1),p_th.vis(1)] = ttest(bias.th(:,4),bias.th(:,5));
[h_r.vis(2),p_r.vis(2)] = ttest(bias.r(:,5),bias.r(:,6));   [h_th.vis(2),p_th.vis(2)] = ttest(bias.th(:,5),bias.th(:,6));
[h_r.vis(3),p_r.vis(3)] = ttest(bias.r(:,4),bias.r(:,6));   [h_th.vis(3),p_th.vis(3)] = ttest(bias.th(:,4),bias.th(:,6));
% combined
[h_r.comb(1),p_r.comb(1)] = ttest(bias.r(:,7),bias.r(:,8)); [h_th.comb(1),p_th.comb(1)] = ttest(bias.th(:,7),bias.th(:,8));
[h_r.comb(2),p_r.comb(2)] = ttest(bias.r(:,8),bias.r(:,9)); [h_th.comb(2),p_th.comb(2)] = ttest(bias.th(:,8),bias.th(:,9));
[h_r.comb(3),p_r.comb(3)] = ttest(bias.r(:,7),bias.r(:,9)); [h_th.comb(3),p_th.comb(3)] = ttest(bias.th(:,7),bias.th(:,9));

% Signed rank
% vestibular
[p_r.ves(1),h_r.ves(1)] = signrank(bias.r(:,1),bias.r(:,2));   [p_th.ves(1),h_th.ves(1)] = signrank(bias.th(:,1),bias.th(:,2));
[p_r.ves(2),h_r.ves(2)] = signrank(bias.r(:,2),bias.r(:,3));   [p_th.ves(2),h_th.ves(2)] = signrank(bias.th(:,2),bias.th(:,3));
[p_r.ves(3),h_r.ves(3)] = signrank(bias.r(:,1),bias.r(:,3));   [p_th.ves(3),h_th.ves(3)] = signrank(bias.th(:,1),bias.th(:,3));
% visual
[p_r.vis(1),h_r.vis(1)] = signrank(bias.r(:,4),bias.r(:,5));   [p_th.vis(1),h_th.vis(1)] = signrank(bias.th(:,4),bias.th(:,5));
[p_r.vis(2),h_r.vis(2)] = signrank(bias.r(:,5),bias.r(:,6));   [p_th.vis(2),h_th.vis(2)] = signrank(bias.th(:,5),bias.th(:,6));
[p_r.vis(3),h_r.vis(3)] = signrank(bias.r(:,4),bias.r(:,6));   [p_th.vis(3),h_th.vis(3)] = signrank(bias.th(:,4),bias.th(:,6));
% combined
[p_r.comb(1),h_r.comb(1)] = signrank(bias.r(:,7),bias.r(:,8)); [p_th.comb(1),h_th.comb(1)] = signrank(bias.th(:,7),bias.th(:,8));
[p_r.comb(2),h_r.comb(2)] = signrank(bias.r(:,8),bias.r(:,9)); [p_th.comb(2),h_th.comb(2)] = signrank(bias.th(:,8),bias.th(:,9));
[p_r.comb(3),h_r.comb(3)] = signrank(bias.r(:,7),bias.r(:,9)); [p_th.comb(3),h_th.comb(3)] = signrank(bias.th(:,7),bias.th(:,9));

%% Response variability and R^2 across sensory modalities (measure of uncertainty)
keepindx = [1 2 7 9 10 12 13 14];
subject = subject_backup(keepindx);
Nsubs = numel(subject);

params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nstim = size(poolindx,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

var_r = []; var_th = []; R2_r = []; R2_th = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlindx = poolindx{i,s};
        
        r_tar = arrayfun(@(x) x.prs.r_tar,subject(i).trials(trlindx));
        th_tar = arrayfun(@(x) x.prs.th_tar,subject(i).trials(trlindx));
        r_sub = arrayfun(@(x) x.prs.r_sub,subject(i).trials(trlindx));
        th_sub = arrayfun(@(x) x.prs.th_sub,subject(i).trials(trlindx));
        
        % radial component
        y_r = r_sub(:);
        X_r = [zeros(numel(r_sub),1) r_tar(:)];
        [~,~,R_r,~,stats_r] = regress(y_r,X_r);
        var_r(i,s) = stats_r(end)/10000; % convert to m^2
        R2_r(i,s) = stats_r(1);
        
        % angular component
        y_th = th_sub(:);
        X_th = [zeros(numel(th_sub),1) th_tar(:)];
        [~,~,R_th,~,stats_th] = regress(y_th,X_th);
        var_th(i,s) = stats_th(end);
        R2_th(i,s) = stats_th(1);
    end
end

% plot variability
figure; 
for s = 1:Nstim
    subplot(1,2,1); hold on;
    bar(s,nanmean(var_r(:,s)),'facecolor',colr(s,:));errorbar(s,nanmean(var_r(:,s)),nanstd(var_r(:,s))./sqrt(Nsubs),'k','capsize',0);
    title('radial variance'); ylabel('variance [m^2]'); xticks(1:Nstim); xticklabels(legend_input);
        
    subplot(1,2,2); hold on;
    bar(s,nanmean(var_th(:,s)),'facecolor',colr(s,:));errorbar(s,nanmean(var_th(:,s)),nanstd(var_th(:,s))./sqrt(Nsubs),'k','capsize',0);
    title('angular variance'); ylabel('variance [deg^2]'); xticks(1:Nstim); xticklabels(legend_input);
end
[h_r,p_r] = ttest(var_r(:,1),var_r(:,2));

[h_th,p_th] = ttest(var_th(:,1),var_th(:,2));

% plot R^2
figure; 
for s = 1:Nstim
    subplot(1,2,1); hold on;
    bar(s,nanmean(R2_r(:,s)),'facecolor',colr(s,:));errorbar(s,nanmean(R2_r(:,s)),nanstd(R2_r(:,s))./sqrt(Nsubs),'k','capsize',0);
    title('radial R^2'); ylabel('R-squared'); xticks(1:Nstim); xticklabels(legend_input);
        
    subplot(1,2,2); hold on;
    bar(s,nanmean(R2_th(:,s)),'facecolor',colr(s,:));errorbar(s,nanmean(R2_th(:,s)),nanstd(R2_th(:,s))./sqrt(Nsubs),'k','capsize',0);
    title('angular R^2'); ylabel('R-squared'); xticks(1:Nstim); xticklabels(legend_input);
end

figure; hold on;
for s = 1:Nstim
    R2_avg = 0.5*(R2_r(:,s)+R2_th(:,s));
    bar(s,nanmean(R2_avg),'facecolor',colr(s,:));errorbar(s,nanmean(R2_avg),nanstd(R2_avg)./sqrt(Nsubs),'k','capsize',0);
    title('average R^2'); ylabel('R-squared'); xticks(1:Nstim); xticklabels(legend_input);
end


%%
%%
%% EYE MOVEMENT ANALYSIS
%%
%%
%% Trajectories of eye movements (screen coordinates)
% smooth signal with gaussian filter to look smoother
params = 'stimtype';
% subject = subject_backup(setdiff(1:end,[2 5 13 15]));
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
dt = 1/60;
si = 1/dt;  ei = 12/dt;
SDthresh = 1;
sac_pad = .3; % seconds
[b_r,b_th] = ErrScatterFit(subject,params,1,0,0);
sac_pad = .3; % seconds
for i = 1:length(subject)
    figure;
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        Zeye = [];  Yeye = [];  Zstart = [];    Ystart = []; Zend = []; Yend = [];   Tend = []; cnt = 0;
        for j = 1:length(indx)-300
            
            if sum(subject(i).trials(indx(j)).continuous.zle) ~=0
            Zeye{j} = (subject(i).trials(indx(j)).continuous.zle + subject(i).trials(indx(j)).continuous.zre)./2;
            Yeye{j} = (subject(i).trials(indx(j)).continuous.yle + subject(i).trials(indx(j)).continuous.yre)./2;
            % Remove blinks
            Zeye{j}(subject(i).trials(indx(j)).continuous.blink) = nan; Yeye{j}(subject(i).trials(indx(j)).continuous.blink) = nan;
            % Remove saccades
            t_sac = subject(i).trials(indx(j)).events.t_sac;
            for t = 1:length(t_sac)
                Zeye{j}(t_sac(t)/dt:t_sac(t)/dt + sac_pad/dt) = nan;    Yeye{j}(t_sac(t)/dt:t_sac(t)/dt + sac_pad/dt) = nan;
            end
            % Check length
            if length(Zeye{j}) < ei
                ind = length(Zeye{j});
            else 
                ind = ei;
            end
            Zeye{j} = Zeye{j}(si:ind);   Zstart(j) = Zeye{j}(1);     Zend(j) = Zeye{j}(end);
            Yeye{j} = Yeye{j}(si:ind);   Ystart(j) = Yeye{j}(1);     Yend(j) = Yeye{j}(end);
            subplot(1,size(poolindx,2),s);  hold on;
            plot(Yeye{j},Zeye{j},'Color',colr(s,:),'LineWidth',.5); 
            end
        end
        plot(Ystart,Zstart,'b.');   plot(Yend,Zend,'k+');
        xlabel('Eye position Y');   ylabel('Eye position Z');   legend(legend_input{s}); 
        axis equal;                 xlim([-40 40]);ylim([-60 20]);  
    end
    sgtitle(subject(i).name);
end
%% Find example trial for target-tracking in Eye coordinates
for i = 1:length(subject)
    figure;
    
   for s = 2:size(poolindx,2) 
       indx = 1:length(eye_movement{i,s}.eyepos.behvcorr.eye_err);
       
       avg_err = eye_movement{i,s}.eyepos.behvcorr.eye_err;
       trl_err = eye_movement{i,s}.eyepos.behvcorr.eye_err_ts.all;
       ts = eye_movement{i,s}.eyepos.behvcorr.ts;
       target_off = eye_movement{i,s}.eyepos.target_off.indx;
       
       for j = [217,220] %120:length(indx)
           plot(eye_movement{i,s}.eyepos.pred.hor_mean.val{j},eye_movement{i,s}.eyepos.pred.ver_mean.val{j},'r--'); hold on;
           plot(eye_movement{i,s}.eyepos.pred.hor_mean.val{j}(1),eye_movement{i,s}.eyepos.pred.ver_mean.val{j}(1),'ro','markersize',12);
           plot(eye_movement{i,s}.eyepos.true.hor_mean.val{j},eye_movement{i,s}.eyepos.true.ver_mean.val{j},'k');
           plot(eye_movement{i,s}.eyepos.true.hor_mean.val{j}(1),eye_movement{i,s}.eyepos.true.ver_mean.val{j}(1),'ko','markersize',12);
           xlabel('HOR [deg]'); ylabel('VER [deg]'); title([subject(i).name ' - Trial ' num2str(j)]);
           axis([-25 25 -50 0])
           hold off;
           
           if 1
           hold on; colr = jet(length(eye_movement{i,s}.eyepos.true.hor_mean.val{j}));
           for n = 1:length(eye_movement{i,s}.eyepos.true.hor_mean.val{j})
               plot(eye_movement{i,s}.eyepos.true.hor_mean.val{j}(n),eye_movement{i,s}.eyepos.true.ver_mean.val{j}(n),'.','color',colr(n,:));
               plot(eye_movement{i,s}.eyepos.pred.hor_mean.val{j}(n),eye_movement{i,s}.eyepos.pred.ver_mean.val{j}(n),'o','markeredgecolor',colr(n,:)*0.5);
           end
           colormap(colr); colorbar;
           c = colorbar;
           c.Label.String = 'time [s]';
           
           figure;
           plot(ts,trl_err(:,j));hold on;
           for  n = 1:length(eye_movement{i,s}.eyepos.true.hor_mean.val{j})
               plot(ts(n),trl_err(n,j),'.','color',colr(n,:));
           end
           xlim([0 12]); ylim([0 60]); vline(target_off(j)*dt,'k--'); hline(avg_err(j));
           xlabel('time [s]'); ylabel('Eye-tracking error [deg]'); title([subject(i).name ' - Trial ' num2str(j)]);
           end

       end
   end
end
%% Predicted vs Observed eye position (when target goes OFF)
params = [];
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
default_prs; upperbnd = 0; BELIEF = 0; 
prs.boots = 0; 
[poolindx,legend_input] = get_poolindx(subject,params);
if 0
[~,eye_movementALL] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
end
plt = 0;
rho = []; pval = [];
Nsamples = 1; % nsamples
dt = 1/60;
for i = 1:length(subject) % [6 12 14]
    if plt; figure; end
    for s = 1:size(poolindx,2)
        N = length(eye_movementALL{i,s}.eyepos.pred.ver_mean.val);
        PredPosZ = [];  PredPosY = [];  TruePosZ = [];  TruePosY = [];
        for j = 1:N
            % target off index
            ind = eye_movementALL{i,s}.eyepos.target_off.indx(j);
            if ind < 1;    ind = 1;    end
            PredPosZ(j) = eye_movementALL{i,s}.eyepos.pred.ver_mean.val{j}(ind);  PredPosY(j) = eye_movementALL{i,s}.eyepos.pred.hor_mean.val{j}(ind);
            TruePosZ(j) = eye_movementALL{i,s}.eyepos.true.ver_mean.val{j}(ind);  TruePosY(j) = eye_movementALL{i,s}.eyepos.true.hor_mean.val{j}(ind);
%             plot(PredPosZ,'b--');  hold on;   plot(TruePosZ,'b');  plot(PredPosY,'k--');  plot(TruePosY,'k');   hold off;
        end
        [rho.z(i,s),pval.z(i,s)] = nancorr(PredPosZ(:),TruePosZ(:));
        [rho.y(i,s),pval.y(i,s)] = nancorr(PredPosY(:),TruePosY(:));       
        if plt
        subplot(2,size(poolindx,2),s); hold on;   axis equal;   xlim([-25 5]);  ylim([-25 5]);   
        plot(PredPosZ,TruePosZ,'.','Color',colr(s,:),'MarkerSize',3);  
        plot(-30:30,-30:30,'k--'); xlabel('Predicted [deg]');  ylabel('Observed [deg]');
        title('Vertical');    legend(legend_input{s});
        vline(0,'k');    hline(0,'k');   
        
        subplot(2,size(poolindx,2),s+size(poolindx,2)); hold on;  axis equal;   xlim([-40 40]);  ylim([-40 40]);    
        plot(PredPosY,TruePosY,'.','Color',colr(s,:),'MarkerSize',3);  
        plot(-30:30,-30:30,'k--');  xlabel('Predicted [deg]');  ylabel('Observed [deg]');
        title('Horizontal');    vline(0,'k');    hline(0,'k');
        end
    end
    if plt;  sgtitle(subject(i).name)  ; end
end
rho.y(rho.y==0) = []; rho.z(rho.z==0) = [];
version_mu = mean(rho.y)
version_sd = std(rho.y)
elevation_mu = mean(rho.z)
elevation_sd = std(rho.z)
%% Time course of lateral version and elevation
            % Lateral version converges to 0 faster for vestibular,
            % consistent with undershooting angular bias
sac_pad = .3; % seconds
si = 1/dt;  ei = 12/dt;
SDthresh = 1;
for i = [1 6 7 15]
    figure;
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        Zeye = [];  Yeye = [];  Zstart = [];    Ystart = []; Zend = []; Yend = [];   Tend = [];     Jindx = [];
        for j = 1:length(indx)-300   
            Zeye{j} = (subject(i).trials(indx(j)).continuous.zle + subject(i).trials(indx(j)).continuous.zre)./2;
            Yeye{j} = (subject(i).trials(indx(j)).continuous.yle + subject(i).trials(indx(j)).continuous.yre)./2;
            % Remove blinks
            Zeye{j}(subject(i).trials(indx(j)).continuous.blink) = nan; Yeye{j}(subject(i).trials(indx(j)).continuous.blink) = nan;
            % Remove saccades
            t_sac = subject(i).trials(indx(j)).events.t_sac;
            for t = 1:length(t_sac)
                Zeye{j}(t_sac(t)/dt:t_sac(t)/dt + sac_pad/dt) = nan;    Yeye{j}(t_sac(t)/dt:t_sac(t)/dt + sac_pad/dt) = nan;
            end
            % Check length
            if length(Zeye{j}) < ei
                ind = length(Zeye{j});
            else 
                ind = ei;
            end
            if ~isempty(subject(i).trials(indx(j)).mc.JS_X_Raw) && any(~isnan(subject(i).trials(indx(j)).mc.JS_X_Raw))
                SteerEnd = find(subject(i).trials(indx(j)).mc.JS_X_Raw,1,'last');
                if SteerEnd < ind;  ind = SteerEnd;     end
            end
            Zeye{j} = Zeye{j}(si:ind);
            if isempty(Zeye{j})
                keyboard
            end
            Zstart(j) = Zeye{j}(1);     
            Zend(j) = Zeye{j}(end);    
            Tend(j) = length(Zeye{j});
            Yeye{j} = Yeye{j}(si:ind);   Ystart(j) = Yeye{j}(1);     Yend(j) = Yeye{j}(end);
            if ~(nanstd(Zeye{j}) <= SDthresh || nanstd(Yeye{j}) <= SDthresh) && ind*dt > 5
                Jindx = [Jindx j];
            subplot(2,size(poolindx,2),s);  hold on;
            plot((1:length(Yeye{j}))*dt,Yeye{j},'Color',colr(s,:),'LineWidth',.1);
            subplot(2,size(poolindx,2),s+3);  hold on;
            plot((1:length(Yeye{j}))*dt,Zeye{j},'Color',colr(s,:),'LineWidth',.1);
            end
        end
        subplot(2,size(poolindx,2),s);  hold on;  hline(0,'k');
        plot(Tend(Jindx)*dt,Yend(Jindx),'k.');   xlabel('time [s]'); ylabel('Lateral version');
        ylim([-40 40]); 
        subplot(2,size(poolindx,2),s+3);  hold on;  hline(0,'k')
        plot(Tend(Jindx)*dt,Zend(Jindx),'k.');   xlabel('time [s]'); ylabel('Elevation');
        ylim([-60 20]);
    end
    sgtitle(subject(i).name);
end

% correlation of observed and predicted eye position over time
default_prs; upperbnd = 0; BELIEF = 0; 
prs.boots = 0; 
params1 = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params1);
if 0
[eye_fixation,eye_movement] = GetEyeAnalysis(subject,prs,params1,upperbnd,BELIEF);
end

for i = 1:length(subject)
   for s = 1:size(poolindx,2) 
    indx = length(eye_movement{i,s}.eyepos.true.ver_mean.val);
       for j = 1:length(indx)
           obshor = eye_movement{i,s}.eyepos.true.hor_mean.val{j};
           obsver = eye_movement{i,s}.eyepos.true.ver_mean.val{j};
           predhor = eye_movement{i,s}.eyepos.pred.hor_mean.val{j};
           predver = eye_movement{i,s}.eyepos.pred.ver_mean.val{j};
           [rhohor,phor] = nancorr(predhor,obshor);
           [rhover,pver] = nancorr(predver,obsver);
       end
   end
end
%% Trial duration CDF
% Each subject separately
params = 'stimtype';
[trldurcdf,trldurs] = trl_duration_cdf(subject,params,plt);

%% Target-Tracking Index when target is OFF************make this a bar plot************
% ind = find(~isnan(eye_movement{i,s}.eyepos.true.ver_mean.val{j}),1);
default_prs;
params = [];
[poolindx,legend_input] = get_poolindx(subject,params);
upperbnd = 0; BELIEF = 0;
prs.boots = 0;

fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'eye_movementALL.mat'), fldcont);
if any(struct_indx)
    if ~exist('eye_movementALL')
        load(fldcont(struct_indx).name);
        disp('.....eye_movementALL Loaded.')
    else; disp('eye_movementALL already exists.')
    end
else
    [~,eye_movementALL] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
    save('eye_movementALL.mat','eye_movementALL','-v7.3')
end
varexp_start = [];  varexp_stop = [];
Subs = 1:length(subject);
% Subs([2 4 5 8 10 11 15]) = [];
subj_tick = [];
figure;hold on;
for i = 1:length(subject) % Subs % [2 5 8  15] are bad
    for s = 1:size(poolindx,2)
        % clip negative varexp to zero
        varexp_start{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained.mu.startaligned);  varexp_start{i,s}(varexp_start{i,s}<0) = 0;
%         sem_start{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained.sem.startaligned);
%         
%         
%         % plot
%         errorbar(i,sqrt(varexp_start{i,s}(60)),sem_start{i,s}(60),sem_start{i,s}(60),'ko',...
%             'Capsize',0,'MarkerFaceColor','k');
        xlabel('Subjects'); ylabel('Target-tracking index'); axis([0 length(subject)+1 0 1]);
        
    end
    
    legend(legend_input{:});
    subj_tick{i} = subject(i).name(1);
    
end
xticks(1:length(subject));    xticklabels([subj_tick]);
suptitle('TTI at target offset');

temp = cellfun(@(x) sqrt(x(60)), varexp_start);
TTI_mu = nanmean(temp);
TTI_sd = nanstd(temp);
%% TTI for modalities
subject = subject_backup;%(choose_subs); % subject_backup(setdiff(1:15,[5 8]));
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
upperbnd = 0;
BELIEF = 0;
rmsac = 0;
prs.boots = 0;

poolsubs = 0;
TTI = []; trlinfo = cell(1,size(poolindx,2));
for s = 1:size(poolindx,2)
    for i = 1:length(subject) % [4 5 8 11] are bad
        indx = poolindx{i,s};
        if upperbnd
            spatialstd = ComputeSpatialError(subject(i).trials(indx));
        else
            spatialstd = [];
        end
        if poolsubs; trlinfo{s} = [trlinfo{s}  subject(i).trials(indx)];
        else; TTI{i,s} = computeTTI(subject(i).trials(indx),spatialstd,prs,BELIEF,rmsac); end
    end
    disp(['......Subject = ' num2str(i)])
end
if poolsubs; TTI = cellfun(@(x) computeTTI(x,spatialstd,prs,BELIEF,rmsac),trlinfo,'uniformoutput',false); end

clear('trlinfo');
type = {'all','ver','hor'};
type = {'all'};
separate_label = 1;
trackdur = plot_TTI(TTI,[],[],type,prs,separate_label);

%% Plot TTE when target goes OFF
params = []; upperbnd = 0; BELIEF = 0;
% subject = subject_backup;
fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'eye_movementALL.mat'), fldcont);
if any(struct_indx)
    if ~exist('eye_movementALL')
        load(fldcont(struct_indx).name);
        disp('.....eye_movementALL Loaded.')
    else; disp('eye_movementALL already exists.')
    end
else
    [~,eye_movementALL] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
    save('eye_movementALL.mat','eye_movementALL','-v7.3')
end
dt = 1/60;
T = 1/dt;
wn = 2;
tte = [];
Nboots = 10000;

% eye_movementALL([4 5 10]) = [];

tte.all = cellfun(@(x) nanmean(x.eyepos.behvcorr.eye_err_ts.all(T-wn:T+wn,:)),eye_movementALL,'uniformoutput',false);
bstrpmat = cellfun(@(x) sort(randi(numel(x),numel(x),Nboots),1), tte.all,'uniformoutput',false);
tte.mu = cellfun(@(x,y) nanmean(nanmean(x(y))), tte.all, bstrpmat);
tte.se = cellfun(@(x,y) nanstd(nanmean(x(y))), tte.all, bstrpmat);

% figure; 
% errorbar(1:length(subject),tte.mu,tte.se,'ks','CapSize',0,'markerfacecolor','k');
% title('target-tracking error at target offset');
% xlabel('subjects'); ylabel('Target-tracking error [deg]'); axis([0 length(subject)+1 5 30]);

figure; hold on;
bar(1,nanmean(tte.mu),'facecolor',[.5 .5 .5]);
plot(1+0.1*randn(numel(tte.mu),1),tte.mu,'ko','markerfacecolor','k'); ylim([5 20]);
ylabel('Target-tracking error [deg]'); title('Target OFF');
%% Plot average TTE over time
params = 'stimtype';
subject = subject_backup(keepindx);
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(3,'Dark2');
cond2 = size(poolindx,2)/3;
if size(poolindx,2) > 3; colr = graded_colormap(colr,cond2); end
default_prs; upperbnd = 0; BELIEF = 0;
prs.boots = 0;
fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'eye_movement.mat'), fldcont);
if any(struct_indx)
    if ~exist('eye_movement')
        load(fldcont(struct_indx).name);
        disp('.....eye_movement Loaded.')
    else; disp('eye_movement already exists.')
    end
else
[~,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
save('eye_movement.mat','eye_movement','-v7.3')
end
dt = 1/60; 

trlength = cellfun(@(x) length(x.eyepos.behvcorr.ts),eye_movement);
minlength = min(trlength(:));
temp.mu = cellfun(@(x) nanmean(x.eyepos.behvcorr.eye_err_ts.all(1:minlength,:),2), eye_movement,'uniformoutput',false);
temp.sd = cellfun(@(x) nanstd(x.eyepos.behvcorr.eye_err_ts.all(1:minlength,:),[],2), eye_movement,'uniformoutput',false);

% pool subjects
% [temp.mu,temp.sd] = poolsubjects(temp.mu,temp.sd);
TTE = [];
figure; hold on;
for s = 1:size(poolindx,2)
    TTE.mu{s} = nanmean(temp.mu{s},2);
    TTE.se{s} = nanstd(temp.mu{s},[],2)./sqrt(length(subject));
    ts = (1:minlength)*dt;
%     if size(poolindx,2)>3; count = floor((s-1)./cond2)+1; else; count = 1; end
%     subplot(1,3,count); hold on;
    shadedErrorBar(ts-1,TTE.mu{s},TTE.se{s},'lineprops',{'linewidth',2,'color',colr(s,:)});
    xlabel('time since target offset [s]'); ylabel('Target-tracking error [deg]'); axis([0 8 5 30]);
%     legend(legend_input{(1:cond2)+cond2*(count-1)});
end

figure; hold on;
for s = 1:size(poolindx,2)
    ts = (1:minlength)*dt;
    subplot(1,size(poolindx,2),s); hold on;
    plot(ts-1,temp.mu{s},'color',colr(s,:));
    xlabel('time since target offset [s]'); ylabel('Target-tracking error [deg]'); axis([0 8 0 30]);
end

%% Steering error vs TT error - 2 taus / peakvel / % of distance cutoff
% Correlation of Behavioural error and Target-Tracking error over various % of Distance
DistPerc = linspace(0,1,101); % DistPerc = [ .4 .5 .6 .7 .75 .8 .85 .9 .95 .99]; DistPerc = [ .7 .75 .8 .85 .9 .95 .99];
params = 'stimtype'; BELIEF = 0;
[rho,pval,rho_sep] = BehvEvsTTECorr_DistancePerc(subject,params,DistPerc,BELIEF);
plevel = 0.05;
errtype = 'inst';
[Rsignif,peakind] = PeakCorrSignificance(rho_sep.(errtype).all.mu,pval.(errtype).all,plevel); nansum(Rsignif)

[rho,pval,rho_sep] = BehvEvsTTECorr_Time(subject,params,BELIEF); % minimum 5 sac/trial

%% Final Saccade Analysis
SaccadeAnalysis;

%% TT error timeseries
%% Modalities
% subjects separately
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
upperbnd = 0;
BELIEF = 0;
% [~,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
if strcmp(params,'stimtau'); colr = repmat(brewermap(size(poolindx,2)/3,'Dark2'),3,1); else; colr = brewermap(size(poolindx,2),'Dark2'); end
for i = 1:length(subject)
    figure; hold on;
    count = 1;
    for s = 1:size(poolindx,2)
        eye_err = eye_movement{i,s}.eyepos.behvcorr.eye_err_ts.mu;
        eye_std = eye_movement{i,s}.eyepos.behvcorr.eye_err_ts.std;
        ts = eye_movement{i,s}.eyepos.behvcorr.ts;
        nanindx = ~isnan(eye_err);
        eye_err = eye_err(nanindx); eye_std = eye_std(nanindx); ts = ts(nanindx);

        if strcmp(params,'stimtau'); subplot(1,3,count); hold on; end
        shadedErrorBar(ts,eye_err,eye_std,'lineprops',{'color',colr(s,:)});
        xlabel('time [s]'); ylabel('Target-Tracking error'); xlim([0 10]); ylim([0 60]);
        
        if mod(s,size(poolindx,2)/3) == 0; count = count + 1; legend(legend_input{s-size(poolindx,2)/3 +1:s}); end
    end
    sgtitle(subject(i).name);
end

%% subject average - Time
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
upperbnd = 0;
BELIEF = 0;
% [eye_fixation,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
if strcmp(params,'stimtau'); colr = repmat(brewermap(size(poolindx,2)/3,'Dark2'),3,1); else; colr = brewermap(size(poolindx,2),'Dark2'); end
AVGeye =[]; SEMeye = []; AVGts = [];
figure; hold on;
EyeLength = [];
for s = 1:size(poolindx,2)
    for i = 1:length(subject)
        EyeLength(i,s) = length(eye_movement{i,s}.eyepos.behvcorr.ts);
    end
    Nt = max(EyeLength(:));
end
count = 1;
for s = 1:size(poolindx,2)
    eye_err = []; ts = [];
    for i = 1:length(subject)
        eye_err(i,:) = [eye_movement{i,s}.eyepos.behvcorr.eye_err_ts.mu ; nan(Nt - length(eye_movement{i,s}.eyepos.behvcorr.ts),1)];
        ts(i,:) = [eye_movement{i,s}.eyepos.behvcorr.ts ; nan(Nt - length(eye_movement{i,s}.eyepos.behvcorr.ts),1)];
    end
    AVGeye(s,:) = nanmean(eye_err,1);
    SEMeye(s,:) = nanstd(eye_err,[],1)/sqrt(length(subject));
    AVGts(s,:) = nanmean(ts,1);
    
    if strcmp(params,'stimtau'); subplot(1,3,count); hold on; end
        shadedErrorBar(AVGts(s,:),AVGeye(s,:)-AVGeye(s,1),SEMeye(s,:),'lineprops',{'color',colr(s,:)});
        xlabel('time [s]'); ylabel('Target-Tracking error'); xlim([0 10]); ylim([0 60]);
        
        if mod(s,size(poolindx,2)/3) == 0; count = count + 1; legend(legend_input{s-size(poolindx,2)/3 +1:s}); end
end
if size(poolindx,2) <= 3; legend(legend_input{:}); end
sgtitle(['All subjects, BELIEF = ' num2str(BELIEF)]);
%% Cross-correlation of observed and predicted eye position
params = 'stimtype';
subject = subject_backup;
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
default_prs; upperbnd = 0; BELIEF = 0; 
prs.boots = 0; 
[poolindx,legend_input] = get_poolindx(subject,params);
if 0
[~,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
end

% subjects separately
figure;
for s = 1:size(poolindx,2) 
    for i = 1:length(subject)
       % vertical
       lag = dt*eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.lag;
       crosscorr = eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.c;
       crosscorr_shuffled = eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.c_shuffled;
       
       subplot(2,size(poolindx,2),s); hold on;
       plot(lag,crosscorr,'color',colr(s,:));
       plot(lag,crosscorr_shuffled,'--','color',colr(s,:)); title('vertical');
       axis([-5 5 -0.5 1]); vline(0,'k');
       xlabel('lag (obs-pred) [s]'); ylabel('cross-correlation (obs,pred)')

       % horizontal
       lag = dt*eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.lag;
       crosscorr = eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.c;
       crosscorr_shuffled = eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.c_shuffled;
       
       subplot(2,size(poolindx,2),s+size(poolindx,2)); hold on;
       plot(lag,crosscorr,'color',colr(s,:));
       plot(lag,crosscorr_shuffled,'--','color',colr(s,:)); title('horizontal');
       axis([-5 5 -0.5 1]); vline(0,'k');
       xlabel('lag (obs-pred) [s]'); ylabel('cross-correlation (obs,pred)')
    end
end
sgtitle(['BELIEF = ' num2str(BELIEF)])

% subject average
figure;
crosscorr = [];
for s = 1:size(poolindx,2) 
    for i = 1:length(subject)
       % vertical
       lag = dt*eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.lag;
       crosscorr.true.ver(i,:) = eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.c;
       crosscorr.shuffled.ver(i,:) = eye_movement{i,s}.eyepos.pred_vs_true.ver_mean.xcorr.c_shuffled;
       % horizontal
       lag = dt*eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.lag;
       crosscorr.true.hor(i,:) = eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.c;
       crosscorr.shuffled.hor(i,:) = eye_movement{i,s}.eyepos.pred_vs_true.hor_mean.xcorr.c_shuffled;
    end
       subplot(2,size(poolindx,2),s); hold on;
       shadedErrorBar(lag,nanmean(crosscorr.true.ver),nanstd(crosscorr.true.ver)./sqrt(length(subject)),'lineprops',{'color',colr(s,:)});
       plot(lag,nanmean(crosscorr.shuffled.ver),'--','color',colr(s,:)); title('vertical');
       axis([-5 5 -0.5 1]); vline(0,'k');
       xlabel('lag (obs-pred) [s]'); ylabel('cross-correlation (obs,pred)'); legend({'data','shuffled'},'location','southeast');

       subplot(2,size(poolindx,2),s+size(poolindx,2)); hold on;
       shadedErrorBar(lag,nanmean(crosscorr.true.hor),nanstd(crosscorr.true.hor)./sqrt(length(subject)),'lineprops',{'color',colr(s,:)});
       plot(lag,nanmean(crosscorr.shuffled.hor),'--','color',colr(s,:)); title('horizontal');
       axis([-5 5 -0.5 1]); vline(0,'k');
       xlabel('lag (obs-pred) [s]'); ylabel('cross-correlation (obs,pred)')
end
sgtitle(['BELIEF = ' num2str(BELIEF)])

%%
%%
%% LATEST APPROACH
%%
%%
%% Extract time-based eye movement analysis
if strcmp(subject_backup(5).name,'Ding'); subject_backup(5) = []; end
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
rmsac = 0;

fldcont = dir(data_folder);

% sensory condition groups
struct_indx = arrayfun(@(x) contains(x.name,'tracking_regular.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_regular')
        load(fldcont(struct_indx).name);
        disp('.....tracking_regular Loaded.')
    else; disp('tracking_regular already exists.')
    end
else
    params = 'stimtype';
    tracking_regular = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular.mat','tracking_regular','-v7.3')
end

% all trials
struct_indx = arrayfun(@(x) contains(x.name,'tracking_regular_all.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_regular_all')
        load(fldcont(struct_indx).name);
        disp('.....tracking_regular_all Loaded.')
    else; disp('tracking_regular_all already exists.')
    end
else    
    params = [];
    tracking_regular_all = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular_all.mat','tracking_regular_all','-v7.3')
end

Nsubs = size(tracking_regular,1);
Nstim = size(tracking_regular,2);
dt = 1/60;
xaxislim = [0 10];  xaxislabel = 'time [s]';
%% Re-sample quantities based on percentage of distance traveled
default_prs; 
align_flag = []; % 'align2targ';
dd_distperc = 0.005;

struct_indx = arrayfun(@(x) contains(x.name,'tracking_distperc.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_distperc')
        load(fldcont(struct_indx).name);
        disp('.....tracking_distperc Loaded.')
    else; disp('tracking_distperc already exists.')
    end
else
    tracking_distperc = EyeAnalysisOverDistancePercentage(tracking_regular,dd_distperc,align_flag,prs);
    save('tracking_distperc.mat','tracking_distperc','-v7.3')
end

xaxislim = [0 1];  xaxislabel = 'distance %';
%% Re-sample quantities based on distance to the end (end aligned)
default_prs; 
align_flag = [];
dd_dist2end = 5;

struct_indx = arrayfun(@(x) contains(x.name,'tracking_dist2end.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_dist2end')
        load(fldcont(struct_indx).name);
        disp('.....tracking_dist2end Loaded.')
    else; disp('tracking_dist2end already exists.')
    end
else
    tracking_dist2end = EyeAnalysisOverDistancetoEnd(tracking_regular,dd_dist2end,align_flag,prs);
    save('tracking_dist2end.mat','tracking_dist2end','-v7.3')
end

xaxislim = [-600 0];  xaxislabel = 'distance to end [cm]';

%% Eye position as a function of time, distance %, or absolute distance

% I = 1;
Plot_EyePositions(tracking_regular(I,:),subject_name(I));

%% Individual trials examples
condition = {'vestibular','visual','combined'};
jsthresh = 0.1;

% check whether button was recorded
subject_backup = buttonRecFlag(subject_backup);


% Individual trials
a = figure(345); a.Position = [12 327 320 575];
for i = 1:Nsubs
    for s = 2:Nstim
        trlindx = poolindx{i,s};
        for k = 1:numel(trlindx)
            ts = subject(i).trials(trlindx(k)).continuous.ts;
            tau = subject(i).trials(trlindx(k)).prs.tau;
            
            try     t_push = subject(i).trials(trlindx(k)).events.push_on(1);
                    t_end = subject(i).trials(trlindx(k)).events.t_end;
                    t_beg = subject(i).trials(trlindx(k)).events.t_beg;
                    
                    % steering end (joystick let-go)
                    jsy = -subject(i).trials(trlindx(k)).mc.JS_X_Raw;           jsx = subject(i).trials(trlindx(k)).mc.JS_Yaw_Raw;
                    vmax = subject(i).trials(trlindx(k)).prs.vmax;            wmax = subject(i).trials(trlindx(k)).prs.wmax;
                    try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                        t_letgo = subject(i).trials(trlindx(k)).continuous.ts(letgo_i);
                    catch;  t_letgo = nan; end
                    
                    % get eye position
                    ye = 0.5*(subject(i).trials(trlindx(k)).continuous.yle + subject(i).trials(trlindx(k)).continuous.yre);
                    ze = 0.5*(subject(i).trials(trlindx(k)).continuous.zle + subject(i).trials(trlindx(k)).continuous.zre);
                    
                    % get target position
                    phi = subject(i).trials(trlindx(k)).continuous.phi;
                    x_monk = subject(i).trials(trlindx(k)).continuous.xmp;	y_monk = subject(i).trials(trlindx(k)).continuous.ymp;
                    fireflyposx = subject(i).trials(trlindx(k)).prs.fireflyposx;   fireflyposy = subject(i).trials(trlindx(k)).prs.fireflyposy;
                    x_fly = fireflyposx - x_monk;     y_fly = fireflyposy - y_monk;
                    R = @(phi) [cosd(phi)   sind(phi) ; -sind(phi)  cosd(phi)];
                    XY_fly = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi, x_fly, y_fly,'uniformoutput',false));
                    x_fly = XY_fly(:,1);    y_fly = XY_fly(:,2);
                    
                    y = y_fly; y(y < 0) = nan;
                    x = x_fly;
                    z = -prs.height;
                    delta = prs.interoculardist;
                    [yle_targ, zle_targ, yre_targ, zre_targ] = world2eye(x,y,z,delta);
                    yt = 0.5*(yle_targ + yre_targ);	zt = 0.5*(zle_targ + zre_targ);
                    
                    % plot
                    clf;
                    subplot(2,1,1); hold on;
                    plot(ts,ye,'k'); plot(ts,yt,'--','color',[.5 .5 .5]); vline(t_push,'r--'); 
                    xlabel('time [s]'); ylabel('HOR [deg]'); title([condition{s} ': Trial = ' num2str(trlindx(k))]);
                    
                    subplot(2,1,2); hold on;
                    plot(ts,ze,'k'); plot(ts,zt,'--','color',[.5 .5 .5]);vline(t_push,'r--'); 
                    xlabel('time [s]'); ylabel('VER [deg]');
                    suptitle(subject(i).name); 
            catch
            end
             
        end
    end
end

%% Check fluctuation of initial tracking errors across trials
% check whether drift correction is needed
wn = 3;
for i = 1:Nsubs
    figure('name',[subject(i).name ': initial error fluctuation across trials'],'numbertitle','off','position',[0 700 950 275]);
    for s = 1:Nstim
        err_hor = tracking_regular{i,s}.eyepos.error.targ.hor(:,1);
        err_ver = tracking_regular{i,s}.eyepos.error.targ.ver(:,1);
        
        subplot(1,Nstim,s); hold on;
        plot(movmean(err_ver,wn),'color',colr(s,:)); 
        plot(movmean(err_hor,wn),'color',colr(s,:)*0.5); xlabel('No. trials'); ylabel('initial TTE [deg]'); legend('vertical','horizontal')
        axis([0 300 -30 30 ]);
        vline(0:30:250,'k:'); hline(0,'k');
    end
end

%% Number of saccades around end of steering
dt = 1/60; 
figure;
for i = keepindx
    temp_indx = arrayfun(@(x) find(abs(x.mc.JS_X_Raw) > 1,1,'last'),subject(i).trials,'un',0);
    temp_indx1 = temp_indx;
    nanindx = cellfun(@(x) isempty(x), temp_indx);
    temp_indx1(nanindx) = {nan};
    JSend_time{i} = cell2mat(temp_indx1)*dt;
   
    ntrls = numel(JSend_time{i});
    trlength = arrayfun(@(x) x.continuous.ts(end), subject(i).trials);
    rand_time =  rand(1,ntrls).*(trlength-1)+1;
    
    t_sac_temp = arrayfun(@(x) x.events.t_sac(x.events.t_sac - x.events.t_beg > 1) - x.events.t_beg, subject(i).trials,'un',0);
    t_sac_rel{i} = arrayfun(@(n) (t_sac_temp{n} - JSend_time{i}(n))', 1:ntrls,'un',0);
    t_sac_rel{i} = cell2mat(t_sac_rel{i});
    t_sac_rel_chance{i} = arrayfun(@(n) (t_sac_temp{n} - rand_time(n))', 1:ntrls,'un',0);
    t_sac_rel_chance{i} = cell2mat(t_sac_rel_chance{i});

    subplot(2,1,1);histogram(t_sac_rel{i},'normalization','probability'); vline(0,'k'); hold on; 
    axis([-20 20 0 0.3]); title('saccade probability around end of steering'); 
    subplot(2,1,2);histogram(t_sac_rel_chance{i},'normalization','probability'); vline(0,'k'); hold on; 
    axis([-20 20 0 0.3]); title('chance');
end

t_rel = cellfun(@(x) nanmean(x), t_sac_rel);
t_rel_chance = cellfun(@(x) nanmean(x), t_sac_rel_chance);

%% End of steering histogram
dt = 1/60; 
figure;
for i = keepindx
    temp_indx = arrayfun(@(x) find(abs(x.mc.JS_X_Raw) > 1,1,'last'),subject(i).trials,'un',0);
    temp_indx1 = temp_indx;
    nanindx = cellfun(@(x) isempty(x), temp_indx);
    temp_indx1(nanindx) = {nan};
    JSend_time{i} = cell2mat(temp_indx1)*dt;
    
    dist = arrayfun(@(x) cumsum(x.continuous.v)./max(cumsum(x.continuous.v)), subject(i).trials,'un',0);
    dist_indx = arrayfun(@(n) dist{n}(temp_indx{n}), 1:numel(dist),'un',0);
    nanindx = cellfun(@(x) isempty(x), dist_indx);
    dist_indx(nanindx) = {nan};
    JSend_dist{i} = cell2mat(dist_indx);

    
    subplot(2,1,1);
    histogram(JSend_time{i},'normalization','probability'); ylim([0 0.7]); hold on; xlabel('time');
    subplot(2,1,2);
    histogram(JSend_dist{i},'normalization','probability'); ylim([0 0.7]); hold on; xlabel('distperc');
end
title('Probability of ending steering')
%% Show differences in across trial variability across different plotting modes

[var_ver_regular,var_hor_regular] = EyePosVar(tracking_regular(keepindx,:));
[var_ver_distperc,var_hor_distperc] = EyePosVar(tracking_distperc(keepindx,:));

Nt_regular = size(var_hor_regular.mu,1);
Nt_distperc = size(var_hor_distperc.mu,1);
figure;
for s = 1:size(var_hor_regular.mu,2)
subplot(2,1,1);hold on;shadedErrorBar(1:Nt_regular,var_hor_regular.mu(:,s),var_hor_regular.se(:,s),'lineprops',{'color','k'}); ylim([0 3]); title('horizontal time');
subplot(2,1,1);hold on;shadedErrorBar((1:Nt_distperc)*5,var_hor_distperc.mu(:,s),var_hor_distperc.se(:,s),'lineprops',{'color','r'});ylim([0 3]); title('horizontal distperc');

subplot(2,1,2);hold on;shadedErrorBar(1:Nt_regular,var_ver_regular.mu(:,s),var_ver_regular.se(:,s),'lineprops',{'color','k'});ylim([0 3]); title('vertical time');
subplot(2,1,2);hold on;shadedErrorBar((1:Nt_distperc)*5,var_ver_distperc.mu(:,s),var_ver_distperc.se(:,s),'lineprops',{'color','r'});ylim([0 3]); title('vertical distperc');
end
suptitle('across trial variability')
%% Show differences in across trial mean eye positions across different plotting modes

[mu_ver_regular,mu_hor_regular] = EyePosMean(tracking_regular);
[mu_ver_distperc,mu_hor_distperc] = EyePosMean(tracking_distperc);

Nt_regular = size(mu_hor_regular.mu,1);
Nt_distperc = size(mu_hor_distperc.mu,1);

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure;
for s = 1:size(mu_hor_regular.mu,2)
    subplot(2,1,1);hold on;
    shadedErrorBar((1:Nt_distperc)/2,mu_hor_distperc.mu(:,s),mu_hor_distperc.se(:,s),'lineprops',{'color',colr(s,:)});
    title('horizontal'); xlabel('distance %'); xticks([0 25 50 75 100]); yticks([0 0.5 1]);
    
    subplot(2,1,2);hold on;
    shadedErrorBar((1:Nt_distperc)/2,mu_ver_distperc.mu(:,s),mu_ver_distperc.se(:,s),'lineprops',{'color',colr(s,:)});
    title('vertical'); xlabel('distance %'); xticks([0 25 50 75 100]); yticks([-2 -1 0]);
end


% for s = 1:size(mu_hor_regular.mu,2)
% subplot(2,1,1);hold on;
% shadedErrorBar(1:Nt_regular,mu_hor_regular.mu(:,s),mu_hor_regular.se(:,s),'lineprops',{'color','k'}); title('horizontal');
% shadedErrorBar((1:Nt_distperc)*5,mu_hor_distperc.mu(:,s),mu_hor_distperc.se(:,s),'lineprops',{'color','r'}); 
% 
% subplot(2,1,2);hold on;shadedErrorBar(1:Nt_regular,mu_ver_regular.mu(:,s),mu_ver_regular.se(:,s),'lineprops',{'color','k'}); title('vertical');
% shadedErrorBar((1:Nt_distperc)*5,mu_ver_distperc.mu(:,s),mu_ver_distperc.se(:,s),'lineprops',{'color','r'}); 
% end
% legend({'time','distperc'});
suptitle('normalized mean eye positions')
%% Regression of actual vs predicted eye positions when target OFF
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
% keepindx = 1:numel(subject_backup);
tracking = tracking_regular_all(keepindx);
Nsubs = numel(tracking);

bV = []; bH = [];
for i = 1:Nsubs
    hortar = cell2mat(cellfun(@(x) x.tarpos.screen.hor_mean(:,1)', tracking(i,:),'un',0));
    vertar = cell2mat(cellfun(@(x) x.tarpos.screen.ver_mean(:,1)', tracking(i,:),'un',0));
    horeye = cell2mat(cellfun(@(x) x.eyepos.screen.hor_mean(:,1)', tracking(i,:),'un',0));
    vereye = cell2mat(cellfun(@(x) x.eyepos.screen.ver_mean(:,1)', tracking(i,:),'un',0));
    bV{i} = regress(vereye(:),[vertar(:) ones(numel(vertar),1)]);
    bH{i} = regress(horeye(:),[hortar(:) ones(numel(hortar),1)]);
    
    if 1
    figure; 
    subplot(1,2,1); hold on;
    plot(hortar,horeye,'.','color',[.5 .5 .5]); plot(-40:40,-40:40,'k--'); plot(-40:40,bH{i}(1)*(-40:40),'r','linewidth',2);
    vline(0,'k'); hline(0,'k'); axis equal; axis([-40 40 -40 40]); xlabel('Predicted [deg]'); ylabel('Actual [deg]'); title('Horizontal');
    
    subplot(1,2,2); hold on;
    plot(vertar,vereye,'.','color',[.5 .5 .5]); plot(-30:0,-30:0,'k--'); plot(-30:0,bV{i}(1)*(-30:0),'r','linewidth',2); 
    vline(0,'k'); hline(0,'k'); axis equal; axis([-30 10 -30 10]); xlabel('Predicted [deg]'); ylabel('Actual [deg]'); title('Vertical');
    end
end
%% Target tracking error when target OFF
subject = subject_backup;
default_prs;
rmsac = 0;
if 0
tracking_regular_all = EyeAnalysisOverTime(subject,[],rmsac,prs);
end
outlier_thresh = 25;
keepindx = [1 2 7 9 10 12 13 14];
subindx = keepindx;
% subindx = 1:numel(subject);
Nsubs = numel(subindx);


Nboots = 100;
all_err = [];
for i = 1:Nsubs
    indx = cell2mat(cellfun(@(x) abs(x.eyepos.error.targ.all(:,1)) < outlier_thresh, tracking_regular_all(subindx(i)),'un',0));
    errors = cell2mat(cellfun(@(x) x.eyepos.error.targ.all(indx,1)', tracking_regular_all(subindx(i)),'un',0));
    
    ntrls = length(errors);
    for n = 1:Nboots
    randindx = randsample(ntrls,ntrls,1);
    all_err(n,i) = nanmean(errors(randindx));
    end
end
allerr.se = arrayfun(@(i) nanstd(all_err(:,i)), 1:Nsubs);
allerr.mu = arrayfun(@(i) nanmean(all_err(:,i)), 1:Nsubs);

% chance level
tracking_all = tracking_regular_all(keepindx,:);
eye_hor_shuffled = cellfun(@(x) x.eyepos.screen.hor_mean(randi(size(x.eyepos.screen.hor_mean,1),size(x.eyepos.screen.hor_mean,1),1),1),tracking_all,'un',0);
eye_ver_shuffled = cellfun(@(x) x.eyepos.screen.ver_mean(randi(size(x.eyepos.screen.ver_mean,1),size(x.eyepos.screen.ver_mean,1),1),1),tracking_all,'un',0);
tar_hor_shuffled = cellfun(@(x) x.tarpos.screen.hor_mean(randi(size(x.tarpos.screen.hor_mean,1),size(x.tarpos.screen.hor_mean,1),1),1),tracking_all,'un',0);
tar_ver_shuffled = cellfun(@(x) x.tarpos.screen.ver_mean(randi(size(x.tarpos.screen.ver_mean,1),size(x.tarpos.screen.ver_mean,1),1),1),tracking_all,'un',0);
err_all_time_shuffled = cell2mat(cellfun(@(x1,x2,y1,y2) nanmean(sqrt((x1-x2).^2 + (y1-y2).^2)),eye_hor_shuffled,tar_hor_shuffled,eye_ver_shuffled,tar_ver_shuffled,'un',0));
% mean error
tracking_all = tracking_distperc(keepindx,:);
errors_mu_ves = nanmean((cellfun(@(x) nanmean(nanmean(x.eyepos.error.targ.all(:,3:end),2)), tracking_all(:,1))));
errors_mu_vis = nanmean((cellfun(@(x) nanmean(nanmean(x.eyepos.error.targ.all(:,3:end),2)), tracking_all(:,2))));

% figure; errorbar(1:Nsubs,allerr.mu,allerr.se,'s','color','k','markerfacecolor','k','linestyle','none','capsize',0); axis([0 Nsubs+1 0 30]);
% xlabel('Subjects'); ylabel('Target-tracking error [deg]'); title('Target OFF');

figure; hold on;
bar(1,nanmean(allerr.mu),'facecolor',[.5 .5 .5]); hline(nanmean(err_all_time_shuffled(:,1)),'k-');
plot(1+0.1*randn(numel(allerr.mu),1),allerr.mu,'ko','markerfacecolor','k'); 
hline(errors_mu_ves,'r'); hline(errors_mu_vis,'b'); ylim([0 30]);
ylabel('Target-tracking error [deg]'); title('Target OFF');

if 0
    keepindx1 = setdiff(1:numel(subject),find(allerr.mu > 10));
end
%% Average target tracking error over time and distance %
all_subs = 0;
% keepindx = [1 3 4 5 8 9 10 11 13];
keepindx = [1 2 7 9 10 12 13 14];
tracking = tracking_regular; % tracking_regular(keepindx,:);
Nsubs = size(tracking,1);
Nstim = size(tracking,2);
minlength = min(cellfun(@(x) size(x.eyepos.error.targ.all,2),tracking),[],'all');
ts = (1:minlength)/60;
for s = 1:Nstim
    err_all_time{s} = cell2mat(cellfun(@(x) nanmean(x.eyepos.error.targ.all(:,1:minlength)),tracking(:,s),'un',0));
end
% chance level
rng(0);
tracking_all = tracking_regular_all;%(keepindx,:);
eye_hor_shuffled = cellfun(@(x) x.eyepos.screen.hor_mean(randi(size(x.eyepos.screen.hor_mean,1),size(x.eyepos.screen.hor_mean,1),1),1:minlength),tracking_all,'un',0);
eye_ver_shuffled = cellfun(@(x) x.eyepos.screen.ver_mean(randi(size(x.eyepos.screen.ver_mean,1),size(x.eyepos.screen.ver_mean,1),1),1:minlength),tracking_all,'un',0);
tar_hor_shuffled = cellfun(@(x) x.eyepos.screen.hor_mean(randi(size(x.tarpos.screen.hor_mean,1),size(x.tarpos.screen.hor_mean,1),1),1:minlength),tracking_all,'un',0);
tar_ver_shuffled = cellfun(@(x) x.eyepos.screen.ver_mean(randi(size(x.tarpos.screen.ver_mean,1),size(x.tarpos.screen.ver_mean,1),1),1:minlength),tracking_all,'un',0);

err_all_time_shuffled = cell2mat(cellfun(@(x1,x2,y1,y2) nanmean(sqrt((x1-x2).^2 + (y1-y2).^2)),eye_hor_shuffled,tar_hor_shuffled,eye_ver_shuffled,tar_ver_shuffled,'un',0));

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure; subplot(1,2,1); hold on;
if all_subs
    arrayfun(@(s) plot(ts,[err_all_time{s}],'color',colr(s,:)),1:Nstim);
else
    arrayfun(@(s) shadedErrorBar(ts,nanmean(err_all_time{s}),nanstd(err_all_time{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}), 1:Nstim);
end
plot(ts,nanmean(err_all_time_shuffled),'-k');
xlabel('time [s]'); ylabel('target tracking error [deg]'); title('Time');
axis([0 12 0 40]); legend({'vestibular','visual','combined'},'location','southeast');

% normalized plot
subplot(1,2,2); hold on;
if all_subs
    arrayfun(@(s) plot(ts,[err_all_time{s}./nanmean(err_all_time_shuffled)],'color',colr(s,:)),1:Nstim)
else
    arrayfun(@(s) shadedErrorBar(ts,nanmean(err_all_time{s}./nanmean(err_all_time_shuffled)),nanstd(err_all_time{s})./sqrt(Nsubs)./nanmean(err_all_time_shuffled),'lineprops',{'color',colr(s,:)}), 1:Nstim);
end
hline(1,'k-');
xlabel('time [s]'); ylabel('normalized TTE'); title('Time');
axis([0 12 0 2]); legend({'vestibular','visual','combined'},'location','southeast');

% distance %
tracking = tracking_distperc(keepindx,:);
Nsubs = size(tracking,1);
Nstim = size(tracking,2);
minlength = min(cellfun(@(x) size(x.eyepos.error.targ.all,2),tracking),[],'all');
ts = (1:minlength)/minlength*100;
err_all_distperc = [];
for s = 1:Nstim
    err_all_distperc{s} = cell2mat(cellfun(@(x) nanmean(x.eyepos.error.targ.all(:,1:minlength)),tracking(:,s),'un',0));
end
% chance level
eye_hor_shuffled = cellfun(@(x) x.eyepos.screen.hor_mean(randi(size(x.eyepos.screen.hor_mean,1),size(x.eyepos.screen.hor_mean,1),1),1:minlength),tracking(:,1),'un',0);
eye_ver_shuffled = cellfun(@(x) x.eyepos.screen.ver_mean(randi(size(x.eyepos.screen.ver_mean,1),size(x.eyepos.screen.ver_mean,1),1),1:minlength),tracking(:,1),'un',0);
tar_hor_shuffled = cellfun(@(x) x.eyepos.screen.hor_mean(randi(size(x.tarpos.screen.hor_mean,1),size(x.tarpos.screen.hor_mean,1),1),1:minlength),tracking(:,1),'un',0);
tar_ver_shuffled = cellfun(@(x) x.eyepos.screen.ver_mean(randi(size(x.tarpos.screen.ver_mean,1),size(x.tarpos.screen.ver_mean,1),1),1:minlength),tracking(:,1),'un',0);

err_all_time_shuffled = cell2mat(cellfun(@(x1,x2,y1,y2) nanmean(sqrt((x1-x2).^2 + (y1-y2).^2)),eye_hor_shuffled,tar_hor_shuffled,eye_ver_shuffled,tar_ver_shuffled,'un',0));

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure;subplot(1,2,1); hold on;
if all_subs
    arrayfun(@(s) plot(ts,[err_all_distperc{s}],'color',colr(s,:)),1:Nstim);
else
    arrayfun(@(s) shadedErrorBar(ts,nanmean(err_all_distperc{s}),nanstd(err_all_distperc{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)}), 1:Nstim);
end
plot(ts,nanmean(err_all_time_shuffled),'-k');
xlabel('distance percentage (%)'); ylabel('target tracking error [deg]'); title('Distance percentage');
axis([0 100 0 40]); legend({'vestibular','visual','combined'},'location','southeast');

% normalized plot
subplot(1,2,2); hold on;
if all_subs
    arrayfun(@(s) plot(ts,[err_all_distperc{s}./nanmean(err_all_time_shuffled)],'color',colr(s,:)),1:Nstim)
else
    arrayfun(@(s) shadedErrorBar(ts,nanmean(err_all_distperc{s}./nanmean(err_all_time_shuffled)),nanstd(err_all_distperc{s})./sqrt(Nsubs)./nanmean(err_all_time_shuffled),'lineprops',{'color',colr(s,:)}), 1:Nstim);
end
hline(1,'k-');
xlabel('distance percentage (%)'); ylabel('normalized TTE'); title('Distance percentage');
axis([0 100 0 2]); legend({'vestibular','visual','combined'},'location','southeast');

% quantify chancel-level crossing
norm_err = cellfun(@(x) x./nanmean(err_all_time_shuffled),err_all_distperc,'un',0);
for s = 1:size(norm_err,2)
   perc_indx{s} = arrayfun(@(i) find(norm_err{s}(i,:) <= 1,1,'last'), 1:size(norm_err{s},1)) ;
end

perc_val = cell2mat(cellfun(@(x) ts(x)',perc_indx,'un',0));
disp('Chance level crossing mean +/- SD:')
mean(perc_val)
std(perc_val)./sqrt(Nsubs)
%% Tracking error timeseries, normalized by average TTE when target OFF 
% or prediction based on regression of actual vs predicted eye position






%% Tracking errors over time/distance % (NORMALIZED)
errtype = 'targ';
[err_ver_regular,err_hor_regular,err_all_regular] = AverageTrackingError(tracking_regular,errtype);
[err_ver_distperc,err_hor_distperc,err_all_distperc] = AverageTrackingError(tracking_distperc,errtype);

Nt_regular = size(err_ver_regular.mu,1);
Nt_distperc = size(err_ver_distperc.mu,1);

figure; 
for s = 1:size(err_ver_regular.mu,2)
subplot(1,2,1);hold on;shadedErrorBar(1:Nt_regular,err_all_regular.mu(:,s),err_all_regular.se(:,s),'lineprops',{'color','k'});  title('error time'); ylim([0 70]);
subplot(1,2,2);hold on;shadedErrorBar(1:Nt_distperc,err_all_distperc.mu(:,s),err_all_distperc.se(:,s),'lineprops',{'color','r'}); title('error distperc'); ylim([0 70]);

% subplot(1,2,1);hold on;plot(cumsum(diff(err_all_regular.mean(:,s))),'k');  title('error time'); %ylim([-1 2]);
% subplot(1,2,2);hold on;plot(cumsum(diff(err_all_distperc.mean(:,s))),'r'); title('error distperc'); ylim([-1 2]);
end
%% Check average tracking errors diference (target vs stop position)

% over time
TrackingErrorsDifference(tracking_regular);

% over distance percentage
TrackingErrorsDifference(tracking_distperc);

% over distance to end
TrackingErrorsDifference(tracking_dist2end);

%% Check vertical component
Nsubs = size(tracking_distperc,1);
for i = 1:Nsubs
    ver = figure('name',[subject(i).name ': vertical eye positions'],'numbertitle','off','position',[0 500 900 200]);
    for s = 1:Nstim
        % time
        ts = tracking_distperc{i,s}.misc.ts;
        % vertical
        figure(ver); 
        subplot(1,Nstim,s); hold on;
        plot(ts,tracking_distperc{i,s}.eyepos.screen.ver_mean','color',colr(s,:),'linewidth',0.1); xlabel('distance %');ylabel('vertical eye [deg]'); axis([0 1 -60 20]);
    end
end

%% Check differences in performance between vertical trackers and non-trackers

lateindx = 60; % percent
min_diff = 5;
% keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
subindx = FindVerticalTrackers(tracking_distperc(keepindx,:),lateindx,min_diff);

[bias,variab] = TrackersVsNonTrackersPerformance(subject(keepindx),subindx);

if 0 % hand-labeled trackers
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
temp = logical([1 0 0 0 0 1 0 0 1 0 1 0 0 1 ; 1 0 0 1 1 1 1 0 1 1 1 1 1 1 ; 1 0 0 1 1 1 1 0 1 1 1 1 1 1]');
Nsubs = numel(keepindx);
subindx = mat2cell(temp,Nsubs,ones(1,Nstim));

[bias,variab] = TrackersVsNonTrackersPerformance(subject(keepindx),subindx);
end

%% Ratio between SPTE and TTE (REMOVED FROM ANALYSIS)!!!!!
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = 1:numel(subject);
% keepindx = [1 3 4 5 8 9 10 11 13];
keepindx = [1 2 7 9 10 12 13 14];
% keepindx = 1:numel(subject);

stability_constant = 0;
[ratio_kernel,ts] = TTEvsSPTE_Ratio(tracking_distperc(keepindx,:),stability_constant);


% divide between vertical trackers and non-trackers
if 0
lateindx = 50;
min_diff = 5;
subindx = FindVerticalTrackers(tracking_distperc(keepindx,:),lateindx,min_diff);

figure('name','Split Trackers: Ratio between SPTE and TTE','numbertitle','off','position',[-10 690 970 290]);
bT = cellfun(@(x,ind) x(ind,:), ratio_kernel, subindx, 'uniformoutput',false);
bNT = cellfun(@(x,ind) x(~ind,:), ratio_kernel, subindx, 'uniformoutput',false);

kernel = bT;
for s = 1:Nstim
    subN = size(kernel{s},1);
    subplot(1,Nstim,s); hold on; plot(ts,kernel{s},'color',colr(s,:));
    shadedErrorBar(ts,mean(kernel{s}),std(kernel{s})./sqrt(subN),'lineprops',{'linewidth',3,'color',colr(s,:)}); axis([xaxislim 0 3]);
    hline(1,'k--'); xlabel(xaxislabel); ylabel('mean ratio kernel (TTE/SPTE)');
end

kernel = bNT;
for s = 1:Nstim
    subN = size(kernel{s},1);
    subplot(1,Nstim,s); hold on; plot(ts,kernel{s},'color',colr(s,:)*0.5);
    shadedErrorBar(ts,mean(kernel{s},1),std(kernel{s},[],1)./sqrt(subN),'lineprops',{'linewidth',3,'color',colr(s,:)*0.5}); axis([xaxislim 0 3]);
    hline(1,'k--'); xlabel(xaxislabel); ylabel('mean ratio kernel (TTE/SPTE)');
end
end
%% Correlation of average tracking error and steering error

tracking_all = tracking_distperc(keepindx,:);
tte_mu_ves = cellfun(@(x) nanmean(nanmean(x.eyepos.error.targ.all(:,1:end-2))'), tracking_all(:,1));
tte_mu_vis = cellfun(@(x) nanmean(nanmean(x.eyepos.error.targ.all(:,1:end-2))'), tracking_all(:,2));
be_mu_ves = cellfun(@(x) nanmean(x.misc.behverrors.all.val), tracking_all(:,1));
be_mu_vis = cellfun(@(x) nanmean(x.misc.behverrors.all.val), tracking_all(:,2));

[rho_ves,p_ves] = corr(be_mu_ves(:),tte_mu_ves(:));
[rho_vis,p_vis] = corr(be_mu_vis(:),tte_mu_vis(:));

figure; hold on; plot(be_mu_ves,tte_mu_ves,'ro','markerfacecolor','r');
plot(be_mu_vis,tte_mu_vis,'bo','markerfacecolor','b');
axis([0 150 0 30])
xlabel('average stopping error [cm]'); ylabel('average tracking-error [deg]');
legend({['\rho vest = ' num2str(rho_ves)],['\rho vis = ' num2str(rho_vis)]},'location','northwest')


% alternative (compute correlation per subject first)
tte_ves = cellfun(@(x) nanmean(x.eyepos.error.targ.all(:,1:end-2),2), tracking_all(:,1),'un',0);
be_ves = cellfun(@(x) x.misc.behverrors.all.val, tracking_all(:,1),'un',0);
[rho_ves,p_ves] = cellfun(@(x,y) corr(x(:),y(:)),be_ves,tte_ves);

tte_vis = cellfun(@(x) nanmean(x.eyepos.error.targ.all(:,1:end-2),2), tracking_all(:,2),'un',0);
be_vis = cellfun(@(x) x.misc.behverrors.all.val, tracking_all(:,2),'un',0);
[rho_vis,p_vis] = cellfun(@(x,y) corr(x(:),y(:)),be_vis,tte_vis);


%% Correlation between tracking error and Steering error

plevel = 0.05;
refplane = 'screen';
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];

% [rho,pval] = SteeringVsTrackingError(tracking_regular(keepindx,:),refplane);
[rho,pval] = SteeringVsTrackingError(tracking_distperc(keepindx,:),refplane);

[Rsignif,peakind] = PeakCorrSignificance1(rho,pval,plevel); nansum(Rsignif)

%% Position gain between eye and target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
tracking = tracking_distperc(keepindx,:);

[g,ts] = EyePositionGain(tracking);

%% Correlation between eye and target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];

[rho,ts] = EyePositionCorr(tracking_distperc(keepindx,:));
[rho,ts] = EyePositionCorr(tracking_regular(keepindx,:));

%% Multiple regression between eye & target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];

[multi,ts] = EyePositionMultiRegr(tracking_distperc(keepindx,:));

%% Generate target belief over time
dt=1/60;
tracking = tracking_distperc;
[multi,ts_distperc] = EyePositionMultiRegr(tracking);
close;

disp('.......Generating target position belief')
% figure;
for i = 1:size(tracking,1)
    disp(['Subject = ' num2str(i)])
    stimtype = arrayfun(@(x) x.prs.stimtype, subject(i).trials);
    
    for j = 1:numel(stimtype)

        v = subject(i).trials(j).continuous.v;
        w = subject(i).trials(j).continuous.w;
        xfp = subject(i).trials(j).prs.fireflyposx;
        yfp = subject(i).trials(j).prs.fireflyposy;
        x_start = subject(i).trials(j).continuous.xmp(1);
        y_start = subject(i).trials(j).continuous.ymp(1);
        
        [y_targ,z_targ,y_stop,z_stop] = egoTarg2Screen({v},{w},xfp,yfp,x_start,y_start,dt,prs);
        y_targ = y_targ{:};     z_targ = z_targ{:};
        y_stop = y_stop{:};     z_stop = z_stop{:};

        d_sub = cumsum(v)*dt;
        distperc = d_sub./max(d_sub);
        Nt = numel(v);
        
        [~,dp_indx] = arrayfun(@(t) min(abs(distperc(t)-ts_distperc)),1:Nt);
        stop_track_ind = find(dp_indx > 0.8*numel(ts_distperc),1);
        
        b = multi.hor.targ{stimtype(i)}(i,dp_indx)';        yb_targ = y_targ.*b;
        b = multi.ver.targ{stimtype(i)}(i,dp_indx)';        zb_targ = z_targ.*b;
        b = multi.hor.stop{stimtype(i)}(i,dp_indx)';        yb_stop = y_stop.*b;
        b = multi.ver.stop{stimtype(i)}(i,dp_indx)';        zb_stop = z_stop.*b;
        
        % generate belief
        ybp = yb_targ + yb_stop;    ybp(stop_track_ind:end) = nan;
        zbp = zb_targ + zb_stop;    zbp(stop_track_ind:end) = nan;
        
        if 0 
            clf;hold on;plot(ybp,zbp); plot(y_targ,z_targ); plot(y_stop,z_stop); 
            axis([-50 50 -90 10])
        end
        
        subject(i).trials(j).continuous.ybp = ybp;
        subject(i).trials(j).continuous.zbp = zbp;
    end
        
end
subject_backup = subject;

%% Check weight sum
tracking = tracking_distperc(keepindx,:); 
[multi,ts_distperc] = EyePositionMultiRegr(tracking);

weight_sum_hor = cellfun(@(t,s) t+s, multi.hor.targ, multi.hor.stop,'un',0);
weight_sum_ver = cellfun(@(t,s) t+s, multi.ver.targ, multi.ver.stop,'un',0);

figure;
for s = 1:Nstim-1
    subplot(2,2,s);hold on;
    plot(ts_distperc, weight_sum_hor{s},'color',colr(s,:))
    xlabel('distance %'); ylabel('weight sum'); title('HOR');hline([0 1],'k');axis([0 0.8 -0.5 1.5]);

    subplot(2,2,s+2);hold on;
    plot(ts_distperc, weight_sum_ver{s},'color',colr(s,:))
    xlabel('distance %'); ylabel('weight sum'); title('VER');hline([0 1],'k'); axis([0 0.8 -0.5 1.5]);
end

%%
%% New Saccade Analysis
%%
%%
if strcmp(subject_backup(5).name,'Ding'); subject_backup(5) = []; end
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = length(subject);
Nstim = size(poolindx,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

disp('Saccade Analysis..........')

fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'saccade_all.mat'), fldcont);

if any(struct_indx)
    if ~exist('saccade_all')
        load(fldcont(struct_indx).name);
        disp('.....saccade_all Loaded.')
    else; disp('tracking_distperc already exists.')
    end
else
    
    for i = 1:Nsubs
        for s = 1:Nstim
            indx = poolindx{i,s};
            saccade_all{i,s} = AnalyseSaccades(subject(i).trials(indx),prs);
        end
        disp(['........Subject = ' num2str(i)])
    end
    save('saccade_all.mat','saccade_all','-v7.3')
    disp('...... Saved saccade_all.mat');
end

%% Saccade Probability
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

SaccadeProbability(saccade);

%% Saccade amplitude for different trial epochs (Target ON, Steering, End of Trial)
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

SaccadeAmplitude_Epochs(saccade);

%% Tracking errors around saccade
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

Nsubs = size(saccade,1);
Nstim = size(saccade,2);

ngroups = 1;
edges = linspace(-0.1,1,ngroups+1);

targ_err = []; stop_err = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlindx = arrayfun(@(n) saccade{i,s}.distperc > edges(n) & saccade{i,s}.distperc < edges(n+1),1:ngroups,'un',0);
        for n = 1:ngroups
            err_hor_temp = saccade{i,s}.eyepos.saccadealigned.error.targ.hor(:,trlindx{n});
            err_ver_temp = saccade{i,s}.eyepos.saccadealigned.error.targ.ver(:,trlindx{n});
            targ_err.hor{s}(n,:,i) = nanmean(abs(err_hor_temp),2);
            targ_err.ver{s}(n,:,i) = nanmean(abs(err_ver_temp),2);
            targ_err.all{s}(n,:,i) = nanmean(sqrt(err_hor_temp.^2 + err_ver_temp.^2),2);
            
            err_hor_temp = saccade{i,s}.eyepos.saccadealigned.error.stop.hor(:,trlindx{n});
            err_ver_temp = saccade{i,s}.eyepos.saccadealigned.error.stop.ver(:,trlindx{n});
            stop_err.hor{s}(n,:,i) = nanmean(abs(err_hor_temp),2);
            stop_err.ver{s}(n,:,i) = nanmean(abs(err_ver_temp),2);
            stop_err.all{s}(n,:,i) = nanmean(sqrt(err_hor_temp.^2 + err_ver_temp.^2),2);
        end
    end
end
ts = saccade{1}.misc.ts;

figure; hold on;
for s = 1:Nstim
    for n = 1:ngroups
    % target position
    tte.mu = squeeze(nanmean(targ_err.all{s},3));
    tte.se = squeeze(nanstd(targ_err.all{s},[],3))./sqrt(Nsubs);
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts,tte.mu(n,:),tte.se(n,:),'lineprops',{'color',colr(s,:)/n});
    title('euclideian error'); ylabel('TTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
%     subplot(2,Nstim,2); hold on;
%     shadedErrorBar(ts,nanmean(targ_err.hor{s}),nanstd(targ_err.hor{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
%     title('horizontal error'); ylabel('TTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
%     subplot(2,Nstim,3); hold on;
%     shadedErrorBar(ts,nanmean(targ_err.ver{s}),nanstd(targ_err.ver{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
%     title('vertical error'); ylabel('TTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
    
    % stop position
    spte.mu = squeeze(nanmean(stop_err.all{s},3));
    spte.se = squeeze(nanstd(stop_err.all{s},[],3))./sqrt(Nsubs);
    subplot(2,Nstim,Nstim+s); hold on;
    shadedErrorBar(ts,spte.mu(n,:),spte.se(n,:),'lineprops',{'color',colr(s,:)/n});
    title('euclideian error'); ylabel('SPTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
%     subplot(2,Nstim,Nstim+2); hold on;
%     shadedErrorBar(ts,nanmean(stop_err.hor{s}),nanstd(stop_err.hor{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
%     title('horizontal error'); ylabel('SPTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
%     subplot(2,Nstim,Nstim+3); hold on;
%     shadedErrorBar(ts,nanmean(stop_err.ver{s}),nanstd(stop_err.ver{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
%     title('vertical error'); ylabel('SPTE'); xlabel('time since saccade onset'); ylim([5 40]); vline(0,'k');
    end
end

%% Tracking errors after saccade as a function of distance percentage
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[rho,pval,beta] = PostSaccadeErrorsVsDistancePerc(saccade);

%% Regress endpoint of saccade against target/stop positions before saccade and plot weights over distance %
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[rho,pval,beta,multi] = SaccadeEndpointRegression(saccade);

%% Regression KERNEL of saccade endpoint vs target/stop positions early and late in trial
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[rho,pval,beta,multi] = SaccadeEndpointRegressionKernel(saccade);
 
%% Regression KERNEL of saccade amplitude vs tracking errors (target, stop, belief) +++++++IN ANALYSIS+++++++
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);
type_regr = 'linear';
quantity = 'multi';
split_grp = 0;

[rho,pval,beta,multi] = SaccadeAmplitudeKernel(saccade,type_regr,quantity,split_grp);

% print peak stats
kernels = fieldnames(multi);
Nkernels = numel(kernels);
for k = 1:Nkernels
Nsubs = size(saccade,1);
disp(' '); disp(['++++++ Peak ' kernels{k} ' kernel mean ' char(177) ' SEM ++++++']); disp(' ')
kernel_mu = cellfun(@(x) nanmean(x,3), multi.(kernels{k}).hor,'un',0);
kernel_se = cellfun(@(x) nanstd(x,[],3)/sqrt(Nsubs), multi.(kernels{k}).hor,'un',0);
[~,peakind] = cellfun(@(x) max(nanmean(x,3)), multi.(kernels{k}).hor,'un',0);
peak_mu = cellfun(@(x,i) x(i), kernel_mu, peakind);
peak_se = cellfun(@(x,i) x(i), kernel_se, peakind);
disp('HORIZONTAL:'); 
disp(['vestibular: ' num2str(peak_mu(1)) char(177) num2str(peak_se(1))])
disp(['visual: ' num2str(peak_mu(2)) char(177) num2str(peak_se(2))])
disp(['combined: ' num2str(peak_mu(3)) char(177) num2str(peak_se(3))]); disp(' ');

kernel_mu = cellfun(@(x) nanmean(x,3), multi.(kernels{k}).ver,'un',0);
kernel_se = cellfun(@(x) nanstd(x,[],3)/sqrt(Nsubs), multi.(kernels{k}).ver,'un',0);
[~,peakind] = cellfun(@(x) max(nanmean(x,3)), multi.(kernels{k}).ver,'un',0);
peak_mu = cellfun(@(x,i) x(i), kernel_mu, peakind);
peak_se = cellfun(@(x,i) x(i), kernel_se, peakind);
disp('VERTICAL:'); 
disp(['vestibular: ' num2str(peak_mu(1)) char(177) num2str(peak_se(1))])
disp(['visual: ' num2str(peak_mu(2)) char(177) num2str(peak_se(2))])
disp(['combined: ' num2str(peak_mu(3)) char(177) num2str(peak_se(3))]); disp(' ');
end

%% Regression KERNEL of saccade amplitude vs tracking errors early and late in trial
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);
type_regr = 'linear';
quantity = 'multi';
split_grp = 1;

[rho,pval,beta,multi] = SaccadeAmplitudeKernel(saccade,type_regr,quantity,split_grp);

%% Paired differences of regression coefficients before saccade across subjects




%% Regress saccade amplitude against TTE and SPTE before saccade and plot weights over distance %
% play a bit here, seems promising
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[rho,pval,beta,multi] = SaccadeAmplitudeRegression(saccade);


%% Ratio of TTE/SPTE against distance %
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[rho,pval,beta] = PostSaccadeErrorsRatioVsDistancePerc(saccade);

%% Ratio of TTE/SPTE around saccade
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

[mu,sd] = ErrorsRatioKernel(saccade);

%% Sum of saccades vs steering errors correlation
keepindx = [1 2 7 9 10 12 13 14];
saccade = saccade_all(keepindx,:);

Nsubs = size(saccade,1);
Nstim = size(saccade,2);
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

Tmin = 1; Tmax = 5;
Dmin = 0; Dmax = 0.5;

angerr_thresh = 50;

r_corr = []; th_corr = []; r_err_pool = cell(1,Nstim); th_err_pool = cell(1,Nstim); ver_sac_amp_pool = cell(1,Nstim); hor_sac_amp_pool = cell(1,Nstim);
for i = 1:Nsubs
for s = 1:Nstim
    
    tsac = saccade{i,s}.time;
    dsac = saccade{i,s}.distperc;
    tsacindx = dsac > Dmin & dsac < Dmax;
    
    trlindx = saccade{i,s}.trlindx;
    valid_trlindx = trlindx(tsacindx);
    trls = unique(valid_trlindx);
    Ntrls = numel(trls);
    
    ver_sac_amp = []; hor_sac_amp = []; r_err = []; th_err = [];
    for n = 1:Ntrls
        
        % steering errors
        r_err(n) = saccade{i,s}.misc.behverrors.r.val(trls(n));
        th_err(n) = saccade{i,s}.misc.behverrors.th.val(trls(n));
        
        % get valid saccades
        ind = trlindx == trls(n);
        validindx = ind & tsacindx;
        ind = find(validindx);
        if ~isempty(ind)
        % vertical component
        ver_sac_amp(n) = -sum(saccade{i,s}.amplitude.true(1,ind));
        % horizontal component
        hor_sac_amp(n) = -sum(saccade{i,s}.amplitude.true(2,ind));
        else
        ver_sac_amp(n) = nan;
        hor_sac_amp(n) = nan;
        end
    end
    % remove bad trials
%     nanindx = (abs(th_err) > angerr_thresh) | (abs(hor_sac_amp) > angerr_thresh);
%     r_err(nanindx)=nan; th_err(nanindx)=nan;
    
    % compute correlations
    [r_corr(i,s),p_r(i,s)] = nancorr(r_err(:),ver_sac_amp(:));
    [th_corr(i,s),p_th(i,s)] = nancorr(th_err(:),hor_sac_amp(:));
    
    % pool subjects
    r_err_pool{s}(end+1:end+Ntrls) = r_err(:);
    th_err_pool{s}(end+1:end+Ntrls) = th_err(:);
    ver_sac_amp_pool{s}(end+1:end+Ntrls) = ver_sac_amp(:);
    hor_sac_amp_pool{s}(end+1:end+Ntrls) = hor_sac_amp(:);
    
    if 0
        subplot(2,Nstim,s); 
        plot(r_err,ver_sac_amp,'.','markersize',0.001,'color',colr(s,:)); hold on; xlabel('linear distance error [m]'); ylabel('ver. saccade sum');
        title(['rho=' num2str(r_corr(i,s)) ', p=' num2str(p_r(i,s))]); axis([-300 300 -30 30]); hold off;
        
        subplot(2,Nstim,s+Nstim);
        plot(th_err,hor_sac_amp,'.','markersize',0.001,'color',colr(s,:)); hold on; xlabel('angular error [m]'); ylabel('hor. saccade sum');
        title(['rho=' num2str(th_corr(i,s)) ', p=' num2str(p_th(i,s))]); axis([-30 30 -30 30]); hold off;
        
        b_r(s) = regress(hor_sac_amp(:),th_err(:));
    end
end
end
disp('------Angular Corr. for subjects separately:')
mean(th_corr)
sum(p_th<0.05)

disp('------Angular Corr. for pooled subjects:')
[r_corr_pool,p_r_pool] = cellfun(@(x,y) nancorr(x(:),y(:)),r_err_pool,ver_sac_amp_pool);
[th_corr_pool,p_th_pool] = cellfun(@(x,y) nancorr(x(:),y(:)),th_err_pool,hor_sac_amp_pool);
th_corr_pool
p_th_pool<0.05

% scatterplot pooled subjects
figure;
for s = 1:Nstim
subplot(2,Nstim,s); hold on;
plot(r_err_pool{s},ver_sac_amp_pool{s},'.','markersize',0.001,'color',colr(s,:)); xlabel('linear distance error [cm]'); ylabel('ver. saccade sum');
title(['rho=' num2str(r_corr_pool(s)) ', p=' num2str(p_r_pool(s))]); axis([-300 300 -30 30]); vline(0,'k');hline(0,'k');

subplot(2,Nstim,s+Nstim); hold on;
plot(th_err_pool{s},hor_sac_amp_pool{s},'.','markersize',0.001,'color',colr(s,:)); xlabel('angular error [deg]'); ylabel('hor. saccade sum');
title(['rho=' num2str(th_corr_pool(s)) ', p=' num2str(p_th_pool(s))]); axis([-30 30 -30 30]); vline(0,'k');hline(0,'k');

b_r(s) = regress(hor_sac_amp_pool{s}(:),th_err_pool{s}(:));
end

% bar graph
figure; hold on;
for s = 1:Nstim
[~,p_ttest(s)] = ttest(th_corr(:,s));
bar(s,mean(th_corr(:,s)),'edgecolor','none','facecolor',colr(s,:)); 
errorbar(s,mean(th_corr(:,s)),std(th_corr(:,s))./sqrt(Nsubs),'k','capsize',0)
xticks(1:Nstim); xticklabels({'vestibular','visual','combined'}); ylabel('corr. coefficient'); ylim([-0.1 0.5]);
title('sum of saccades vs steering errors');

if p_ttest(s) < 0.05; text(s,0.45,'*');
elseif p_ttest(s) < 0.01; text(s,0.45,'**');
elseif p_ttest(s) < 0.001; text(s,0.45,'***');
end

end

%% Early saturation of angular errors in vestibular condition (MIGHT NEED FOR SUPPLEMENTAL!!!) +++++++++++++++++++++++++
keepindx = [1 2 7 9 10 12 13 14];
subject = subject_backup(keepindx);
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = length(subject);
Nstim = size(poolindx,2);

maxlength = round(15/dt); % 15 seconds
straight_thresh = 10;
nanwin = 150;
nanthresh = 0.3;

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure;
th_ratio = []; h = []; p = [];
for i = 1:Nsubs
    clf;
    for s = 1:Nstim-1
        trlindx = poolindx{i,s};
        % endpoints
        x_sub = arrayfun(@(x) x.continuous.xmp(end)-x.continuous.xmp(1), subject(i).trials(trlindx));
        y_sub = arrayfun(@(x) x.continuous.ymp(end)-x.continuous.ymp(1), subject(i).trials(trlindx));
        [r_sub,th_sub] = cart2polarY(x_sub,y_sub);
        th_sub(abs(th_sub)<straight_thresh) = nan; % remove straight ahead targets
        
        % trajectories
        x_sub_traj = arrayfun(@(x) [x.continuous.xmp-x.continuous.xmp(1) ; nan(maxlength-numel(x.continuous.ts),1)], subject(i).trials(trlindx),'un',0);
        y_sub_traj = arrayfun(@(x) [x.continuous.ymp-x.continuous.ymp(1) ; nan(maxlength-numel(x.continuous.ts),1)], subject(i).trials(trlindx),'un',0);
        [r_sub_traj,th_sub_traj] = cellfun(@(x,y) cart2polarY(x,y),x_sub_traj,y_sub_traj,'un',0);
        th_sub_traj = cellfun(@(x,xe) nanify(x, find(abs(x(1:nanwin)) > nanthresh*xe)), th_sub_traj, num2cell(th_sub),'un',0); % correct for artifacts at trial start
        th_sub_traj = cellfun(@(x,xe) nanify(x, find(abs(x) > 1.2*xe)), th_sub_traj, num2cell(th_sub),'un',0); % clip for rest of artifacts
        
        % remove bad trials
        rmindx = cellfun(@(x) sum(isnan(x))/numel(x) > 0.4, th_sub_traj);
        th_sub_traj(rmindx) = [];
        th_sub(rmindx) = [];
        
        % take ratio
        th_ratio_tmp = cellfun(@(x,xe) abs(x(1:maxlength)/xe), th_sub_traj, num2cell(th_sub),'un',0);
    
        ts = (1:maxlength)*dt;
        
        hold on;
        shadedErrorBar(ts,nanmean([th_ratio_tmp{:}],2),nanstd([th_ratio_tmp{:}],[],2),'lineprops',{'color',colr(s,:)})
        xlabel('time [s]'); ylabel('\theta_{t}\\\theta_{final}');
        
        th_ratio{i,s} = th_ratio_tmp;
    end
    suptitle(subject(i).name);
end

figure; hold on; plot(th_ratio_50(:,1),th_ratio_50(:,2),'+','color',[.5 .5 .5]); plot(0:10,0:10,'k--');
xlabel('vestibular'); ylabel('visual'); title('rotation_{50}');axis equal;axis([2 8 2 8]);

th_ratio_mu = cellfun(@(x) nanmean([x{:}],2),th_ratio,'un',0);
th_ratio_50_indx = cellfun(@(x) find(x > 0.5, 1), th_ratio_mu);
th_ratio_50 = ts(th_ratio_50_indx);

[h,p] = ttest(th_ratio_50(:,1),th_ratio_50(:,2));

%%
%%
%% Repeat analysis without saccades
%%
%%
%% Remove saccades

if strcmp(subject_backup(5).name,'Ding'); subject_backup(5) = []; end
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
rmsac = 1;
%
% subject_nosac = SaccadeRemover(subject,prs);
%

% Analyze saccade-free eye movements
disp('Saccade Removal Analysis..........')

fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'tracking_regular_nosac.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_regular_nosac')
        load(fldcont(struct_indx).name);
        disp('.....tracking_regular_nosac Loaded.')
    else; disp('tracking_regular_nosac already exists.')
    end
else
    params = 'stimtype';
    tracking_regular_nosac = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular_nosac.mat','tracking_regular_nosac','-v7.3')
end

Nsubs = size(tracking_regular_nosac,1);
Nstim = size(tracking_regular_nosac,2);
dt = 1/60;
xaxislim = [0 10];  xaxislabel = 'time [s]';

%

% for i = 1:Nsubs
%     for s = 1:Nstim
%         indx = poolindx{i,s};
%         saccade_nosac{i,s} = AnalyseSaccades(subject_nosac(i).trials(indx),prs);
%     end
%     disp(['........Subject = ' num2str(i)])
% end

%% Re-sample quantities based on percentage of distance traveled
default_prs; 
align_flag = []; % 'align2targ';
dd_distperc = 0.005;


fldcont = dir(data_folder);
struct_indx = arrayfun(@(x) contains(x.name,'tracking_distperc_nosac.mat'), fldcont);
if any(struct_indx)
    if ~exist('tracking_distperc_nosac')
        load(fldcont(struct_indx).name);
        disp('.....tracking_distperc_nosac Loaded.')
    else; disp('tracking_distperc_nosac already exists.')
    end
else
    params = 'stimtype';
    tracking_distperc_nosac = EyeAnalysisOverDistancePercentage(tracking_regular_nosac,dd_distperc,align_flag,prs);
    save('tracking_distperc_nosac.mat','tracking_distperc_nosac','-v7.3')
end

xaxislim = [0 1];  xaxislabel = 'distance %';
%% Check average tracking errors difference (target vs stop position)

% over time
TrackingErrorsDifference(tracking_regular_nosac);
suptitle('Saccades Removed');

% over distance percentage
TrackingErrorsDifference(tracking_distperc_nosac);
suptitle('Saccades Removed');

%% Check differences in performance between vertical trackers and non-trackers

lateindx = 60; % percent
min_diff = 5;
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
subindx = FindVerticalTrackers(tracking_distperc_nosac(keepindx,:),lateindx,min_diff);

[bias,variab] = TrackersVsNonTrackersPerformance(subject(keepindx),subindx);

if 0 % hand-labeled trackers
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
temp = logical([1 0 0 0 0 1 0 0 1 0 1 0 0 1 ; 1 0 0 1 1 1 1 0 1 1 1 1 1 1 ; 1 0 0 1 1 1 1 0 1 1 1 1 1 1]');
Nsubs = numel(keepindx);
subindx = mat2cell(temp,Nsubs,ones(1,Nstim));

[bias,variab] = TrackersVsNonTrackersPerformance(subject(keepindx),subindx);
end

%% Ratio between SPTE and TTE
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);

stability_constant = 0;
[ratio_kernel,ts] = TTEvsSPTE_Ratio(tracking_distperc_nosac(keepindx,:),stability_constant);
suptitle('Saccades Removed');

%% Correlation between tracking error and Steering error

plevel = 0.1;
refplane = 'screen';
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];

[rho,pval] = SteeringVsTrackingError(tracking_distperc_nosac(keepindx,:),refplane);
suptitle('Saccades Removed');

[Rsignif,peakind] = PeakCorrSignificance1(rho,pval,plevel); nansum(Rsignif)

%% Position gain between eye and target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];
tracking = tracking_distperc_nosac(keepindx,:);

[g,ts] = EyePositionGain(tracking);
suptitle('Saccades Removed');

%% Compare position gain between target and actual/saccade-free eye positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);
keepindx = [1 2 7 9 10 12 13 14];

tracking = tracking_distperc(keepindx,:);
[g,ts,xaxis] = EyePositionGain(tracking); 
close;

tracking_nosac = tracking_distperc_nosac(keepindx,:);
[g_nosac,ts,xaxis] = EyePositionGain(tracking_nosac); 
close;

xaxis.lim = [0 .5];
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
figure('name','Position Gain w.r.t. target comparison (actual vs saccade-free)','numbertitle','off');
for s = 1:Nstim
    % horizontal
    subplot(2,Nstim,s); hold on;
    shadedErrorBar(ts,mean(g.hor.targ{s}),std(g.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts,mean(g_nosac.hor.targ{s}),std(g_nosac.hor.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    legend('actual eye','saccade-free');
    title('horizontal');
    % vertical
    subplot(2,Nstim,s+Nstim); hold on; 
    shadedErrorBar(ts,mean(g.ver.targ{s}),std(g.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts,mean(g_nosac.ver.targ{s}),std(g_nosac.ver.targ{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.5});
    hline(1,'k--'); xlabel(xaxis.label); ylabel('position gain'); axis([xaxis.lim 0 2]);
    legend('actual eye','saccade-free');
    title('vertical');
end

%% Correlation between eye and target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);

[rho,ts] = EyePositionCorr(tracking_distperc_nosac(keepindx,:));
suptitle('Saccades Removed');

%% Multiple regression between eye & target/stop positions

keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);

[multi,ts] = EyePositionMultiRegr(tracking_distperc_nosac(keepindx,:));
suptitle('Saccades Removed');

%% Regress endpoint of saccade against target/stop positions before saccade and plot weights over distance %
keepindx = setdiff(1:numel(subject_backup),[3 5 8 11]);

[rho,pval,beta,multi] = SaccadeEndpointRegression(saccade_nosac(keepindx,:));
suptitle('Saccades Removed');

%% Comparison of SE and TE correlation w and w/o saccades
thresh=0.7; % average up to 70% of trial
tracking = tracking_distperc(keepindx,:);
tracking_nosac = tracking_distperc_nosac(keepindx,:);
Nsubs = size(tracking,1);
Nstim = size(tracking,2);
Nt = numel(tracking{1}.eyepos.errcorr.targ.screen.r); %round(9/dt);

for i = 1:Nsubs
   for s = 1:Nstim 
       corr_sac(i,s) = max(tracking{i,s}.eyepos.errcorr.targ.screen.r(1:thresh*Nt));
       corr_nosac(i,s) = max(tracking_nosac{i,s}.eyepos.errcorr.targ.screen.r(1:thresh*Nt));
   end
end

% unity line plot
figure; hold on;
for s = 1:Nstim-1
    plot(corr_sac(:,s),corr_nosac(:,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none');plot(-1:1,-1:1,'k--');
end
xlabel({'Corr. coef Average','(actual eye)'}); ylabel({'Corr. coef Average','(saccade-free)'});title(['0-' num2str(thresh*100) '% distance'])
axis([0 0.9 0 0.9])

% open/closed bar plot
% figure; hold on;
% for s = 1:Nstim-1
%     bar(s,mean(corr_sac(:,s)),'facecolor','none','edgecolor',[.5 .5 .5]); errorbar(s,mean(corr_sac(:,s)),std(corr_sac(:,s))./sqrt(Nsubs),'color',[.5 .5 .5],'capsize',0);
%     bar(s,mean(corr_nosac(:,s)),'facecolor',colr(s,:),'edgecolor','none'); errorbar(s,mean(corr_nosac(:,s)),std(corr_nosac(:,s))./sqrt(Nsubs),'k','capsize',0);
% end

% bar plot of correlation difference
N=100000;
figure; hold on;
for s = 1:Nstim-1
    bar(s,mean(corr_sac(:,s)-corr_nosac(:,s)),'facecolor',colr(s,:),'edgecolor','none'); errorbar(s,mean(corr_sac(:,s)-corr_nosac(:,s)),std(corr_sac(:,s)-corr_nosac(:,s))./sqrt(Nsubs),'k','capsize',0);
    % bootstrap
    tmp = corr_sac(:,s)-corr_nosac(:,s);
    randindx = randi(Nsubs,[Nsubs N]);
    tmp = nanmean(tmp(randindx));
    p(s) = sum(tmp <= 0)/N;
end
ylabel('corr. coef. difference'); xticks(1:2);xticklabels({'vestibular','visual'});
text(1.2,0.15,'*');
text(2.2,0.15,'**');

% connecting lines
figure; 
for s = 1:Nstim-1
    subplot(1,2,s);hold on;
    plot(1:2,[corr_sac(:,s) corr_nosac(:,s)],'o-','color',colr(s,:),'markerfacecolor',colr(s,:),'markeredgecolor','none'); 
    axis([0 3 0 1]);xticks(1:2); xticklabels({'actual eye','saccade-free'});ylabel('corr. coef. average');
end

%%
%%
%%
%% Relative tracking index (RTI)
figure('name','All Subjects: Relative Tracking Index','numbertitle','off','position',[0 700 900 250]);
RTI = [];
for i = 1:Nsubs
    for s = 1:Nstim
        ts = tracking{i,s}.misc.ts;
        TTE = tracking{i,s}.eyepos.error.targ.all;
        SPTE = tracking{i,s}.eyepos.error.stop.all;
        
        RTI_temp = TTE./(SPTE+TTE);
        RTI{s}(i,:) = nanmean(RTI_temp); % RTI_sd = nanstd(RTI_temp);        
    end
end

for s = 1:Nstim
    subplot(1,Nstim,s); hold on;
    plot(ts,RTI{s},'color',colr(s,:)); 
    shadedErrorBar(ts,nanmean(RTI{s}),nanstd(RTI{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    axis([xaxislim 0 1]); hline(0.5,'k--');
    xlabel(xaxislabel); ylabel(['Relative Tracking Index']);
    
end

%% Target-tracking index comparison (target vs stop position)

default_prs;
[varexp_start,varexp_stop] = TargetVsStopTTI(tracking_regular,prs);

%% Is it percentage or absolute distance to end when subjects stop tracking?
% check variability of data series
SD = []; ts = [];
for i = 1:Nsubs
    ver = figure('name',[subject(i).name ': Vertical data series variability'],'numbertitle','off','position',[0 550 900 400]);
    hor = figure('name',[subject(i).name ': Horizontal data series variability'],'numbertitle','off','position',[0 50 900 400]);
    for s = 1:Nstim
        Tlength = 12*60;
        ts.time = tracking{i,s}.misc.ts(1:Tlength);
        ts.distperc = tracking_distperc{i,s}.misc.ds;
        D2Elength = floor(500/dd_dist2end) - 1;
        ts.dist2end = tracking_dist2end{i,s}.misc.ds(end-D2Elength:end);

        % vertical
        SD.ver.time = nanstd(tracking{i,s}.eyepos.screen.ver_mean(:,1:Tlength));
        SD.ver.distperc = nanstd(tracking_distperc{i,s}.eyepos.screen.ver_mean);
        SD.ver.dist2end = nanstd(tracking_dist2end{i,s}.eyepos.screen.ver_mean(:,end-D2Elength:end));
        figure(ver);
        subplot(1,Nstim,s); hold on;
        plot((ts.time)./max(ts.time),SD.ver.time,'color',colr(s,:),'linewidth',2); ylim([0 20]);
%         subplot(3,Nstim,s+Nstim); hold on;
        plot(ts.distperc./max(ts.distperc),SD.ver.distperc,'color',colr(s,:)*0.6,'linewidth',2); ylim([0 20]);
%         subplot(3,Nstim,s+2*Nstim); hold on;
        plot(flip(ts.dist2end./min(ts.dist2end)),SD.ver.dist2end,'color',colr(s,:)*0,'linewidth',2); ylim([0 20]);
        
        % horizontal
        targ_sign = sign(tracking{i,s}.eyepos.screen.hor_mean(:,1));
        SD.hor.time = nanstd(tracking{i,s}.eyepos.screen.hor_mean(:,1:Tlength).*targ_sign);
        SD.hor.distperc = nanstd(tracking_distperc{i,s}.eyepos.screen.hor_mean.*targ_sign);
        SD.hor.dist2end = nanstd(tracking_dist2end{i,s}.eyepos.screen.hor_mean(:,end-D2Elength:end).*targ_sign);
        
        figure(hor);
        subplot(1,Nstim,s); hold on;
        plot(ts.time./max(ts.time),SD.hor.time,'color',colr(s,:),'linewidth',2); ylim([0 20]);
%         subplot(3,Nstim,s+Nstim); hold on;
        plot(ts.distperc./max(ts.distperc),SD.hor.distperc,'color',colr(s,:)*0.6,'linewidth',2); ylim([0 20]);
%         subplot(3,Nstim,s+2*Nstim); hold on;
        plot(flip(ts.dist2end./min(ts.dist2end)),SD.hor.dist2end,'color',colr(s,:)*0,'linewidth',2); ylim([0 20]);


    end
end

%% Saccade rate as a function of time, distance %, and distance to end
sacrate = [];
for i = 1:Nsubs
    % figure('name',[subject(i).name ': Saccade Rate'],'numbertitle','off','position',[0 450 900 500]);
    for s = 1:Nstim
        ntrls = numel(unique(tracking_regular{i,s}.saccade.trl));
        numsac = numel(tracking_regular{i,s}.saccade.trl);
        % time
        nt = -0.2:0.2:20;
        [ny,nt] = hist(tracking_regular{i,s}.saccade.time,nt);
        ny = ny/numsac;
        sacrate.time.val{s}(i,:) = ny(2:end-1); sacrate.time.bins{s} = nt(2:end-1);  
        % distance percentage
        ndp = -0.005:0.01:1.005;
        [ny,ndp] = hist(tracking_distperc{i,s}.saccade.time,ndp);
        ny = ny/numsac;
        sacrate.distperc.val{s}(i,:) = ny(2:end-1) ; sacrate.distperc.bins{s} = ndp(2:end-1);  
        % distance to end
        nd2e = -600:6:3;
        [ny,nd2e] = hist(tracking_dist2end{i,s}.saccade.time,nd2e);
        ny = ny/numsac;
        sacrate.dist2end.val{s}(i,:) = ny(2:end-1) ; sacrate.dist2end.bins{s} = nd2e(2:end-1);  
        
        if 0
        subplot(3,Nstim,s); hold on;
        plot(sacrate.time.bins{s},sacrate.time.val{s}(i,:),'color',colr(s,:));
        xlabel('time [s]'); ylabel('Prob. of saccades'); ylim([0 0.2]);
        subplot(3,Nstim,s+Nstim); hold on;
        plot(sacrate.distperc.bins{s},sacrate.distperc.val{s}(i,:),'color',colr(s,:)*0.7);
        xlabel('distance %'); ylabel('Prob. of saccades'); ylim([0 0.2]);
        subplot(3,Nstim,s+2*Nstim); hold on;
        plot(sacrate.dist2end.bins{s},sacrate.dist2end.val{s}(i,:),'color',colr(s,:)*0.3);
        xlabel('distance to end [cm]'); ylabel('Prob. of saccades'); ylim([0 0.2]);
        end
    end
end

figure('name','Saccade Rate','numbertitle','off','position',[0 450 900 500]);
for s = 1:Nstim
    
    subplot(3,Nstim,s); hold on;
    shadedErrorBar(sacrate.time.bins{s},mean(sacrate.time.val{s}),std(sacrate.time.val{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)});
    xlabel('time [s]'); ylabel('Prob. of saccades'); ylim([0 0.2]);
    
    subplot(3,Nstim,s+Nstim); hold on;
    shadedErrorBar(sacrate.distperc.bins{s},mean(sacrate.distperc.val{s}),std(sacrate.distperc.val{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.7});
    xlabel('distance %'); ylabel('Prob. of saccades'); ylim([0 0.2]);
    
    subplot(3,Nstim,s+2*Nstim); hold on;
    shadedErrorBar(sacrate.dist2end.bins{s},mean(sacrate.dist2end.val{s}),std(sacrate.dist2end.val{s})./sqrt(Nsubs),'lineprops',{'color',colr(s,:)*0.3});
    xlabel('distance to end [cm]'); ylabel('Prob. of saccades'); ylim([0 0.2]);

    
end


%%
%%
%% RNN info for Kaushik
%%
%% Correlation between target and response
nperm = 1000;
maxrewardwin = 600;

tau = []; auc = [];
for i = 1:Nsubs
    for s = 1:Nstim
        
        tau{i} = arrayfun(@(x) x.prs.tau, subject(i).trials);
        velctrlindx{i} = tau{i} < 1.5;
        trlindx = intersect(poolindx{i,s},find(velctrlindx{i}));
        
        th_sub = arrayfun(@(x) x.prs.th_sub, subject(i).trials(trlindx));
        th_tar = arrayfun(@(x) x.prs.th_tar, subject(i).trials(trlindx));
        rho_th(i,s) = nancorr(th_tar(:),th_sub(:));
        
        r_sub = arrayfun(@(x) x.prs.r_sub, subject(i).trials(trlindx));
        r_tar = arrayfun(@(x) x.prs.r_tar, subject(i).trials(trlindx));
        rho_r(i,s) = nancorr(r_tar(:),r_sub(:));
        
        % ROC
        x_tar = arrayfun(@(x) x.prs.fireflyposx, subject(i).trials(trlindx));
        y_tar = arrayfun(@(x) x.prs.fireflyposy, subject(i).trials(trlindx));
        x_sub = arrayfun(@(x) x.continuous.xmp(end) - x.continuous.xmp(1), subject(i).trials(trlindx));
        y_sub = arrayfun(@(x) x.continuous.ymp(end) - x.continuous.ymp(1), subject(i).trials(trlindx));
        b_x = regress(x_sub(:),x_tar(:));
        b_y = regress(y_sub(:),y_tar(:));
                
        X_monk = [x_sub(:)/b_x   y_sub(:)/b_y];
        X_fly = [x_tar(:)   y_tar(:)];
        [rewardwin, pCorrect, pCorrect_shuffled_mu] = ComputeROCFirefly(X_fly,X_monk,maxrewardwin,nperm);
%         auc(i,s) = sum(pCorrect)*diff(rewardwin(1:2));
        auc(i,s) = sum(0.5 * (pCorrect(2:end) + pCorrect(1:end-1)) .* diff(pCorrect_shuffled_mu));
    end
end

disp('+++++ rho_theta +++++')
for s = 1:Nstim
disp([legend_input{s} ': ' num2str(mean(rho_th(:,s))) ' ' char(177) ' ' num2str(std(rho_th(:,s))/sqrt(Nsubs))])
end

disp(' ')
disp('+++++ rho_r +++++')
for s = 1:Nstim
disp([legend_input{s} ': ' num2str(mean(rho_r(:,s))) ' ' char(177) ' ' num2str(std(rho_r(:,s))/sqrt(Nsubs))])
end

% Plot AUC

figure; hold on;
for s = 1:Nstim
    bar(s,nanmean(auc(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(auc(:,s)),nanstd(auc(:,s))/sqrt(Nsubs),'k','capsize',0);
end
xticks([1:Nstim]); xticklabels(legend_input);ylabel('AUC');

disp(['+++++ AUC mean ' char(177) ' SEM +++++'])
for s = 1:Nstim
disp([legend_input{s} ': ' num2str(mean(auc(:,s))) ' ' char(177) ' ' num2str(std(auc(:,s))/sqrt(Nsubs))])
end


%%
%%
%% ADDITIONAL ANALYSIS OF SACCADE RELATIONSHIP WITH END OF STEERING
%%
%%

StopTrackingAnalysis;

