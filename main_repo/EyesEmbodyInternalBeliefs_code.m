%%
%%
%% Eyes embody dynamic internal beliefs during visual and inertial navigation (Stavropoulos, Lakshminarasimhan, Angelaki)
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

default_prs;

% Load data
fldcont = dir(data_folder);
dataindx = arrayfun(@(x) contains(x.name,'EyeMovementDataset.mat'), fldcont);
if any(dataindx)
    tic;
    disp('.... Loading EyeMovementDataset.mat')
    load(fldcont(dataindx).name);
    toc;
else
    error('dataset was not found.')
end

%% Extract time-based eye movement analysis
rmsac = 0;

% sensory condition groups
if ~exist('tracking_regular')
    params = 'stimtype';
    tracking_regular = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular.mat','tracking_regular','-v7.3');
end

% all trials
if ~exist('tracking_regular_all')
    params = [];
    tracking_regular_all = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular_all.mat','tracking_regular_all','-v7.3');
end

Nsubs = size(tracking_regular,1);
Nstim = size(tracking_regular,2);
dt = 1/60;
xaxislim = [0 10];  xaxislabel = 'time [s]';

%% Re-sample quantities based on percentage of distance traveled
align_flag = []; % 'align2targ';
dd_distperc = 0.005;

if ~exist('tracking_distperc')
    tracking_distperc = EyeAnalysisOverDistancePercentage(tracking_regular,dd_distperc,align_flag,prs);
    save('tracking_distperc.mat','tracking_distperc','-v7.3')
end

xaxislim = [0 1];  xaxislabel = 'distance %';

%%
%% BEHAVIORAL PERFORMANCE FIGURES
%%
%% Bias Detection and Multiplicative Model Fit (+ compare 1st-button-push bias with end-of-trial bias)

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
for i = 1%:length(subject)
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

%% Regression of actual vs predicted eye positions when target OFF
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
outlier_thresh = 25;
subindx = keepindx;
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
    arrayfun(@(s) shadedErrorBar(ts,nanmedian(err_all_distperc{s}./nanmean(err_all_time_shuffled)),nanstd(err_all_distperc{s})./sqrt(Nsubs)./nanmean(err_all_time_shuffled),'lineprops',{'color',colr(s,:)}), 1:Nstim);
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

%% Correlation between tracking error and Steering error

plevel = 0.05;
refplane = 'screen';

% [rho,pval] = SteeringVsTrackingError(tracking_regular(keepindx,:),refplane);
[rho,pval] = SteeringVsTrackingError(tracking_distperc(keepindx,:),refplane);

[Rsignif,peakind] = PeakCorrSignificance1(rho,pval,plevel); nansum(Rsignif)

%% Multiple regression between eye & target/stop positions

[multi,ts] = EyePositionMultiRegr(tracking_distperc(keepindx,:));

% Average of both components' coefficients
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple

figure('name','Multiple Regression (average)','numbertitle','off');
for s = 1:Nstim
    % average
    Nt = numel(ts);
    b_targpos_mu = movmedian(mean([mean(multi.hor.targ{s}) ; mean(multi.ver.targ{s})]),10); b_targpos_mu(Nt-15:Nt) = nan;
    b_targpos_se = mean([std(multi.hor.targ{s})/sqrt(Nsubs) ; std(multi.ver.targ{s})/sqrt(Nsubs)]); b_targpos_se(Nt-15:Nt) = nan;
    
    b_stoppos_mu = movmedian(mean([mean(multi.hor.stop{s}) ; mean(multi.ver.stop{s})]),10); b_stoppos_mu(Nt-15:Nt) = nan;
    b_stoppos_se = mean([std(multi.hor.stop{s})/sqrt(Nsubs) ; std(multi.ver.stop{s})/sqrt(Nsubs)]); b_stoppos_se(Nt-15:Nt) = nan;

    subplot(1,Nstim,s); hold on;
    shadedErrorBar(ts,b_targpos_mu,b_targpos_se,'lineprops',{'color',colr(s,:)});
    shadedErrorBar(ts,b_stoppos_mu,b_stoppos_se,'lineprops',{'color',colr(s,:)*0.5});
    xlabel(xaxislabel); ylabel('regression coef.'); axis([xaxislim -0.2 1.2]);   hline(1,'k--');
    legend('target position','stop position');
    title('average');
end

%% Generate target belief over time
dt=1/60;
tracking = tracking_distperc(keepindx,:); %tracking_distperc;
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

%%
%% New Saccade Analysis
%%
%%
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
    else; error('"tracking_distperc" does not exist.')
    end
end

%% Saccade Probability
saccade = saccade_all(keepindx,:);

SaccadeProbability(saccade);

%% Saccade amplitude for different trial epochs (Target ON, Steering, End of Trial)
saccade = saccade_all(keepindx,:);

SaccadeAmplitude_Epochs(saccade);

%% Regression KERNEL of saccade amplitude vs tracking errors (target, stop, belief) +++++++IN ANALYSIS+++++++
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

params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
rmsac = 1;

% Analyze saccade-free eye movements

if ~exist('tracking_regular_nosac')
    params = 'stimtype';
    tracking_regular = EyeAnalysisOverTime(subject,params,rmsac,prs);
    save('tracking_regular_nosac.mat','tracking_regular_nosac','-v7.3');
end

Nsubs = size(tracking_regular_nosac,1);
Nstim = size(tracking_regular_nosac,2);
dt = 1/60;
xaxislim = [0 10];  xaxislabel = 'time [s]';

%% Re-sample quantities based on percentage of distance traveled
default_prs; 
align_flag = []; % 'align2targ';
dd_distperc = 0.005;


if ~exist('tracking_distperc_nosac')
    tracking_distperc = EyeAnalysisOverDistancePercentage(tracking_regular,dd_distperc,align_flag,prs);
    save('tracking_distperc_nosac.mat','tracking_distperc_nosac','-v7.3')
end

xaxislim = [0 1];  xaxislabel = 'distance %';


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
