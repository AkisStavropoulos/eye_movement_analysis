%% Dynamical tau Preliminary analysis script
% extract all data from current folder


% clear
data_folder = 'C:\Users\ges6\Documents\MATLAB\Data\Single firefly with motion cuing\Subjects\New\';
cd(data_folder);
folders = dir(data_folder);   folders = folders([folders.isdir]);

indx = cell2mat(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders,'uniformoutput',false));
subject_name = {folders(indx).name};

for i = 1:numel(subject_name)
    
    cd([data_folder subject_name{i}]);

    default_prs;
    subject(i).trials = AddTrials2Behaviour(prs); % includes some fixes
    if strcmp(subject_name{i},'Crystal')
        FaultySpike2PositionFix(subject(i).trials); % file 559 is corrupt after trial 30
    end
    
    saveinput = input('Save? 0 or 1:    ');
    if saveinput
        trials = subject(i).trials;
        save(subject_name{i},'trials'); % re-run Aaron, deleted by accident !
        disp(['Data saved: ' subject_name{i}]);
    end
end


% % Sorting trials based on conditions (stim type, JS coefficient, intersection of both)
% limtype = 'hard';
% bins = 3;
% [stimtype,tau,stimtau] = condsorting_bintau(trials,bins,limtype);
% stnames = fieldnames(stimtype); % stimtype
% taunames = fieldnames(tau); % taus
% % save distance and angular error for every trial
% [r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,0,0);
% trials = save_errors(trials,r_tar,r_sub,theta_tar,theta_sub);
% Save data
%% load data
subject_name = 'Eleo';
DataFolderRead;
cd([data_folder subject_name]);
load(subject_name);
%% check GIA error and moog position / plot trajectories
thresh = 0.07; % find trials over the threshold
plt = 0;
indx_err = find_trl_GIA_error(trials,thresh,plt);
%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
scatterDistAng(trials,0,1);

% different color for each stimtype
scatterDistAng2(trials,stimtype,0,1);

scatterDistAng3(trials,tau,0,1);

scatterDistAng4(trials,stimtau,0,1);
%% Plot radial and angular errors as a function of TAU
% all trials
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,0,0);
errors_f_tau(trials,[],r_tar,r_sub,theta_tar,theta_sub);
% stimtype
errors_vs_tau(trials,stimtype); % run save_errors first
% [r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,0,0);
% errors_f_tau(trials,stimtype,r_tar,r_sub,theta_tar,theta_sub);
%% final error as a function of the initial target DISTANCE 
% no sign - all trials
dist2target(trials,1,0); 
% with sign - all trials
errors_vs_dist(trials,[]);

% no sign - stimtype
dist2target2(trials,stimtype,0,0);
% sign - stimtype
errors_vs_dist(trials,stimtype);

% no sign - tau
[r_tar,err] = dist2target3(trials,tau,0,0);
% sign - tau
errors_vs_dist(trials,tau);

% no sign - stimtau
[r_tar,err] = dist2target4(trials, stimtau,0,0);
% sign - stimtau
errors_vs_dist(trials,stimtau);

%% final error as a function of the initial target ANGLE
% all trials
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,0,0);
errors_vs_ang(trials,[],r_tar,r_sub,theta_tar,theta_sub)

% stimtype
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,0,0);
errors_vs_ang(trials,stimtype,r_tar,r_sub,theta_tar,theta_sub)

% tau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng3(trials,tau,0,0);
errors_vs_ang(trials,tau,r_tar,r_sub,theta_tar,theta_sub)

% stimtau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng4(trials,stimtau,0,0);
errors_vs_ang(trials,stimtau,r_tar,r_sub,theta_tar,theta_sub)

%% final error as a function of the initial target EXPERIENCE
% all trials
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,0,0);
errors_vs_exp(trials,[],r_tar)

% stimtype
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,0,0);
errors_vs_exp(trials,stimtype,r_tar)

% tau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng3(trials,tau,0,0);
errors_vs_exp(trials,tau,r_tar)

% stimtau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng4(trials,stimtau,0,0);
errors_vs_exp(trials,stimtau,r_tar)

%% plot the mean UNSIGNED error
% choose one
[r_tar,err] = dist2target2(trials,stimtype,0,0);
[r_tar,err] = dist2target3(trials,tau,0,0);
[r_tar,err] = dist2target4(trials, stimtau, 0,0);

binsize = 75;
[mean_err,std_err] = meanbin_unsigned_error(trials,r_tar,err,binsize,1);

%% mean SIGNED error
binsize = 75;
binsize_th = 5;
% all trials
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,0,0);
[r_err,theta_err] = signed_error(r_tar,r_sub,theta_tar,theta_sub);
[mean_err,std_err] = meanbin_signed_error(trials,r_tar,r_err,binsize,1); % distance
[mean_err,std_err] = meanbinerror_th(trials,theta_tar,theta_err,binsize_th,1); % angle
% stimtype
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,0,0);
[r_err,theta_err] = signed_error(r_tar,r_sub,theta_tar,theta_sub);
[mean_err,std_err] = meanbin_signed_error(trials,r_tar,r_err,binsize,1);
[mean_err,std_err] = meanbinerror_th(trials,theta_tar,theta_err,binsize_th,1);

% tau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng3(trials,tau,0,0);
[r_err,theta_err] = signed_error(r_tar,r_sub,theta_tar,theta_sub);
[mean_err,std_err] = meanbin_signed_error(trials,r_tar,r_err,binsize,1);
[mean_err,std_err] = meanbinerror_th(trials,theta_tar,theta_err,binsize_th,1);

% stimtau
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng4(trials,stimtau,0,0);
[r_err,theta_err] = signed_error(r_tar,r_sub,theta_tar,theta_sub);
[mean_err,std_err] = meanbin_signed_error(trials,r_tar,r_err,binsize,1);
[mean_err,std_err] = meanbinerror_th(trials,theta_tar,theta_err,binsize_th,1);

%% epsilon: error normalized by initial distance
% for modalities

% [r_tar,err] = dist2target2(trials,stimtype,0,0);
% eps = epsilonerror1(err,r_tar);
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,0,0);
eps = epsilonerror(r_sub,r_tar);
% for taus

% [r_tar,err] = dist2target3(trials,tau,0,0);
% eps = epsilonerror1(err,r_tar);
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng3(trials,tau,0,0);
eps = epsilonerror(r_sub,r_tar);
% stimtau

% [r_tar,err] = dist2target4(trials, stimtau, 0,0);
% eps = epsilonerror1(err,r_tar);
[r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng4(trials,stimtau,0,0);
eps = epsilonerror(r_sub,r_tar);
%% bias and variance for each condition
[errvec,tarloc] = errors2d(trials,stimtype);
[errvec,tarloc] = errors2d(trials,tau);
[errvec,tarloc] = errors2d(trials,stimtau);

% choose the targets that correspond to the conditions you are interested

% compute bias and variance across conditions
binsize = 100;
[meanbias,sigma,~] = biasvardecomp(trials,errvec,tarloc,binsize,1);
%% Plot trajectories for all conditions and modalities

plottraj1(trials,stimtype)
plottraj2(trials,tau)
%% History analysis
% requires to run first the section "save distance and angular error for every trial"
history_depth = 3;
bins = 5;
[trials,delta,deltas] = comp_delta(trials,bins,history_depth);
errors_f_delta(trials,deltas); % 1 trial history

% multilinear regression
plot_mean_std_history_error(trials,stimtype);


%% Joystick analysis
JoystickAnalysis;

plotcontrol2d(trials,tau);

plotvel2D(trials,tau);

plotallovel2D(trials,tau);




%% Check GIA error with performance

stim_err_behv_err_plot(trials,stimtype);

%% fit slow prior model and plot data

[opt_coefs,opt_fval,models] = fit_slow_prior(trials,stimtype,'explore');

[r_tar,r_sub,theta_tar,theta_sub] = plot_opt_scatter(trials,models,stimtype,1,1);

%% plot biases (slope) for data and model
[~,~,~,~,b_r,b_th] = scatterDistAng2(trials,stimtype,1);
[~,~,~,~,b_r_model,b_th_model] = plot_opt_scatter(trials,models,stimtype,1);

plot_bias(trials,b_r,b_th,b_r_model,b_th_model);
%% Linear regression fitting of distance and angle scatterplots
% will probably need to be fixed
for n = 1:length(stnames)-1
    clear x;clear y;
    x = r_tar.(stnames{n})';
    y = r_sub.(stnames{n})';
    [b_r.(stnames{n}), a_r.(stnames{n}), bint_r.(stnames{n}), aint_r.(stnames{n}), r_r.(stnames{n}), p_r.(stnames{n})]...
        = regress_perp(x,y,0.05,2);
end

[b_theta, a_theta, bint_theta, aint_theta, r_theta, p_theta] = regress_perp(theta_tar,theta_sub,0.05,1);
