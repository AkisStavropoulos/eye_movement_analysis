%% Extract data
% extract all data from current folder
clear
default_prs;
trials = AddTrials2Behaviour(prs);

% FaultySpike2PositionFix;

% TestMoogDelivery;
%% Sorting trials based on conditions (stim type, JS coefficient, intersection of both)
[stimtype,jscoef,stimjs] = condsorting(trials);
stnames = fieldnames(stimtype); % stimtype
coefnames = fieldnames(jscoef); % JS coefficient
condnames = fieldnames(stimjs); % intersection

%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
scatterDistAng(trials);

% different color for each stimtype
scatterDistAng2(trials,stimtype);

% different color for each JS coefficient
scatterDistAng3(trials,jscoef);

% different color for each modality given condition
scatterDistAng4(trials,stimjs,0,1);

%% final error as a function of the initial distance of the target
dist2target(trials,1);

dist2target2(trials,stimtype,1);

dist2target3(trials,jscoef,1);

dist2target4(trials,stimjs,1);
%% plot the mean error
% choose one
[r_tar,err] = dist2target2(trials,stimtype);
[mean_err,std_err] = meanbinerror(r_tar,err,1);

[r_tar,err] = dist2target3(trials,jscoef);
[mean_err,std_err] = meanbinerror(r_tar,err,1);

[r_tar,err] = dist2target4(trials,stimtype,jscoef);
[mean_err,std_err] = meanbinerror(r_tar,err,1);

%% epsilon: error normalized by initial distance
% for modalities
[r_tar,err] = dist2target2(trials,stimtype);
[mean_err,std_err] = meanbinerror(r_tar,err);
eps = epsilonerror(err,r_tar);
% for coefficients
[r_tar,err] = dist2target3(trials,jscoef);
[mean_err,std_err] = meanbinerror(r_tar,err);
eps = epsilonerror(err,r_tar);
% for all conditions
[r_tar,err] = dist2target4(trials,stimtype,jscoef);
[mean_err,std_err] = meanbinerror(r_tar,err);
eps = epsilonerror(err,r_tar);
close all;
%% bias and variance for each condition
[errvec,tarloc] = errors2d(trials,stimjs);
% choose the targets that correspond to the conditions you are interested

% compute bias and variance across conditions
[meanbias,sigma,~] = biasvardecomp(errvec,tarloc,1);


%% Plot trajectories for all conditions and modalities

plottraj1(trials,stimtype)

plottraj2(trials,jscoef)

%% fit slow prior model and plot data

[opt_coefs,opt_fval,models] = fit_slow_prior(trials,stimjs,'explore');

[r_tar,r_sub,theta_tar,theta_sub] = plot_opt_scatter(trials,models,stimjs);

%% plot biases (slope) for data and model
[~,~,~,~,b_r,b_th] = scatterDistAng4(trials,stimjs,1);
[~,~,~,~,b_r_model,b_th_model] = plot_opt_scatter(trials,models,stimjs,1);

plot_bias(trials,b_r,b_th,b_r_model,b_th_model);

%% Linear regression fitting of distance and angle scatterplots
% will probably need to be fixed
for n = 1:length(condnames)-1
    clear x;clear y;
    x = r_tar.(condnames{n})';
    y = r_sub.(condnames{n})';
    [b_r.(condnames{n}), a_r.(condnames{n}), bint_r.(condnames{n}), aint_r.(condnames{n}), r_r.(condnames{n}), p_r.(condnames{n})]...
        = regress_perp(x,y,0.05,2);
end

[b_theta, a_theta, bint_theta, aint_theta, r_theta, p_theta] = regress_perp(theta_tar,theta_sub,alpha,option);

%% Linear regression and ROC analysis

% stats = AnalyseBehaviour(trials,prs);

%% fit kausik's curve on radius and scatterplot
scatterDistAng(trials,1);

% different color for each stimtype
scatterDistAng2(trials,stimtype,1);

% different color for each JS coefficient
scatterDistAng3(trials,jscoef,1);

% different color for each modality given condition
scatterDistAng4(trials,stimtype,jscoef,1);

%% deconvolve linear and angular velocity

DeconvolveLinearJoystickInput;

DeconvolveAngularJoystickInput;

