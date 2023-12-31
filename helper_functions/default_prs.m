% function prs = default_prs(monk_id,session_id)

% if nargin<2, session_id = 1; end

%% session specific parameters
% monkeyInfoFile_joysticktask;
% monkeyInfo = monkeyInfo([monkeyInfo.session_id]==session_id & [monkeyInfo.monk_id]==monk_id);
% prs.filepath_behv = ['C:\Users\jklakshm\Documents\Data\firefly-monkey\' monkeyInfo.folder '\behavioural data\'];
% prs.maxchannels = max(monkeyInfo.channels);
% prs.coord = monkeyInfo.coord;
% prs.units = monkeyInfo.units;
% prs.comments = monkeyInfo.comments;

%% data acquisition parameters
prs.recalibrate = 0; % recalibrate eye positions offline
prs.rec_thresh = 3; % variability of eye signal threshold
prs.blink_thresh_pos = 60;
prs.blink_thresh_vel = 25;
prs.nanpadding = 100;
prs.fixateduration = .75;
prs.fs_smr = 5000/6; % sampling rate of smr file
prs.filtwidth = 10; % width in samples (10 samples @ fs_smr = 10x0.0012 = 12 ms)
prs.filtsize = 10*prs.filtwidth; % size in samples
prs.factor_downsample = 13; % select every nth sample
prs.dt = 5/(prs.fs_smr);
prs.screendist = 32.5; %cm
prs.height = 70; %cm
prs.interoculardist = 6; %cm
prs.framerate = 60; %(sec)^-1
prs.x0 = 0; % x-position at trial onset (cm)
prs.y0 = -32.5; %y-position at trial onset (cm)
prs.jitter_marker = 0.25; % variability in marker time relative to actual event time (s)
prs.minpeakprominence.monkpos = 10; % expected magnitude of change in monkey position during teleportation (cm)
prs.minpeakprominence.flypos = 1; % expected magnitude of change in fly position in consecutive trials (cm)
prs.eyechannels = [2 2];

%% static stimulus parameters
prs.monk_startpos = [0 -45];
prs.fly_ONduration = 1; % 1 sec ??
prs.fixation_start = 1; % time that subjects are already fixating at the target

%% data analysis parameters
% behavioural analysis
prs.alldeg = 1; % extract both tracking and behavioural error in degrees
prs.TauAdjust = 0; % adjust frame cutoffs to varying duration of stopping time caused by varying taus
prs.Ntaus = 2; % number of taus for which cut-off happens
prs.PeakVel_cutoff = 0; % Target-tracking error cut-off based on Peak velocity
prs.Dist_cutoff = 0; % Target-tracking error cut-off based on distance percentage
prs.Dist2End = 100; % Target-tracking error cut-off based on Distance to the end
prs.DistPerc = .7; % percentile of Distance at which Target-tracking error cut-off happens
prs.boots = 0; % bootstrap 1 or not 0
prs.bootstrap_trl = 100; % number of trials to bootstrap
prs.norm = 0; % normalize TTI by initial value of TTI
prs.StopPos = 0;
prs.ngroups = 2; % accuracy magnitude groups
prs.mintrialsforstats = 50; % need at least 100 trials for stats to be meaningful
prs.npermutations = 50; % number of permutations for trial shuffled estimates

prs.medfiltwidth = 100; % width of median filter removing one sample spikes
prs.sgfiltwidth = 61; % width of Savitsky-Golay filter smoothing eye signal
prs.saccade_extraction = 0;
prs.saccade_thresh = 42; % deg/s
prs.trl_cutoff = .5;  % seconds
prs.lambda = 0.7; %1.5; % lambda: acceleration threshold, how many std's away from mean
prs.min_sac_mag = 0.5; % degrees, minimum saccade magnitude
prs.saccade_duration = [0.05 0.3]; %seconds [min max]
prs.sac_refr_period = .05; % seconds
prs.overlap_pct = [0 100]; % overlap of pre- and post-saccade percentiles of eye position
prs.peri_sac_dur = 0.24; % length of pre- and post-saccade position calculation (seconds)
prs.max_MSE_cal_fit = 45;
prs.v_thresh = 5; % cm/s
prs.v_time2thresh = 0.05; % (s) approx time to go from zero to threshold or vice-versa
prs.ncorrbins = 100; % 100 bins of data in each trial
prs.pretrial = 0.4; % (s)
prs.posttrial = 0.4; % (s)
prs.min_intersaccade = 0.1; % (s) minimum inter-saccade interval
prs.position_thresh = 60; % threshold for maximum predicted eye positions

prs.velthresh = 10; % velocity threshold for creating saccade-free target-tracking index

prs.presac_t = 0.1; % s, pre-saccadic eye position to calculate saccade magnitude and direction
prs.postsac_t = 0.15; % s, post-saccadic eye position to calculate saccade magnitude and direction

% % time window for psth of event aligned responses
% prs.temporal_binwidth = 0.02; % time binwidth for neural data analysis (s)
% prs.spkkrnlwidth = 0.05; % width of the gaussian kernel convolved with spike trains (s)
% prs.spkkrnlwidth = prs.spkkrnlwidth/prs.temporal_binwidth; % width in samples
% prs.spkkrnlsize = round(10*prs.spkkrnlwidth);
% prs.ts.move = -0.5:prs.temporal_binwidth:3.5;
% prs.ts.target = -0.5:prs.temporal_binwidth:3.5;
% prs.ts.stop = -3.5:prs.temporal_binwidth:0.5;
% prs.ts.reward = -3.5:prs.temporal_binwidth:0.5;
% prs.peaktimewindow = [-0.5 0.5]; % time-window around the events within which to look for peak response
% prs.minpeakprominence = 2; % minimum height of peak response relative to closest valley (spk/s)
% 
% % time-rescaling analysis
% prs.ts_shortesttrialgroup.move = -0.5:prs.temporal_binwidth:1.5;
% prs.ts_shortesttrialgroup.target = -0.5:prs.temporal_binwidth:1.5;
% prs.ts_shortesttrialgroup.stop = -1.5:prs.temporal_binwidth:0.5;
% prs.ts_shortesttrialgroup.reward = -1.5:prs.temporal_binwidth:0.5;
% prs.ntrialgroups = 5; % number of groups based on trial duration
% 
% % correlograms
% prs.duration_zeropad = 0.05; % zeros to pad to end of trial before concatenating (s)
% prs.corr_lag = 1; % timescale of correlograms +/-(s)
% 
% % computing standard errors
% prs.nbootstraps = 100; % number of bootstraps for estimating standard errors
% 
% % define no. of bins for tuning curves by binning method
% prs.tuning.nbins1d_binning = 10; % bin edges for tuning curves by 'binning' method
% prs.tuning.nbins2d_binning = [10;10]; % define bin edges for 2-D tuning curves by 'binning' method
% % define no. of nearest neighbors for tuning curves by k-nearest neighbors method
% prs.tuning.k_knn = @(x) round(sqrt(x)); % k=sqrt(N) where N is the total no. of observations
% prs.tuning.nbins1d_knn = 100; prs.tuning.nbins2d_knn = [100 ; 100];
% % define kernel type for tuning curves by Nadayara-Watson kernel regression
% prs.tuning.kernel_nw = 'Gaussian'; % choose from 'Uniform', 'Epanechnikov', 'Biweight', 'Gaussian'
% prs.tuning.bandwidth_nw = []; prs.tuning.bandwidth2d_nw = [];
% prs.tuning.nbins_nw = []; prs.tuning.nbins2d_nw = [];
% % define kernel type for tuning curves by local linear regression
% prs.tuning.kernel_locallinear = 'Gaussian'; % choose from 'Uniform', 'Epanechnikov', 'Biweight', 'Gaussian'
% prs.tuning.bandwidth_locallinear = [];
% 
% % range of stimulus values [min max]
% prs.binrange.v = [0 ; 200]; %cm/s
% prs.binrange.w = [-90 ; 90]; %deg/s
% prs.binrange.r_targ = [0 ; 400]; %cm
% prs.binrange.d = [0 ; 400]; %cm
% prs.binrange.phi = [-90 ; 90]; %deg
% prs.binrange.eye_ver = [-25 ; 5]; %deg
% prs.binrange.eye_hor = [-40 ; 40]; %deg
% prs.binrange.target = [-0.24 ; 0.48];
% prs.binrange.move = [-0.36 ; 0.36];
% prs.binrange.stop = [-0.48 ; 0.24];
% prs.binrange.reward = [-0.36 ; 0.36];
% 
% % fitting models to neural data
% prs.neuralfiltwidth = 10;
% prs.nfolds = 3; % number of folds for cross-validation
% 
% % Generalised additive model for feature tuning
% prs.GAM_beta = 5e1; % hyperparameter to penalise unsparse coupling
% prs.GAM_alpha = 0.05; % significance level for model comparison
% 
% % Gradient descent - parameters
% prs.GD_alpha = 1;
% prs.GD_niters = 200;
% prs.GD_featurescale = false;
% prs.GD_modelname = 'LR'; % name of model to fit: linear regression == 'LR'

% %% hash table to map layman terms to variable names
% prs.varlookup = containers.Map;
% prs.varlookup('target') = 't_targ';
% prs.varlookup('move') = 't_move';
% prs.varlookup('stop') = 't_stop';
% prs.varlookup('reward') = 't_rew';
% 
% %% plotting parameters
% prs.binwidth_abs = prs.temporal_binwidth; % use same width as for the analysis
% prs.binwidth_warp = 0.01;
% prs.trlkrnlwidth = 50; % width of the gaussian kernel for trial averaging (number of trials)
% prs.maxtrls = 5000; % maximum #trials to plot at once.
% prs.rewardwin = 65; % size of reward window (cm)
prs.maxrewardwin = 400; % maximum reward window for ROC analysis
% 
% %% list of analyses to perform
% % *specify methods and variables for analyses (fewer => faster obvisously)*
% % traditional methods
% prs.tuning_events = {'move','target','stop','reward'}; % discrete events - choose from elements of event_vars (above)
% prs.tuning_continuous = {'v','w','d','phi'}; % continuous variables - choose from elements of continuous_vars (above)
% prs.tuning_method = 'binning'; % choose from (increasing computational complexity): 'binning', 'k-nearest', 'nadaraya-watson', 'local-linear'
% % GAM fitting
% prs.GAM_varname = {'v','w','d','phi','r_targ','eye_ver','eye_hor'}; % list of variable names to include in the generalised additive model
% prs.GAM_vartype = {'1D','1D','1D','1D','1D','1D','1D'}; % type of variable: '1d', '1dcirc'
% prs.GAM_linkfunc = 'log'; % choice of link function: 'log','identity','logit'
% prs.GAM_nbins = {10,10,10,10,10,10,10}; % number of bins for each variable
% prs.GAM_lambda = {5e1,5e1,5e1,5e1,5e1,5e1,5e1}; % hyperparameter to penalise rough weight profiles
% prs.GAM_varchoose = [0,0,0,0,0,0,0]; % set to 1 to always include a variable, 0 to make it optional
% % population analysis
% prs.canoncorr_vars = {'v','w','d','phi'}; % list of variables to include in the task variable matrix
% prs.simulate_vars = {'v','w','d','phi'}; % list of variables to use as inputs in simulation
% prs.popreadout_continuous = {'v','w','d','phi','r_targ','alpha','beta'};
% 
% % ****which analyses to do****
% % behavioural
prs.split_trials = true; % split trials into different stimulus conditions
prs.regress_behv = true; % regress response against target position
prs.regress_eye = false; % regress eye position against target position
% % traditional methods
% prs.evaluate_peaks = false; % evaluate significance of event-locked responses
% prs.compute_tuning = false; % compute tuning functions
% % GAM fitting
% prs.fitGAM_tuning = true; % fit generalised additive models to single neuron responses using both task variables + events as predictors
% prs.fitGAM_coupled = false; % fit generalised additive models to single neuron responses with cross-neuronal coupling
% % population analysis
% prs.compute_canoncorr = false; % compute cannonical correlation between population response and task variables
% prs.regress_popreadout = false; % regress population activity against individual task variables
% prs.simulate_population = false; % simulate population activity by running the encoding models
