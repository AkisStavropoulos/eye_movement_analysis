%% Import data
data = ImportSMR('kl015.smr'); % vmax = 50 cm/s, wmax = 90 deg/s, NO FEEDBACK
data = ImportSMR('kl018.smr'); % vmax = 90 cm/s, wmax = 90 deg/s, NO FEEDBACK

data = ImportSMR('kl036.smr'); % vmax = 120 cm/s, wmax = 80 deg/s, FEEDBACK

data = ImportSMR('kl054.smr'); % vmax = 50 cm/s, wmax = 90 deg/s, VESTIBULAR, problematic indicator

data = ImportSMR('kl057.smr'); % vmax = 153 cm/s, wmax = 90 deg/s, ACC.COEF.=0.9
data = ImportSMR('kl065.smr'); % vmax = 100 cm/s, wmax = 90 deg/s, ACC.COEF.=0.9
data = ImportSMR('kl066.smr'); % vmax = 100 cm/s, wmax = 90 deg/s, ACC.COEF.=0.99
% INTERLEAVED TRIALS
data = ImportSMR('kl078.smr'); % vmax = 90 cm/s, wmax = 90 deg/s, ACC.COEF.=0

%
data = ImportSMR('kaushik_10-Jul-2018_18-42_smr.smr');
data = ImportSMR('m51c0843.smr');
%% Import MC variables
data = dlmread('MC_Variables601');
mcvariables = AddMCData(data,prs);
%% Extract data
% extract all data from current folder
default_prs;
trials = AddTrials2Behaviour(prs);
% extract one session of the folder, after importing the file
% data = ImportSMR('kl078.smr'); 
[trials,ch] = AddSMRData(data,prs); % need both outputs for the moment(trl,ch), otherwise only trl
%
dt = prs.dt; % 0.0012
ts = dt:dt:dt*length(ch.tsi);

indxst = find(data(chno.('mrk')).imp.mrk == 115);
indxfb = find(data(chno.('mrk')).imp.mrk == 102);


starts = ch.mrk(indxst);
fb = ch.mrk(indxfb);

marker = zeros(size(ts));
for i = 1:length(fb) % fb or starts
    indx = find(ts(fb(i) <= ts),1,'last');
    marker(indx) = 6;
end

figure; plot(ts,ch.tsi);
hold on;
plot(ts,marker);
plot(ts,0.01*ch.xmp);

%% Plot all channels
chnames = fieldnames(ch);
for i = 1:length(ch_title)
    if ~any(strcmp(chnames{i},'mrk'))
        d = length(ch.(chnames{i})) - length(ts);
        figure;plot(ts,ch.(chnames{i})(1:end-d));title(chnames{i});xlabel('time (s)');vline([t.beg]);
    end
end
%% check gaussian and butter filtering
check_filters;

%% check sensors
check_sensors;

%% check sensor integrals over all trials and correlations
sensor_corr;    
    
%% plot firefly position and trajectory
% better use the indices produced below
for i = 235 %round(length(trials)/2):length(trials) %indx_s(end:-10:end-40)
    figure;plot(trials(i).continuous.xmp,trials(i).continuous.ymp);axis equal;
    hold on;plot(trials(i).prs.fireflyposx,trials(i).prs.fireflyposy,'r*');
    title(['trial No. ' num2str(i)]);
end

% plot joystick input too, to see whether subjects moved
for i = 235 %round(length(trials)/2):length(trials) %indx_s(end:-10:end-40)
    figure;plot(trials(i).mc.timestamp,[trials(i).mc.JS_Yaw_Raw,trials(i).mc.JS_X_Raw]);
%     hold on;plot(trials(i).continuous.ts,[trials(i).continuous.w trials(i).continuous.v]);
    title(['trial No. ' num2str(i)]);
end


%% choose trials based on Joystick Coefficient
indx_a = [];
a = 0.95; % choose 0, .95, .975, .99
for i = 1:length(trials)
    if trials(i).prs.js_coef == a
        indx_a = [indx_a i];
    end
end
%% Choose trials based on Stimulus Type
indx_s = [];
s = 1; % choose 1, 2, 3
for i = 1:length(trials)
    if trials(i).prs.stimtype == s
        indx_s = [indx_s i];
    end
end
%% Choose trials based on Stim Type AND Joystick coefficient 
indx_sa = intersect(indx_s, indx_a);

%% Create struct for condition sorting
% choose trials based on conditions
s = [1 2 3];
a = [0 .95 .975 .99];
% create structs for conditions
stimtype = [];
for i = 1:length(s)
    stimtype = setfield(stimtype,['s' num2str(s(i))],[]);
end
jscoef = [];
for i = 1:length(a)
    co = num2str(a(i));
    if length(co)>1
        co = co(3:end);
    end
    jscoef = setfield(jscoef,['a' co],[]);
end
%
for j = 1:length(trials)
    indx_s = find(s == trials(j).prs.stimtype);
    if indx_s
    stimtype.(['s' num2str(s(indx_s))]) = [stimtype.(['s' num2str(s(indx_s))]) j];
    end
    
    indx_a = find(a == trials(j).prs.js_coef);
    if indx_a
        co = num2str(a(indx_a));
        if length(co)>1
            co = co(3:end);
        end
        jscoef.(['a' co]) = [jscoef.(['a' co]) j];
    end
end
stimtype.all = 1:length(trials);
jscoef.all = 1:length(trials);     
%
stnames = fieldnames(stimtype);
coefnames = fieldnames(jscoef);
%% calculate distance from target at the end of trial

for j = 1:length(stnames)
    d = [];
    for i = 1:length(stimtype.(stnames{j}))
        indx = stimtype.(stnames{j})(i);
        dtemp = sqrt((trials(indx).continuous.xmp(end) - trials(indx).prs.fireflyposx)^2 +...
            (trials(indx).continuous.ymp(end) - trials(indx).prs.fireflyposy)^2);
        d = [d dtemp];
    end
    distance.(stnames{j})  = mean(d);
    
    for k = 1:length(coefnames)-1
        d2 = [];
        indx_sa = intersect(stimtype.(stnames{j}),jscoef.(coefnames{k}));
        for i = 1:length(indx_sa)
            indx = indx_sa(i);
            d2temp = sqrt((trials(indx).continuous.xmp(end) - trials(indx).prs.fireflyposx)^2 +...
                (trials(indx).continuous.ymp(end) - trials(indx).prs.fireflyposy)^2);
            d2 = [d2 d2temp];
        end
        distance.([stnames{j} '_' coefnames{k}]) = mean(d2);
    end
end
% 
% for k = 1:length(coefnames)
%     d3 = [];
%     for i = 1:length(jscoef.(coefnames{k}))
%         indx = jscoef.(coefnames{k})(i);
%         d3temp = sqrt((trials(indx).continuous.xmp(end) - trials(indx).prs.fireflyposx)^2 +...
%             (trials(indx).continuous.ymp(end) - trials(indx).prs.fireflyposy)^2);
%         d3 = [d3 d3temp];
%     end
%     distance.(coefnames{k})  = mean(d3);
% end
% 
%% check targets distribution
for i = 1:length(trials)
    fireflyposx(i) = trials(i).prs.fireflyposx - trials(i).continuous.xmp(1);
    fireflyposy(i) = trials(i).prs.fireflyposy - trials(i).continuous.ymp(1);
end
% FOV
xx = -tand(38)*700:0;
yy = -tand(38)*xx;
y = sqrt(-xx.^2 + 700^2);

figure;plot(fireflyposx,fireflyposy,'.k');axis equal;
hold on;
plot(xx,yy,'r');plot(-xx,yy,'r');
plot(xx,y,'r');plot(-xx,y,'r');
title('Firefly distribution in FOV');ylabel('y (cm)');xlabel('x (cm)');
%% Plot trajectories
% different colors for each stimulus type
condition = [{'vestibular'} {'visual'} {'combined'}];
for i = 1:length(stnames)-1
    colr = ['k' 'r' 'b'];
    figure;
    for j = 1:length(stimtype.(stnames{i}))
        indx = stimtype.(stnames{i})(j);
        plot(trials(indx).continuous.xmp,trials(indx).continuous.ymp,colr(i));hold on;
    end
    title(['trajectories: ' condition{i}]);hold off; axis equal;
end
% same color for all trials
figure;
for i = 1:length(trials)
    plot(trials(i).continuous.xmp,trials(i).continuous.ymp, 'k');hold on;
end
hold off;axis equal;

% different colors for each JS coefficient
condition = [{'a = 0'} {'a = 0.975'} {'a = 0.99'}];
for i = 1:length(coefnames)-1
    colr = ['k' 'r' 'm'];
    figure;
    for j = 1:length(jscoef.(coefnames{i}))
        indx = jscoef.(coefnames{i})(j);
        plot(trials(indx).continuous.xmp,trials(indx).continuous.ymp,colr(i));hold on;
    end
    title(['trajectories: ' condition{i}]);hold off; axis equal;
end

%% Plot INDIVIDUAL trajectories
N = 1:length(trials); % can change it to include the trials you want, or even random
% N = [13 33 47];
for i = N
    figure;plot(trials(i).continuous.xmp,trials(i).continuous.ymp);hold on; 
    plot(trials(i).prs.fireflyposx,trials(i).prs.fireflyposy,'r*');
    title(['trial no. ' num2str(i)]);ylabel('y (cm)');xlabel('x (cm)');
end
    
%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
% polar transformation
theta_sub = [];
r_sub = [];
theta_tar = [];
r_tar = [];
for i = 1:length(trials)
    y_0 = trials(i).continuous.ymp(1);
    y_s = trials(i).continuous.ymp(end);
    y_t = trials(i).prs.fireflyposy;
    x_0 = trials(i).continuous.xmp(1);
    x_s = trials(i).continuous.xmp(end);
    x_t = trials(i).prs.fireflyposx;
    %     if (x_s - x_0) < 0 && (y_s - y_0) < 0
    %         break;
    %     else
    theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
    %     end
    theta_sub = [theta_sub theta_s];
    
    theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
    theta_tar = [theta_tar theta_t];
    
    r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
    r_sub = [r_sub r_s];
    
    r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
    r_tar = [r_tar r_t];
end
%
x = -1000:1000;
y = x;
%
figure;
subplot(1,2,1); plot(r_tar,r_sub,'.b');title('scatterplot of radius');
hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;

subplot(1,2,2); plot(theta_tar,theta_sub,'.k');title('scatterplot of \theta');
hold on;plot(x,y,'r--');plot(x,-y,'r--');
ylim([-180 180]);xlim([-180 180]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;

%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
% different color for each condition
% polar transformation
x = -1000:1000;
y = x;
% figure;
for j = 1:length(stnames)-1
    theta_sub = [];
    r_sub = [];
    theta_tar = [];
    r_tar = [];
    colr = ['k' 'r' 'g'];
    for i = 1:length(stimtype.(stnames{j}))
        indx = stimtype.(stnames{j})(i);
        y_0 = trials(indx).continuous.ymp(1);
        y_s = trials(indx).continuous.ymp(end);
        y_t = trials(indx).prs.fireflyposy;
        x_0 = trials(indx).continuous.xmp(1);
        x_s = trials(indx).continuous.xmp(end);
        x_t = trials(indx).prs.fireflyposx;
        %     if (x_s - x_0) < 0 && (y_s - y_0) < 0
        %         break;
        %     else
        theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
        %     end
        theta_sub = [theta_sub theta_s];
        
        theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
        theta_tar = [theta_tar theta_t];
        
        r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
        r_sub = [r_sub r_s];
        
        r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
        r_tar = [r_tar r_t];
    end
    figure;

    subplot(1,2,1); hold on;plot(r_tar,r_sub,[colr(j) '.']);title('scatterplot of radius');
    hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
    
    subplot(1,2,2); hold on;plot(theta_tar,theta_sub,[colr(j) '.']);title('scatterplot of \theta');
    hold on;plot(x,y,'r--');plot(x,-y,'r--');
    ylim([-180 180]);xlim([-180 180]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
end

%% max and min distance from targets

for i = 1:length(trials)
    dfp(i) = sqrt((trials(i).prs.fireflyposy - trials(i).continuous.ymp(1))^2 + (trials(i).prs.fireflyposx - trials(i).continuous.xmp(1))^2);
end
max_d = max(dfp);
min_d = min(dfp);
indx_d_max = find(dfp==max(dfp));
%% max and min angle from targets

for i = 1:length(trials)
    afp(i) = atan2d((trials(i).prs.fireflyposx - trials(i).continuous.xmp(1)), (trials(i).prs.fireflyposy - trials(i).continuous.ymp(1)));
end
max_a = max(afp);
min_a = min(afp);
indx_a_max = find(afp==max(afp));
%% Find bad trials / check trial length
trial_length = [];
for i = 1:length(trials)
    trial_l = length(trials(i).continuous.ts);
    trial_length = [trial_length trial_l];
end

short_trials = find(trial_length < 43);
%% Find bad trials / check whether subject has not moved
travel = [];
for i = 1:length(trials)
    travel_temp = sqrt((trials(i).continuous.ymp(end)-trials(i).continuous.ymp(1))^2 + (trials(i).continuous.xmp(end)-trials(i).continuous.xmp(1))^2);
    travel = [travel travel_temp];
end

imob_trls = find(travel < 50); % 50 cm
% trials(imob_trls) = [];

%% Remove bad trials
while min(trial_length)<42 % found that max duration of bad trials is 41 (pulse interval)
    trial_length = [];
    for i = 1:length(trials)
        trial_l = length(trials(i).continuous.ts);
        trial_length = [trial_length trial_l];
    end
%     min(trial_length)
    wrongtrl = find(trial_length==min(trial_length));
    trials(wrongtrl) = [];
end

wrtrls = find(trial_length<=41)
%% Check origin of subject position
figure;
origin_x = [];
origin_y = [];
for i = 1:length(trials)
    origin_x = [origin_x trials(i).continuous.xmp(1)];
    origin_y = [origin_y trials(i).continuous.ymp(1)];
    plot(trials(i).continuous.xmp(1),trials(i).continuous.ymp(1),'.');hold on;title('origin position of every trial')
end
xlabel('x');ylabel('y');axis equal;hold off;

%% distance from target at start of trial
dmonk = []; dorigin = [];
for i = 1:length(trials)
    
    dtemp1 = sqrt((trials(i).continuous.xmp(1) - trials(i).prs.fireflyposx)^2 +...
        (trials(i).continuous.ymp(1) - trials(i).prs.fireflyposy)^2);
    dmonk = [dmonk dtemp1];
    
    dtemp2 = sqrt((trials(i).prs.fireflyposx)^2 + (trials(i).prs.fireflyposy)^2);
    dorigin = [dorigin dtemp2];
    
end

yfp = []; xfp = [];
for i = 1:length(trials)
  
    xfptemp = trials(i).prs.fireflyposx;
    xfp = [xfp xfptemp];
  
    yfptemp = trials(i).prs.fireflyposy;
    yfp = [yfp yfptemp];
  
end

xorigin = []; yorigin = [];
for i = 1:length(trials)
    
    xtemp = trials(i).continuous.xmp(1);
    xorigin = [xorigin xtemp];
    
    ytemp = trials(i).continuous.ymp(1);
    yorigin = [yorigin ytemp];
end

xyorigin = sqrt((xorigin).^2 + (yorigin).^2);

[~,dmonk_max_ind] = max(dmonk);
[~,dorigin_max_ind] = max(dorigin);
[~,xorigin_min_ind] = min(xorigin);
[~,yorigin_min_ind] = min(yorigin);
[~,xyorigin_max_ind] = max(xyorigin);



