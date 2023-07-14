%% debug start indicator script
%% fish for faulty sessions out of all included in current folder
flist_smr=dir('*.smr');
for i=1:length(flist_smr)
    fnum_smr(i) = str2num(flist_smr(i).name(end-6:end-4)); 
end
nfiles = length(flist_smr);

flist_log=dir('*.log'); 
for i=1:length(flist_log)
    fnum_log(i) = str2num(flist_log(i).name(end-6:end-4)); 
end

default_prs;
bad_trls_name = [];
countb = 0;
trl_discrep_name = [];
countd = 0;
for i=1:nfiles
    fprintf(['... reading ' flist_smr(i).name '\n']);
    % read .log file
    [trials_log,totaltrls] = AddLOGData(flist_log(i).name);
    % read .smr file
    trials_smr = [];
    data_smr = ImportSMR(flist_smr(i).name);
    [trl] = AddSMRData(data_smr,prs);
    %     trials_smr = [trials_smr AddSMRData(data_smr,prs)];
    % check for bad trials
    trial_length = [];
    for j = 1:length(trl)
        trial_l = length(trl(j).continuous.ts);
        trial_length = [trial_length trial_l];
    end
    
    bad_trials = find(trial_length < 43);
    % PRINT WHICH SESSIONS HAVE BAD TRIALS
    if ~isempty(bad_trials)
        countb = countb+1;
        bad_trls_name{countb} = [flist_smr(i).name];
        fprintf(['bad trials in ' flist_smr(i).name '\n']);
    end
    % check for total trials discrepancy
    if totaltrls ~= length(trl)
        countd = countd+1;
        trl_discrep_name{1,countd} = [flist_log(i).name];   % name
        trl_discrep_name{2,countd} = totaltrls;             % log file trials
        trl_discrep_name{3,countd} = length(trl);           % smr file trials
        fprintf(['trial discrepancy in ' flist_log(i).name '\n']);
     end

end



%% for one session at a time

default_prs;
% extract one session of the folder, after importing the file
file = 'kl081';
data = ImportSMR([file '.smr']); 
[trl,ch,t,dt] = AddSMRData(data,prs); % need both outputs for the moment(trl,ch), otherwise only trl

%% check log/smr trials discrepancy

[logtrls,totaltrls] = AddLOGData([file '.log']);

if totaltrls ~= length(trl)
    fprintf(['trial discrepancy in ' file '\n']);
end

    %% Find bad trials
trial_length = [];
for i = 1:length(trl)
trial_l = length(trl(i).continuous.ts);
trial_length = [trial_length trial_l];
end

bad_trials = find(trial_length < 43);

% 
% if isempty(bad_trials)
%     clear;
% % end
% else
    %% check TSI
    %% check channel headers
    nch = length(data);
    ch_title = cell(1,nch);
    hdr = {data.hdr};
    for i=1:nch
        if ~isempty(hdr{i})
            ch_title{i} = hdr{i}.title;
        else
            ch_title{i} = 'nan';
        end
    end
    ts = dt:dt:dt*length(ch.tsi);
    
    %% Plot tsi and position
    figure; plot(ts,ch.tsi);hold on; plot(ts,.01*ch.xmp);plot(ts,.01*ch.ymp);hold off;
% end
%
%% find other small trials
a = find(trial_length < 1/dt); % trials that lasted less than the flash of the target

% find time of pulse based on the a vector above
b = find(diff(ch.tsi) == 5);
timesofsmalltrials = b(a)*dt;

if (length(trial_length)-totaltrls) == length(a)
    trl(a) = [];
end

%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
% polar transformation
theta_sub = [];
r_sub = [];
theta_tar = [];
r_tar = [];
for i = 1:length(logtrls)
    y_0 = trl(i).continuous.ymp(1);
    y_s = trl(i).continuous.ymp(end);
    y_t = logtrls(i).prs.fireflyposy;
    x_0 = trl(i).continuous.xmp(1);
    x_s = trl(i).continuous.xmp(end);
    x_t = logtrls(i).prs.fireflyposx;
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
suptitle(['file: ' file]);
% clear;