function [trl,ch] = AddSMRData(data,prs)
% this version filters the raw data before dividing into trials
%% check channel headers (done)
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

%% channel titles (done)
chno.mrk = find(strcmp(ch_title,'Keyboard')); % keyboard as the marker
chno.tsi = find(strcmp(ch_title,'TrialSta')); % trial start indicator as the marker

chno.yle = find(strcmp(ch_title,'LDy')); 
chno.zle = find(strcmp(ch_title,'LDz'));
chno.yre = find(strcmp(ch_title,'RDy')); 
chno.zre = find(strcmp(ch_title,'RDz'));

chno.phi = find(strcmp(ch_title,'MonkeyYa')); % monkey yaw position
chno.xmp = find(strcmp(ch_title,'MonkeyX')); 
chno.ymp = find(strcmp(ch_title,'MonkeyY'));
chno.v = find(strcmp(ch_title,'ForwardV')); 
chno.w = find(strcmp(ch_title,'AngularV'));
chno.mtr = find(strcmp(ch_title,'MotorCom')); % yaw motor command
chno.xac = find(strcmp(ch_title,'AccX')); % acceleration x-axis
chno.yac = find(strcmp(ch_title,'AccY')); % acceleration y-axis

chno.vrol = find(strcmp(ch_title,'VelRoll')); % roll velocity 
chno.vyaw = find(strcmp(ch_title,'VelYaw')); % yaw velocity 
chno.vpit = find(strcmp(ch_title,'VelPitch')); % pitch velocity 
% chno.key = find(strcmp(ch_title,'Keyboard')); % keyboard , check what it is about


%% scale (done)
scaling.t = data(chno.mrk).hdr.tim.Scale*data(chno.mrk).hdr.tim.Units; % for markers
scaling.tsi = data(chno.tsi).hdr.adc.Scale; offset.tsi = data(chno.tsi).hdr.adc.DC;

scaling.yle = data(chno.yle).hdr.adc.Scale; offset.yle = data(chno.yle).hdr.adc.DC;
scaling.yre = data(chno.yre).hdr.adc.Scale; offset.yre = data(chno.yre).hdr.adc.DC; 
scaling.zle = data(chno.zle).hdr.adc.Scale; offset.zle = data(chno.zle).hdr.adc.DC; 
scaling.zre = data(chno.zre).hdr.adc.Scale; offset.zre = data(chno.zre).hdr.adc.DC;

scaling.phi = data(chno.phi).hdr.adc.Scale; offset.phi = data(chno.phi).hdr.adc.DC;
scaling.xmp = data(chno.xmp).hdr.adc.Scale; offset.xmp = data(chno.xmp).hdr.adc.DC;
scaling.ymp = -data(chno.ymp).hdr.adc.Scale; offset.ymp = -data(chno.ymp).hdr.adc.DC; % !!!!!SIGN
scaling.v = data(chno.v).hdr.adc.Scale; offset.v = data(chno.v).hdr.adc.DC;
scaling.w = data(chno.w).hdr.adc.Scale; offset.w = data(chno.w).hdr.adc.DC;
scaling.mtr = data(chno.mtr).hdr.adc.Scale; offset.mtr = data(chno.mtr).hdr.adc.DC;
scaling.xac = data(chno.xac).hdr.adc.Scale; offset.xac = data(chno.xac).hdr.adc.DC;
scaling.yac = data(chno.yac).hdr.adc.Scale; offset.yac = data(chno.yac).hdr.adc.DC;

scaling.vrol = data(chno.vrol).hdr.adc.Scale; offset.vrol = data(chno.vrol).hdr.adc.DC;
scaling.vyaw = data(chno.vyaw).hdr.adc.Scale; offset.vyaw = data(chno.vyaw).hdr.adc.DC;
scaling.vpit = data(chno.vpit).hdr.adc.Scale; offset.vpit = data(chno.vpit).hdr.adc.DC;

%% load relevant channels (done)
chnames = fieldnames(chno); MAX_LENGTH = inf; dt = [];
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},'mrk'))
        ch.(chnames{i}) = double(data(chno.(chnames{i})).imp.adc)*scaling.(chnames{i}) + offset.(chnames{i});
        dt = [dt prod(data(chno.(chnames{i})).hdr.adc.SampleInterval)];
        MAX_LENGTH = min(length(ch.(chnames{i})),MAX_LENGTH);
    end
end

ch.('mrk') = double(data(chno.('mrk')).imp.tim)*scaling.t;

if length(unique(dt))==1
    dt = dt(1);
else
   error('channels must all have identical sampling rates');
end
%% event markers (done)
markers = data(chno.mrk).imp.mrk(:,1); % from 'keyboard'
%% event times (done)
ts = dt:dt:dt*MAX_LENGTH;
tsi = ch.tsi;
tsi(tsi>=1)= 5;
tsi(tsi<=0.5)= 0;
pulse_max = max(tsi); 
if find(markers == 104) % if FEEDBACK
    indx_beg = find(diff(tsi) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
    indx_end = find(diff(tsi) == pulse_max); % choose the start of the pulse, since it is the command to reset position
    t.beg = ts(indx_beg(1:2:end));
    t.end = ts(indx_end(1:2:end));
    t.end(1) = [];
    t.feedback = ts(indx_beg(2:2:end));
    t.beg = t.beg(1:length(t.end));
    t.events = sort([ts(indx_beg) ts(indx_end)]);
else 
    indx_beg = find(diff(tsi) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
    indx_end = find(diff(tsi) == pulse_max); % choose the start of the pulse, since it is the command to reset position
  
    t.events = sort([ts(indx_beg) ts(indx_end)]);
    t.events(1) = [];     t.events(end) = []; % see tsi plot to check events, 1st event start, last event end
   
    t.beg = t.events(1:2:end);
    t.end = t.events(2:2:end);
    t.beg = t.beg(1:length(t.end));
end
ch.tsi = tsi;
% indxon = find(round(diff(ch.tsi)) == pulse_max);
% indxoff = find(round(diff(ch.tsi)) == -pulse_max);
% pulse_duration  = indxoff - indxon;

%% define gaussian filter (done)
sig = prs.filtwidth; %filter width
sz = prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
%% define Butterworth filter
SR = 1/dt;
[b,a] = butter(2,30/(SR/2),'low'); % for the sensors
%% filter position, velocity and sensor channels (done)
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},{'tsi','mrk','yle','yre','zle','zre'}))
        if any(strcmp(chnames{i},{'xac','yac','vrol','vyaw','vpit'}))
            ch.(chnames{i}) = filtfilt(b,a,ch.(chnames{i})); % butterowrth filter for the sensors
        else
            ch.(chnames{i}) = conv(ch.(chnames{i})(1:MAX_LENGTH),h,'same'); % gaussian filter for data
            %         ch.(chnames{i}) = ch.(chnames{i})(sz/2+1:end);
        end
    end
end
ch.yle = ch.yle(1:MAX_LENGTH);
ch.yre = ch.yre(1:MAX_LENGTH);
ch.zle = ch.zle(1:MAX_LENGTH);
ch.zre = ch.zre(1:MAX_LENGTH);
%% detect saccade times
% take derivative of eye position = eye velocity
if (var(ch.zle) > var(ch.zre)) % use the eye with a working eye coil
    dze = diff(ch.zle);
    dye = diff(ch.yle);
else
    dze = diff(ch.zre);
    dye = diff(ch.yre);
end
de = sqrt(dze.^2 + dye.^2); % speed of eye movement

% apply threshold on eye speed
saccade_thresh = prs.saccade_thresh;
thresh = saccade_thresh/prs.fs_smr; % threshold in units of deg/sample
indx_thresh = de>thresh;
dindx_thresh = diff(indx_thresh);
t_saccade = find(dindx_thresh>0)/prs.fs_smr;

% remove duplicates by applying a saccade refractory period
min_isi = prs.min_intersaccade;
t_saccade(diff(t_saccade)<min_isi) = [];
t.saccade = t_saccade;

%% replace the broken eye coil (if any) with NaNs
if var(ch.zle) < 10 || var(ch.zle) > 1000
    ch.zle(:) = nan;
    ch.yle(:) = nan;
end
if var(ch.zre) < 10 || var(ch.zre) > 1000
    ch.zre(:) = nan;
    ch.yre(:) = nan;
end

%% detect start-of-movement and end-of-movement times for each trial
v_thresh = prs.v_thresh;
v_time2thresh = prs.v_time2thresh;
v = ch.v;
for j=1:length(t.end)
   % start-of-movement
   if j==1, t.move(j) = t.beg(j); % first trial is special because there is no pre-trial period
   else
       indx = find(v(ts>t.end(j-1) & ts<t.end(j)) > v_thresh,1); % first upward threshold-crossing
       if ~isempty(indx), t.move(j) = t.end(j-1) + indx*dt;
       else, t.move(j) = t.beg(j); end % if monkey never moved, set movement onset = target onset
   end
   % end-of-movement
   indx = find(v(ts>t.move(j) & ts<t.end(j)) < v_thresh,1); % first downward threshold-crossing
   if ~isempty(indx), t.stop(j) = t.move(j) + indx*dt;
   else, t.stop(j) = t.end(j); end % if monkey never stopped, set movement end = trial end
   % if monkey stopped prematurely, set movement end = trial end
   if (t.stop(j)<t.beg(j) || t.stop(j)-t.move(j)<0.5), t.stop(j) = t.end(j); end
end

%% extract trials (and downsample for storage)
dt = dt*prs.factor_downsample;
for j=1:length(t.end)
    % define pretrial period
%     pretrial = prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
%     posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
    for i=1:length(chnames)
        if ~any(strcmp(chnames{i},'mrk'))
            trl(j).continuous.(chnames{i}) = ch.(chnames{i})(ts>t.beg(j) & ts<t.end(j));
%             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
        end
    end
    trl(j).continuous.ts = (dt:dt:length(trl(j).continuous.(chnames{2}))*dt)';
    trl(j).continuous.firefly = trl(j).continuous.ts>=0 & trl(j).continuous.ts<(0+prs.fly_ONduration);
    trl(j).events.t_beg = t.beg(j);
    trl(j).events.t_end = t.end(j);
%     trl(j).events.t_move = t.move(j);
%     trl(j).events.t_stop = t.stop(j);
    % saccade time
    trl(j).events.t_sac = t.saccade(t.saccade>(t.beg(j)) & t.saccade<t.end(j));
    % reward time
%     if any(t.reward>t.beg(j) & t.reward<t.end(j))
%         trl(j).logical.reward = true;
%         trl(j).events.t_rew = t.reward(t.reward>t.beg(j) & t.reward<t.end(j));
%     else
%         trl(j).logical.reward = false;
%         trl(j).events.t_rew = nan;
%     end
     % ptb time
%     if any(t.ptb>t.beg(j) & t.ptb<t.end(j))
%         trl(j).logical.ptb = true;
%         trl(j).events.t_ptb = t.ptb(t.ptb>t.beg(j) & t.ptb<t.end(j));
%     else
%         trl(j).logical.ptb = false;
%         trl(j).events.t_ptb = nan;
%     end
end

%% set position values prior to target onset to nan
for j=1:length(trl)
    for i=1:length(chnames)
        if any(strcmp(chnames{i},{'xmp','ymp'}))
            trl(j).continuous.(chnames{i})(trl(j).continuous.ts<0) = nan; % target onset happens exactly at t_beg ?????
        end
    end
end

%% timestamps referenced relative to exp_beg
% exp_beg = t.events(1);
% exp_end = t.events(end);
% 
% for i=1:length(trl)
%     trl(i).events.t_beg = trl(i).events.t_beg - exp_beg;
%     trl(i).events.t_end = trl(i).events.t_end - exp_beg - trl(i).events.t_beg;    
%     trl(i).events.t_sac = trl(i).events.t_sac - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_move = trl(i).events.t_move - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_stop = trl(i).events.t_stop - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_ptb = trl(i).events.t_ptb - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_targ = 0;
% end

%% downsample continuous data
% for i=1:length(chnames)
%     if ~any(strcmp(chnames{i},'mrk'))
%         ch.(chnames{i}) = ch.(chnames{i})(ts>exp_beg & ts<exp_end);
%         ch.(chnames{i}) = downsample(ch.(chnames{i}),prs.factor_downsample);
%     end
% end
% ts = ts(ts>exp_beg & ts<exp_end) - exp_beg;
% ch.ts = downsample(ts,prs.factor_downsample); ch.ts = ch.ts(:);
ch.ntrls = length(trl);
