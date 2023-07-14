%% GIA error distribution
non_visual = sort([stimtype.s1 stimtype.s3]);

GIAerrorX = [];
GIAerrorY = [];
for i = 1:length(non_visual)
GIAerrorX = [GIAerrorX; trials(non_visual(i)).mc.Rat_X_GIAerror];
GIAerrorY = [GIAerrorY; trials(non_visual(i)).mc.Rat_Y_GIAerror];
end


figure;
subplot(2,1,1);histogram(100*GIAerrorX);xlabel('X GIA error (cm/s^2)');ylabel('total timepoints')
xlim([-10 10]);ylim([0 20000]);title('Distribution of GIA error over all trials');

subplot(2,1,2);histogram(100*GIAerrorY);xlabel('Y GIA error (cm/s^2)');ylabel('total timepoints')
xlim([-10 10]);ylim([0 10000]);title('Distribution of GIA error over all trials');

suptitle(trials(1).prs.subject);
%% cumulative GIA error distribution
non_visual = sort([stimtype.s1 stimtype.s3]);
% find longest trial
dt = trials(1).continuous.ts(2) - trials(1).continuous.ts(1);
for i = 1:length(trials)
    trial_dur(i) = length(trials(i).continuous.ts)*dt;
end
% color trials according to trial duration
colr = jet(ceil(max(trial_dur)));

GIAerrXcum = [];
GIAerrYcum = [];
ind_rem = [];
figure;
for i = 1:length(non_visual)
    X_temp = [];
    Y_temp = [];
    X_temp = cumsum(trials(non_visual(i)).mc.Rat_X_GIAerror);
    % plot GIAerror of all trials
    plot((1:length(X_temp))*1/60,X_temp,'Color',colr(ceil(trial_dur(non_visual(i))),:));hold on;
    %
    Y_temp = cumsum(trials(non_visual(i)).mc.Rat_Y_GIAerror);
    if ~isempty(X_temp)
        GIAerrXcum = [GIAerrXcum; X_temp(end)];
        GIAerrYcum = [GIAerrYcum; Y_temp(end)];
    else 
        ind_rem = [ind_rem i];
    end
end
title('Cumulative GIA error for each trial');xlabel('time (s)');ylabel('GIA error (m/s^2)');
suptitle(trials(1).prs.subject);

figure;
subplot(2,1,1);histogram(100*GIAerrXcum,50);xlabel('X GIA error (cm/s^2)');ylabel('total trials')
title('Distribution of cumulative GIA error over all trials');

subplot(2,1,2);histogram(100*GIAerrYcum,50);xlabel('Y GIA error (cm/s^2)');ylabel('total trials')
title('Distribution of cumulative GIA error over all trials');

suptitle(trials(1).prs.subject);
%% cumulative GIA error as a function of distance
non_visual(ind_rem) = [];
travel_d = [];
for i = 1:length(trials)
    dist = sqrt((trials(i).continuous.xmp(end) - trials(i).continuous.xmp(1)).^2 + (trials(i).continuous.ymp(end) - trials(i).continuous.ymp(1)).^2);
    travel_d = [travel_d dist];
end
target_d = target_dist(trials(non_visual));
travel_d = travel_d(non_visual);


x = target_d';
X = [ones(length(x),1) x];
y = GIAerrXcum;
b = regress(y,X);

figure;plot(target_d,100*GIAerrXcum,'.');hold on;
plot(x,100*(x*b(2) + b(1)));title('Cumulative trial GIA error as a function of target distance');
ylabel('GIA error (cm/s^2)');xlabel('target distance (m)');grid on;
suptitle(trials(1).prs.subject);


x = travel_d';
X = [ones(length(x),1) x];
y = GIAerrXcum;
b = regress(y,X);


figure;plot(travel_d,100*GIAerrXcum,'.');hold on;plot(x,100*(x*b(2) + b(1)))
title('Cumulative trial GIA error as a function of distance traveled');
ylabel('GIA error (cm/s^2)');xlabel('distance traveled (m)');grid on;
suptitle(trials(1).prs.subject);
%% GIA error over time
T = trials(1).prs.T + 1;
dt = 1/60;
sf = 1/dt;
GIAerrX = [];
GIAerrY = [];
travel_d = [];
for i = 1:length(trials)
    GIAerrX_temp = [];
    GIAerrY_temp = [];
    GIAerrX_temp = trials(i).mc.Rat_X_GIAerror;
    GIAerrY_temp = trials(i).mc.Rat_Y_GIAerror;
    if ~isempty(GIAerrX_temp)
        GIAerrX(i,:) = interp1(1:length(GIAerrX_temp),GIAerrX_temp',1:(floor(T/dt)));
        GIAerrY(i,:) = interp1(1:length(GIAerrY_temp),GIAerrY_temp',1:(floor(T/dt)));
    end
    dist = sqrt((trials(i).continuous.xmp(end) - trials(i).continuous.xmp(1)).^2 + (trials(i).continuous.ymp(end) - trials(i).continuous.ymp(1)).^2);
    travel_d = [travel_d dist];
end
GIAerrX(isnan(GIAerrX)) = 0;
GIAerrY(isnan(GIAerrY)) = 0;


meanGIAerrX = mean(GIAerrX,2);
meanGIAerrY = mean(GIAerrY,2);
non_visual = sort([stimtype.s1 stimtype.s3]);
non_visual(ind_rem) = [];
target_d = target_dist(trials(non_visual));
travel_d = travel_d(non_visual);


x = target_d';
X = [ones(length(x),1) x];
y = meanGIAerrX(non_visual);
b = regress(y,X);

figure;plot(target_d,100*meanGIAerrX(non_visual),'.');hold on;
plot(x,100*(x*b(2) + b(1)));title('mean trial GIA error as a funciton of target distance');
ylabel('GIA error (cm/s^2)');xlabel('target distance (m)');grid on;

x = travel_d';
X = [ones(length(x),1) x];
y = meanGIAerrX(non_visual);
b = regress(y,X);


figure;plot(travel_d,100*meanGIAerrX(non_visual),'.');hold on;
plot(x,100*(x*b(2) + b(1)));title('mean trial GIA error as a funciton of distance traveled');
ylabel('GIA error (cm/s^2)');xlabel('distance traveled (m)');grid on;

%% cumulative GIA error distribution
non_visual = sort([stimtype.s1 stimtype.s3]);
% find target distance and traveled distance for each trial
target_d = target_dist(trials);
travel_d = traveled_dist(trials);

% color trials according to trial duration
colr = jet(ceil(max(travel_d) - min(travel_d) + 1));

figure;
for i = 1:length(non_visual)
    X_temp = [];
    Y_temp = [];
    X_temp = cumsum(trials(non_visual(i)).mc.Rat_X_GIAerror);
    % plot GIAerror of all trials
    plot((1:length(X_temp))*1/60,X_temp,'Color',colr(ceil(travel_d(non_visual(i)) - min(travel_d) + 1),:));hold on;
    %
    Y_temp = cumsum(trials(non_visual(i)).mc.Rat_Y_GIAerror);
end
title('Cumulative GIA error, colored for TRAVELED distance');xlabel('time (s)');ylabel('GIA error (m/s^2)');
suptitle(trials(1).prs.subject);
