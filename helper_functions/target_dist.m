function dist = target_dist(trials,trialnum)
% Calculate initial distance between subject and target

% trialnum: indices of trials for which you want to find the distance
% can be indices of conditions (stimtype, jscoef, stimjs)

if nargin < 2
    trialnum = 1:length(trials);
end

theta_sub = [];
r_sub = [];
theta_tar = [];
r_tar = [];

for i = 1:length(trialnum)
    y_0 = trials(trialnum(i)).continuous.ymp(1);
    y_s = trials(trialnum(i)).continuous.ymp(end);
    y_t = trials(trialnum(i)).prs.fireflyposy;
    x_0 = trials(trialnum(i)).continuous.xmp(1);
    x_s = trials(trialnum(i)).continuous.xmp(end);
    x_t = trials(trialnum(i)).prs.fireflyposx;
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

dist = r_tar;