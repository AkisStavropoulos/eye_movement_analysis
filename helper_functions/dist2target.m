function [r_tar,err] = dist2target(trials, plt, fit)
%% Scatterplot of error as a function of initial distance from target
% fit: yes 1 or no 0
% plt: plot 1 or not 0
if nargin < 2
    plt = 0;
    fit = 0;
elseif nargin < 3
    fit = 0;
end
% for all trials
% polar transformation
x = -1000:1000;
y = x;

err1 = [];
r_tar1 = [];
for i = 1:length(trials)
    y_0 = trials(i).continuous.ymp(1); % monkey Y original position y0
    x_0 = trials(i).continuous.xmp(1); % monkey X original position x0
    y_s = trials(i).continuous.ymp(end); % monkey Y final position
    x_s = trials(i).continuous.xmp(end); % monkey X final position
    y_t = trials(i).prs.fireflyposy; % target position
    x_t = trials(i).prs.fireflyposx;
    
    %         theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
    %         theta_sub1 = [theta_sub1 theta_s];
    %         theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
    %         theta_tar1 = [theta_tar1 theta_t];
    %         r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
    %         r_sub1 = [r_sub1 r_s];
    
    r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
    r_tar1 = [r_tar1 r_t];
    
    err_temp = sqrt((x_t - x_s)^2 + (y_t - y_s)^2);
    err1 = [err1 err_temp];
    
end
% save data
r_tar = r_tar1;
err = err1;


% plot data
if plt
    if fit
        y_opt_err = [];
        
        [params_err,fitline_err,y_err] = linmodel(r_tar1,err1);
        y_opt_err = fitline_err(params_err);
        
        y_err = y_err(x,params_err);
        
        figure;         suptitle([trials(i).prs.subject ' - all trials']);
        
        hold on;plot(r_tar,err,'b.');
        title('scatterplot of error over distance of targets');
        plot(x,y_err,'r');ylabel('error (cm)');xlabel('r_t (cm)');
        hline(y_err(x == 800));hold off;
        %         hold on;plot(x,y,'r--');plot(r_tar1,y_opt_err,'.r');
        ylim([0 800]);xlim([0 800]);
        
    else
        figure;         suptitle([trials(i).prs.subject ' - all trials']);
        
        hold on;plot(r_tar,err,'k.');
        title('scatterplot of ERROR over distance of targets');
        ylabel('error (cm)');xlabel('r_t (cm)');hold off;
        %         hold on;plot(x,y,'r--');
        ylim([0 800]);xlim([0 800]);
        
    end
end
