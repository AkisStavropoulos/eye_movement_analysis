function [r_tar,r_sub_but,theta_tar,theta_sub_but] = scatterDistAng_but(trials,but_indx,fit,plt)
%% Scatterplot of distance and angle of target vs FIRST BUTTON PUSH (arena coordinates)
% fit: yes 1 or no 0
if nargin < 2
fit = 0;
end

% polar transformation

theta_sub_but = [];
r_sub_but = [];
theta_tar = [];
r_tar = [];
for i = 1:length(trials)
    y_0 = trials(i).continuous.ymp(1);
    x_0 = trials(i).continuous.xmp(1);
    y_t = trials(i).prs.fireflyposy;
    x_t = trials(i).prs.fireflyposx;
    
    if ~isempty(but_indx{i})
        
    y_s = trials(i).continuous.ymp(but_indx{i});
    x_s = trials(i).continuous.xmp(but_indx{i});

    theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
    r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
    else
        theta_s = nan;
        r_s = nan;
    end
    
    theta_sub_but = [theta_sub_but theta_s];
    r_sub_but = [r_sub_but r_s];

    theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
    theta_tar = [theta_tar theta_t];
    
    r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
    r_tar = [r_tar r_t];
end
%
x = -1000:1000;
y = x;

% plot data
if plt
    if fit
        
        [params_r,fitline_r] = linmodel(r_tar,r_sub_but);
        [params_th,fitline_th] = linmodel(theta_tar,theta_sub_but);
        y_opt_r = fitline_r(params_r);
        y_opt_th = fitline_th(params_th);
        
        figure;         suptitle([trials(i).prs.subject ' - all trials']);grid on;
        subplot(1,2,1); plot(r_tar,r_sub_but,'.b');title('scatterplot of radius');hold on;plot(r_tar,y_opt_r,'.r');
        hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        subplot(1,2,2); plot(theta_tar,theta_sub_but,'.k');title('scatterplot of \theta');hold on;plot(theta_tar,y_opt_th,'.r');
        hold on;plot(x,y,'r--');plot(x,-y,'r--');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    else
        
        figure; suptitle([trials(i).prs.subject ' - all trials']);
        subplot(1,2,1); plot(r_tar,r_sub_but,'.b');title('scatterplot of radius');grid on;
        hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        subplot(1,2,2); plot(theta_tar,theta_sub_but,'.k');title('scatterplot of \theta');
        hold on;plot(x,y,'r--');plot(x,-y,'r--');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    end
end