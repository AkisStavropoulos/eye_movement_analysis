function [r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(trials,button,fit,plt)
%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
% fit: yes 1 or no 0
% button: 1 to get position at first button push, 0 to get position at end of trial

if nargin < 3
fit = 0;
plt = 0;
end

% polar transformation

theta_sub = [];
r_sub = [];
theta_tar = [];
r_tar = [];
for i = 1:length(trials)
    if button
        
        butpush = find(trials(i).mc.flag2,1);
        
        if ~isnan(trials(i).mc.flag2(butpush))
            
            indx = find(trials(i).continuous.ts >= trials(i).mc.timestamp(butpush),1);

            y_0 = trials(i).continuous.ymp(1);
            y_s = trials(i).continuous.ymp(indx);
            y_t = trials(i).prs.fireflyposy;
            x_s = trials(i).continuous.xmp(indx);
            x_0 = trials(i).continuous.xmp(1);
            x_t = trials(i).prs.fireflyposx;
            
            theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
            
            theta_sub = [theta_sub theta_s];
            
            theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
            theta_tar = [theta_tar theta_t];
            
            r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
            r_sub = [r_sub r_s];
            
            r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
            r_tar = [r_tar r_t];
            
        else
            y_0 = trials(i).continuous.ymp(1);
            y_s = nan;
            y_t = trials(i).prs.fireflyposy;
            x_s = nan;
            x_0 = trials(i).continuous.xmp(1);
            x_t = trials(i).prs.fireflyposx;
            
            theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
            
            theta_sub = [theta_sub theta_s];
            
            theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
            theta_tar = [theta_tar theta_t];
            
            r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
            r_sub = [r_sub r_s];
            
            r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
            r_tar = [r_tar r_t];
            
        end
        
    else
        y_0 = trials(i).continuous.ymp(1);
        y_s = trials(i).continuous.ymp(end);
        y_t = trials(i).prs.fireflyposy;
        x_s = trials(i).continuous.xmp(end);
        x_0 = trials(i).continuous.xmp(1);
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
end
%
x = -1000:1000;
y = x;

% plot data
if plt
    if fit
        
        [params_r,fitline_r] = linmodel(r_tar,r_sub);
        [params_th,fitline_th] = linmodel(theta_tar,theta_sub);
        y_opt_r = fitline_r(params_r);
        y_opt_th = fitline_th(params_th);
        
        figure;         suptitle([trials(i).prs.subject ' - all trials']);grid on;
        subplot(1,2,1); plot(r_tar,r_sub,'.b');title('scatterplot of radius');hold on;plot(r_tar,y_opt_r,'.r');
        hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        subplot(1,2,2); plot(theta_tar,theta_sub,'.k');title('scatterplot of \theta');hold on;plot(theta_tar,y_opt_th,'.r');
        hold on;plot(x,y,'r--');plot(x,-y,'r--');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    else
        
        figure; suptitle([trials(i).prs.subject ' - all trials']);
        subplot(1,2,1); plot(r_tar,r_sub,'.b');title('scatterplot of radius');grid on;
        hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        subplot(1,2,2); plot(theta_tar,theta_sub,'.k');title('scatterplot of \theta');
        hold on;plot(x,y,'r--');plot(x,-y,'r--');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    end
end