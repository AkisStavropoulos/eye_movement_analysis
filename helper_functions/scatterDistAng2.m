function [r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng2(trials,stimtype,fit,plt)
%% Scatterplot of distance and angle of target and end of trajectory (arena coordinates)
% fit: yes 1 or no 0
if nargin < 3
fit = 0;
end
% for different Stimulus type conditions (vestibular, visual, combined)
% polar transformation
stnames = fieldnames(stimtype);
x = -1000:1000;
y = x;
condition = [{'vestibular'} {'visual'} {'combined'} {'joystick'}];

for j = 1:length(stnames)-1
    theta_sub1 = [];
    r_sub1 = [];
    theta_tar1 = [];
    r_tar1 = [];
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
        theta_sub1 = [theta_sub1 theta_s];
        
        theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
        theta_tar1 = [theta_tar1 theta_t];
        
        r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
        r_sub1 = [r_sub1 r_s];
        
        r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
        r_tar1 = [r_tar1 r_t];
    end
    % save data
    theta_sub.(stnames{j}) = theta_sub1;
    theta_tar.(stnames{j}) = theta_tar1;
    r_sub.(stnames{j}) = r_sub1;
    r_tar.(stnames{j}) = r_tar1;
    
    % plot data
    if stnames{j} == 's1'
        cond = condition{1};
    elseif stnames{j} == 's2'
        cond = condition{2};
    elseif stnames{j} == 's3'
        cond = condition{3};
    elseif stnames{j} == 's4'
        cond = condition{4};
    end
if plt
    colr = colormap(jet);
    if fit
        y_opt_r = [];
        y_opt_th = [];
 
        [params_r,fitline_r,y_r] = linmodel(r_tar1,r_sub1);
        [params_th,fitline_th,y_th] = linmodel(theta_tar1,theta_sub1);
        y_opt_r = fitline_r(params_r);
        y_opt_th = fitline_th(params_th);
        
        y_r = y_r(x,params_r);
        y_th = y_th(x,params_th);
        
        figure;
        suptitle([trials(i).prs.subject ' - ' cond ' condition']);
        
        subplot(1,2,1); hold on;plot(r_tar.(stnames{j}),r_sub.(stnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');
        title('scatterplot of radius');grid on;
        hold on;plot(x,y,'r--');plot(r_tar1,y_opt_r,'.r');plot(x,y_r,'r');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        
        subplot(1,2,2); hold on;plot(theta_tar.(stnames{j}),theta_sub.(stnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');
        title('scatterplot of \theta');grid on;
        hold on;plot(x,y,'r--');plot(x,-y,'r--');plot(theta_tar1,y_opt_th,'.r');plot(x,y_th,'r');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    else
        figure;
        suptitle([trials(i).prs.subject ' - ' cond ' condition']);
        
        subplot(1,2,1); hold on;plot(r_tar.(stnames{j}),r_sub.(stnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');
        title('scatterplot of radius');grid on;
        hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        
        
        subplot(1,2,2); hold on;plot(theta_tar.(stnames{j}),theta_sub.(stnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');
        title('scatterplot of \theta');grid on;
        hold on;plot(x,y,'r--');plot(x,-y,'r--');
        ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        
    end
    % obtain figures to rearrange them on the screen later
    fig(j) = gcf;
    fig(j).Units = 'centimeters';
end
    
end
if plt
% place figures on Right screen in vertical order
count = 0;
for i = 1:length(fig)
    figdim = fig(i).Position(3:4);
   fig(i).Position = [-1.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
    count = count + 1;    
end
end
