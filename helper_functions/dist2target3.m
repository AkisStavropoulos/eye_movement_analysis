function [r_tar,err] = dist2target3(trials, tau, plt, fit)
%% Scatterplot of error as a function of initial distance from target
% fit: yes 1 or no 0
% plt: plot 1 or not 0
if nargin < 3
    plt = 0;
    fit = 0;
elseif nargin < 4
    fit = 0;
end
% for different Joystick Coefficient conditions
% polar transformation
taunames = fieldnames(tau);
ind_lim = [];
ind_bin = [];
for i = 1:length(taunames)
    if strcmp(taunames{i}(1),'l')
        ind_lim = [ind_lim i];
    elseif strcmp(taunames{i}(1),'b')
        ind_bin = [ind_bin i];
    end
end
taulimnames = taunames(ind_lim);
taubinnames = taunames(ind_bin);

for i = 1:length(taulimnames)
    condition{i} = {[trials(1).prs.subject ' - ' '\tau : ' '[' num2str(tau.(taulimnames{i})) ']']};
    colr = lines(length(taulimnames));
end

x = -1000:1000;
y = x;


for j = 1:length(taubinnames)
    err1 = [];
    r_tar1 = [];
    for i = 1:length(tau.(taubinnames{j}))
        indx = tau.(taubinnames{j})(i);
        y_0 = trials(indx).continuous.ymp(1); % monkey Y original position y0
        x_0 = trials(indx).continuous.xmp(1); % monkey X original position x0
        y_s = trials(indx).continuous.ymp(end); % monkey Y final position
        x_s = trials(indx).continuous.xmp(end); % monkey X final position
        y_t = trials(indx).prs.fireflyposy; % target position
        x_t = trials(indx).prs.fireflyposx;
        
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
    r_tar.(taubinnames{j}) = r_tar1;
    err.(taubinnames{j}) = err1;
    
    
    % plot data
    if plt
        if fit
            y_opt_err = [];
            
            [params_err,fitline_err,y_err] = linmodel(r_tar1,err1);
            y_opt_err = fitline_err(params_err);
            
            y_err = y_err(x,params_err);
            
            figure;
            suptitle([condition{j}]);
            
            hold on;plot(r_tar.(taubinnames{j}),err.(taubinnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');title('scatterplot of error over distance of targets');
            plot(x,y_err,'r');ylabel('error (cm)');xlabel('r_t (cm)');
            hline(y_err(x == 800));hold off;
            %         hold on;plot(x,y,'r--');plot(r_tar1,y_opt_err,'.r');
            ylim([0 800]);xlim([0 800]);
            
        else
            figure;
            suptitle([condition{j}]);
            
            hold on;plot(r_tar.(taubinnames{j}),err.(taubinnames{j}),'Color',colr(j,:),'Marker','.','LineStyle','none');title('scatterplot of error over distance of targets');
            ylabel('error (cm)');xlabel('r_t (cm)');hold off;
            %         hold on;plot(x,y,'r--');
            ylim([0 800]);xlim([0 800]);
            
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
        fig(i).Position = [-0.5+count*(figdim(1)) 1.5 figdim];  % -50: RL position , 4: top-down position (down-screen)
        count = count + 1;
    end
end
