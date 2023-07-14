function [r_tar,r_sub,theta_tar,theta_sub,b_r,b_th] = scatterDistAng4(trials,stimtau,fit,plt)
% scatterplot of radius and theta, subject over target
%% plot different color for each modality given condition
% fit: yes 1 or no 0
if nargin < 3
    fit = 0;
    plt = 0;
elseif nargin < 4
    plt = 0;
end

x = -1000:1000;
y = x;
% for different Tau conditions
condnames1 = fieldnames(stimtau);
ind_keep = [];
for i = 1:length(condnames1)
    if ~strcmp(condnames1{i}(1),'l')
        ind_keep = [ind_keep i];
    end
end
condnames = condnames1(ind_keep);
tau_bins = (length(condnames)-1)/3;

for i = 1:length(condnames)-1
    if strcmp(condnames{i}(1:2),'s1')
        conditionS{i} = 'vestibular';
    elseif strcmp(condnames{i}(1:2),'s2')
        conditionS{i} = 'visual';
    elseif strcmp(condnames{i}(1:2),'s3')
        conditionS{i} = 'combined';
    end
    for j = 1:tau_bins
        if strcmp(condnames{i}(3:end),['bin' num2str(j)])
            conditionJS{i} = ['\tau : ' '[' num2str(stimtau.(['lim' num2str(j)])) ']'];
        end
    end
end

%% create new conditions
count = 0;
for n = 1:length(condnames)-1
        % create condition names
%     clear conds; clear condjs;
    if any(strcmp(condnames{n}(1:2),'s1'))
        conds = conditionS{n};
    elseif any(strcmp(condnames{n}(1:2),'s2'))
        conds = conditionS{n};
    elseif any(strcmp(condnames{n}(1:2),'s3'))
        conds = conditionS{n};
    end
    
    for j = 1:tau_bins
        if strcmp(condnames{n}(3:end),['bin' num2str(j)])
        condjs = conditionJS{j};
        end
    end
   
    condition{n} = [trials(1).prs.subject ' - ' conds ' and ' condjs];
    % create colors for conditions
    colr_temp = lines(tau_bins)';
    colr = repmat(colr_temp,1,length(conditionS));
    colr = colr';
    % provide initial positions of figures on screen
    figposx =  [-1.5 -1.5 -1.5 ...
                -1.5 -1.5 -1.5...
                50 50 50 ...
                50 50 50];
    figposy =  [15 15 15 ...
                1 1 1 ...
                15 15 15 ...
                3 3 3];
        
        theta_sub1 = [];
        r_sub1 = [];
        theta_tar1 = [];
        r_tar1 = [];
        indx = [];
        
        indx = stimtau.(condnames{n});
        for k = 1:length(indx)
            % create scatters
            y_0 = trials(indx(k)).continuous.ymp(1);
            y_s = trials(indx(k)).continuous.ymp(end);
            y_t = trials(indx(k)).prs.fireflyposy;
            x_0 = trials(indx(k)).continuous.xmp(1);
            x_s = trials(indx(k)).continuous.xmp(end);
            x_t = trials(indx(k)).prs.fireflyposx;

            theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
            theta_sub1 = [theta_sub1 theta_s];
            
            theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
            theta_tar1 = [theta_tar1 theta_t];
            
            r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
            r_sub1 = [r_sub1 r_s];
            
            r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
            r_tar1 = [r_tar1 r_t];
        end
        count = count + 1;
        % save data
        theta_sub.(condnames{n}) = theta_sub1;
        theta_tar.(condnames{n}) = theta_tar1;
        r_sub.(condnames{n}) = r_sub1;
        r_tar.(condnames{n}) = r_tar1;
        % plot data
        
        if fit
            % distance
            [b_r1, a_r, bint_r, aint_r, r_r, p_r] = regress_perp(r_tar1',r_sub1',0.05,2);
            b_r.b.(condnames{n}) = b_r1;
            b_r.int.(condnames{n}) = bint_r;
            r_fit = b_r1*x;
            
            if plt
                figure;
                suptitle(condition{n});
                subplot(1,2,1); hold on;grid on;plot(r_tar1,r_sub1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of radius');
                hold on;plot(x,y,'r--');plot(x,r_fit,'r','LineWidth',3);ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
            end
            % theta
            [b_th1, a_th, bint_th, aint_th, r_th, p_th] = regress_perp(theta_tar1',theta_sub1',0.05,2);
            b_th.b.(condnames{n}) = b_th1;
            b_th.int.(condnames{n}) = bint_th;
            th_fit = b_th1*x;
            if plt
                subplot(1,2,2); hold on;grid on;plot(theta_tar1,theta_sub1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of \theta');
                hold on;plot(x,y,'r--');plot(x,-y,'r--');plot(x,th_fit,'r','LineWidth',2);
                ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
            end
        else
            if plt
                figure;
                suptitle(condition{n});
                
                subplot(1,2,1); hold on;grid on;plot(r_tar1,r_sub1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of radius');
                hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
                
                
                subplot(1,2,2); hold on;grid on;plot(theta_tar1,theta_sub1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of \theta');
                hold on;plot(x,y,'r--');plot(x,-y,'r--');
                ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
            end
        end
        % obtain figures to rearrange them on the screen later
        if plt
            fig(n) = gcf;
            fig(n).Units = 'centimeters';
        end
end

% place figures on Left screen and Top Center screen
if plt
    count1 = 0;
    for n = 1:length(fig)
        figdim = fig(n).Position(3:4);
        fig(n).Position = [figposx(n)+(count1*(figdim(1)-2)) figposy(n) figdim]; % figposx: initial LR position , figposy: initial top-down position
%         fig(i).Position = [-1.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
        count1 = count1 + 1;
        if ~rem(n,3); count1 = 0; end
    end
    
end