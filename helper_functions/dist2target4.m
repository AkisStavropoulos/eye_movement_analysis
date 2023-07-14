function [r_tar,err] = dist2target4(trials, stimtau, plt, fit)
%% Scatterplot of error as a function of initial distance from target
% fit: yes 1 or no 0
% plt: plot 1 or not 0
if nargin < 3
    plt = 0;
    fit = 0;
elseif nargin < 4
    fit = 0;
end
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
%%
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
            
    err1 = [];
    r_tar1 = [];
    indx = [];
    
    indx = stimtau.(condnames{n});
    for k = 1:length(indx)
        y_0 = trials(indx(k)).continuous.ymp(1); % monkey Y original position y0
        x_0 = trials(indx(k)).continuous.xmp(1); % monkey X original position x0
        y_s = trials(indx(k)).continuous.ymp(end); % monkey Y final position
        x_s = trials(indx(k)).continuous.xmp(end); % monkey X final position
        y_t = trials(indx(k)).prs.fireflyposy; % target position
        x_t = trials(indx(k)).prs.fireflyposx;
        
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
    count = count + 1;
    % save data
    r_tar.(condnames{n}) = r_tar1;
    err.(condnames{n}) = err1;
    
    
    % plot data
    if plt
        if fit
            clear x;clear y;
            x = r_tar.(condnames{n})';
            y = err.(condnames{n})';
            [b_r, a_r, bint_r, aint_r, r_r, p_r] = regress_perp(x,y,0.05,2);
            
            x = -1000:1000;
            y = b_r*x;
            figure;
            suptitle(condition{n});
            
            hold on;plot(r_tar1,err1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of error over distance of targets');
            plot(x,y,'r');ylabel('error (cm)');xlabel('r_t (cm)');
            hline(y(x == 800));hold off;
            %         hold on;plot(x,y,'r--');plot(r_tar1,y_opt_err,'.r');
            ylim([0 800]);xlim([0 800]);
            
        else
            figure;
            suptitle(condition{n});
            
            hold on;plot(r_tar1,err1,'Color',colr(count,:),'Marker','.','LineStyle','none');title('scatterplot of error over distance of targets');
            ylabel('error (cm)');xlabel('r_t (cm)');hold off;
            %         hold on;plot(x,y,'r--');
            ylim([0 800]);xlim([0 800]);
            
        end
        % obtain figures to rearrange them on the screen later
        fig(n) = gcf;
        fig(n).Units = 'centimeters';
    end
    
end
% place figures on Right screen in vertical order
if plt
    count1 = 0;
    for j = 1:length(fig)
        figdim = fig(j).Position(3:4);
        fig(j).Position = [figposx(j)+(count1*(figdim(1)-0.5)) figposy(j) figdim]; % figposx: initial LR position , figposy: initial top-down position
        count1 = count1 + 1;
        if ~rem(j,3); count1 = 0; end
    end
end
