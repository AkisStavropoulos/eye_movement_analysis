function [r_tar,r_sub,theta_tar,theta_sub,b_r,b_th] = plot_opt_scatter(trials,models,condition,fit,plt)
% input trials
% input models: includes optmized parameters and optimized function values
% input conditions: stimtype, jscoef, stimjs

% models.optprs: either 1 x 4 vector, or M x 4 (depends on the M number of initial values tested)
% models.fval: 1 x 1 or M x 1 (same as for optprs)
% this function chooses the optprs that produced the lowest fval.
% For plotting all optprs obtained by all the different initial values
% tested, you'll have to tweak it.
% fit: to fit enter 1, no fit 0
% plt: plot if 1, no plot if 0
%% check inputs
if nargin < 3
    condition = false;
    fit = 0;
    plt = 0;
elseif nargin < 4 || isempty(fit)
    fit = 0;
    plt = 0;
elseif nargin < 5
    plt = 0;
end

if exist('condition')
    condnames = fieldnames(condition);
    for n = 1:length(condnames)-1
        indx{n} = condition.(condnames{n});
    end
else
    indx = {1:length(trials)};
end

clear global;
global x_f y_f x_m0 y_m0 speed;

if exist('condition')
    condnames = fieldnames(condition);
    params = [trials.prs];
    % for different Joystick Coefficient conditions
    coefs = unique([params.js_coef]);
    if length(coefs) > 1
        for i = 1:length(coefs)
            conditionJS(i) = {['a = ' num2str(coefs(i))]};
        end
    else conditionJS  = {['a = ' num2str(coefs)]};
    end
    % for different Stimulus type conditions (vestibular, visual, combined)
    conditionS = [{'vestibular'} {'visual'} {'combined'}  {'joystick'}];
end
% create condition names
if length(condnames) > 3
    for n = 1:length(condnames)-1
        clear conds; clear condjs;
        if any(strcmp(condnames{n}(1:2),'s1'))
            conds = conditionS{1};
        elseif any(strcmp(condnames{n}(1:2),'s2'))
            conds = conditionS{2};
        elseif any(strcmp(condnames{n}(1:2),'s3'))
            conds = conditionS{3};
                elseif any(strcmp(condnames{n}(1:2),'s4'))
            conds = conditionS{4};
        end
        
        if any(strcmp(condnames{n}(3:end),'a0'))
            condjs = conditionJS{1};
        elseif any(strcmp(condnames{n}(3:end),'a99'))
            condjs = conditionJS{3};
        else
            condjs = conditionJS{2};
        end        
        cond{n} = [conds ' and ' condjs];       
    end
    
elseif length(condnames) <= 3
    if any(strcmp(condnames{n}(1:2),'s1'))
        cond = conditionS;
    elseif any(strcmp(condnames{n}(3:end),'a0'))
        cond = conditionJS;
    end
end
% provide initial positions of figures on screen, ONLY FOR ALL 9 CONDITIONS
figposx = [-50 -50 -50 -50 -50 -50 -0.5 -0.5 -0.5];
figposy = [18 18 18 3 3 3 18 18 18];

for n = 1:length(indx)
    clear trls;
    trls = trials(indx{n});
    %% speeds and positions
    
    ntrls = length(trls);
    x_f = NaN(ntrls,1); y_f = NaN(ntrls,1);
    x_m0 = NaN(ntrls,1); y_m0 = NaN(ntrls,1);
    speed = struct('w',cell(ntrls,1),'v',cell(ntrls,1));
    for i=1:ntrls
        x_f(i) = mean(trls(i).prs.fireflyposx);
        y_f(i) = mean(trls(i).prs.fireflyposy);
        x_m0(i) = trls(i).continuous.xmp(1);
        y_m0(i) = trls(i).continuous.ymp(1);
        speed(i).w = trls(i).continuous.w;
        speed(i).v = trls(i).continuous.v;
    end
    %% save and plot data
    % extract trajectory from optimized parameters
    if exist('condition')
        optprs = models.(condnames{n}).optprs;
        fval = models.(condnames{n}).fval;
    else
        optprs = models.optprs;
        fval = models.fval;
    end
    
    % if optprs contains params for different initial points (M x 4), choose the most optimal
    if size(optprs,1) > 1
        [~,ind] = min(fval);
        optprs = optprs(ind,:);
        fval = fval(ind);
    end
    
    [~,sim_traj] = errfun_bayesianmodel(optprs);
    b.sim = sim_traj;
    % create scatterplot for distance and angle
    theta_sub1 = [];r_sub1 = [];theta_tar1 = [];r_tar1 = [];
    for j = 1:length(b.sim)
        y_0 = b.sim(j).y(1);
        y_s = b.sim(j).y(end);
        y_t = y_f(j);
        x_0 = b.sim(j).x(1);
        x_s = b.sim(j).x(end);
        x_t = x_f(j);
        theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
        theta_sub1 = [theta_sub1 theta_s];
        
        theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
        theta_tar1 = [theta_tar1 theta_t];
        
        r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
        r_sub1 = [r_sub1 r_s];
        
        r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
        r_tar1 = [r_tar1 r_t];
    end
        % save data
    theta_sub.(condnames{n}) = theta_sub1;
    theta_tar.(condnames{n}) = theta_tar1;
    r_sub.(condnames{n}) = r_sub1;
    r_tar.(condnames{n}) = r_tar1;
    % plot data
    x = -1000:1000;    y = x;
    if fit
        % distance
        [b_r1, a_r, bint_r, aint_r, r_r, p_r] = regress_perp(r_tar1',r_sub1',0.05,2);
        b_r.b.(condnames{n}) = b_r1;
        b_r.int.(condnames{n}) = bint_r;
        r_fit = b_r1*x;
        if plt
            figure;
            suptitle([cond{n} ', model sim for opt prs [a_\omega b_\omega a_v b_v] = [' num2str(optprs(1),2) '  ' num2str(optprs(2),2) '  ' num2str(optprs(3),2) '  ' num2str(optprs(4),2) ']']);
            subplot(1,2,1); plot(r_tar1,r_sub1,'.b');title('scatterplot of radius');hold on;
            hold on;plot(x,y,'r--');plot(x,r_fit,'r','LineWidth',3);ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
        end
        % theta
        [b_th1, a_th, bint_th, aint_th, r_th, p_th] = regress_perp(theta_tar1',theta_sub1',0.05,2);
        b_th.b.(condnames{n}) = b_th1;
        b_th.int.(condnames{n}) = bint_th;
        th_fit = b_th1*x;
        if plt
            subplot(1,2,2); plot(theta_tar1,theta_sub1,'.k');title('scatterplot of \theta');hold on;
            hold on;plot(x,y,'r--');hline(0);plot(x,th_fit,'r','LineWidth',2);
            ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        end
    else
        if plt
            figure;
            x = -1000:1000;    y = x;
            suptitle([cond{n} ', model sim for opt prs [a_\omega b_\omega a_v b_v] = [' num2str(optprs(1),2) '  ' num2str(optprs(2),2) '  ' num2str(optprs(3),2) '  ' num2str(optprs(4),2) ']']);
            subplot(1,2,1); plot(r_tar1,r_sub1,'.b');title('scatterplot of radius');hold on;
            hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
            
            subplot(1,2,2); plot(theta_tar1,theta_sub1,'.k');title('scatterplot of \theta');hold on;
            hold on;plot(x,y,'r--');hline(0);
            ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
        end
    end
    
    % obtain figures to rearrange them on the screen later
    fig(n) = gcf;
    fig(n).Units = 'centimeters';
    
end
% place figures on Left screen and Top Center screen
% count1 = 0;
% for n = 1:length(fig)
%     figdim = fig(n).Position(3:4);
%     fig(n).Position = [figposx(n)+(count1*(figdim(1)-0.5)) figposy(n) figdim]; % figposx: initial LR position , figposy: initial top-down position
%     count1 = count1 + 1;
%     if ~rem(n,3); count1 = 0; end
% end

