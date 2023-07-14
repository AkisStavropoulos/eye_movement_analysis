function [opt_coefs,opt_fval,models] = fit_slow_prior(trials,condition,initiation,plt)
%% Test the slow-speed prior model
% models = []
% condition: input condition (stimtype, jscoef, stimjs) or skip (fit all trials)
% initiation: 'explore' starts from different initial points, 'fixed' uses fixed values
% plt: plot if 1, not if 0
%% check inputs
tic;
if nargin == 3
    plt = 0;
elseif nargin == 2
    plt = 0;
    initiation = 'explore';
    disp('No initiation input. Will explore different initial points. ');
elseif nargin == 1
    plt = 0;
    initiation = 'explore';
    disp('No initiation input. Will explore different initial points. ');
    condition = false;
end

if strcmp(initiation,'explore')
    b_w0 = logspace(-1,2,10);
    b_v0 = logspace(-1,2,10);
elseif strcmp(initiation,'fixed')
    b_w0 = 46;
    b_v0 = 46; % these initial points have been found to produce the most optimal values (46.4159)
else
    error('initiation input not recognized. Use ''explore'' or ''fixed''.');
end

if exist('condition')
    condnames = fieldnames(condition);
    for n = 1:length(condnames)-1
        indx{n} = condition.(condnames{n});
    end
else
    indx = {1:length(trials)};
end

modelname = 'prior';

clear global;
global x_f y_f x_m0 y_m0 speed;
opt_coefs = [];
opt_fval = [];
%% iterate through conditions and trials
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
    %% optimize
    count = 1;
    k_w0 = 1; k_v0 = 1; % in the paper, k is beta
    prs0 = [k_w0 k_v0];
    LB = [-2 -2]; UB = [2 2];
    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
    optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
    k_w = optprs(1);
    k_v = optprs(2);
    
    % [xt,yt] = gen_traj(speed(i).w,speed(i).v,x_m0(i),y_m0(i));
    % plot(xt,yt);
    %% optimize for many initial points, to choose the best ones
    for k = 1:length(b_w0)
        for i = 1:length(b_v0)
            model(count).name = modelname;
            if exist('condition')
                disp(['Condition ' condnames{n}]);
            end
            disp(['Iterating through k = ' num2str(k) ' and i = ' num2str(i)]);
            
            xt = [];yt = [];
            a_w0(i) = (k_w-1)/b_w0(k);
            a_v0(i) = (k_v-1)/b_v0(i);
            prs0 = [a_w0(i) b_w0(k) a_v0(i) b_v0(i)];
            LB = [-10 0.05 -10 0.05]; UB = [10 inf 10 inf];
            opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
            [model(count).optprs,model(count).fval,model(count).exitflag,...
                model(count).output,model(count).lambda] = fmincon(@errfun_bayesianmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
            
            % save data
            if exist('condition')
                models.(condnames{n}) = model;
                opt_coefs.(condnames{n})(count,:) = [model(count).optprs]; % [a_w  b_w  a_v  b_v  ;  ... ]
                opt_fval.(condnames{n})(count,:) = [model(count).fval];
            else
                models = model;
                opt_coefs = [opt_coefs ; model(count).optprs]; % [a_w  b_w  a_v  b_v  ;  ... ]
                opt_fval = [opt_fval ; model(count).fval];
                
            end
            % plot data
            if plt
                % extract trajectory from optimized parameters
                [~,sim_traj] = errfun_bayesianmodel(model(count).optprs);
                b(i).sim = sim_traj;
                % create scatterplot for distance and angle
                theta_sub = [];r_sub = [];theta_tar = [];r_tar = [];
                for j = 1:length(b(i).sim)
                    y_0 = b(i).sim(j).y(1);
                    y_s = b(i).sim(j).y(end);
                    y_t = y_f(j);
                    x_0 = b(i).sim(j).x(1);
                    x_s = b(i).sim(j).x(end);
                    x_t = x_f(j);
                    theta_s = atan2d((x_s - x_0),(y_s - y_0)); % inverted x and y
                    theta_sub = [theta_sub theta_s];
                    
                    theta_t = atan2d((x_t - x_0),(y_t - y_0)); % inverted x and y
                    theta_tar = [theta_tar theta_t];
                    
                    r_s = sqrt((x_s - x_0)^2 + (y_s - y_0)^2);
                    r_sub = [r_sub r_s];
                    
                    r_t = sqrt((x_t - x_0)^2 + (y_t - y_0)^2);
                    r_tar = [r_tar r_t];
                end
                
                figure;
                x = -1000:1000;    y = x;
                suptitle(['Stimtype ' condnames{n} ', initial params: b_w = ' num2str(b_w0(k)) ',  b_v = ' num2str(b_v0(i)) ' [' num2str(model(count).optprs) ']']);
                subplot(1,2,1); plot(r_tar,r_sub,'.b');title('scatterplot of radius');hold on;
                hold on;plot(x,y,'r--');ylim([0 800]);xlim([0 800]);ylabel('r_s (cm)');xlabel('r_t (cm)');hold off;
                
                subplot(1,2,2); plot(theta_tar,theta_sub,'.k');title('scatterplot of \theta');hold on;
                hold on;plot(x,y,'r--');plot(x,-y,'r--');
                ylim([-100 100]);xlim([-90 90]);ylabel('\theta_s (degrees)');xlabel('\theta_t (degrees)');hold off;
            end
            count = count + 1;
        end
    end
    
end
% % plot trajectories
% figure;plot(trials(20).continuous.xmp,trials(20).continuous.ymp);hold on;plot(b(count).sim(20).x,b(count).sim(20).y,'r');
% title('trajectory comparison');plot(x_f(20),y_f(20),'rx');legend('actual trajectory','simulated trajectory','firefly');

toc;
    function [C,Ceq] = nonlincon(x)
        C(1) = x(1)*x(2);
        C(2) = -x(1)*x(2)-1;
        C(3) = x(3)*x(4);
        C(4) = -x(3)*x(4)-1;
        Ceq = [];
        