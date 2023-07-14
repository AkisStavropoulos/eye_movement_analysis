function [models,conditions] = fit_model(conditions,trls,models,expnum,modelname)
% exp 1/3: null, prior, leaky, prior_leaky
% exp 2: fixedprior, freeprior, fixedleak, freeleak
clear global;
global x_f y_f x_m0 y_m0 speed;

%% speeds and positions
% trls = conditions(end).trls;
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

%% optimise
switch expnum
    case {'1','2','3','4'}
        switch modelname
            case 'null'
                k=1; models(k).name = modelname;
                b_w0 = 1000; b_v0 = 1000;
                prs0 = [b_w0 b_v0];
                LB = [1 1]; UB = [inf inf];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                [models(k).optprs,models(k).fval,models(k).exitflag,...
                    models(k).output,models(k).lambda] = fmincon(@errfun_chancemodel2,prs0,[],[],[],[],LB,UB,[],opts);
            case 'prior'
                k=2; models(k).name = modelname;
                k_w0 = 1; k_v0 = 1;
                prs0 = [k_w0 k_v0];
                LB = [-2 -2]; UB = [2 2];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
                k_w = optprs(1);
                k_v = optprs(2);
                b_w0 = logspace(-1,2,10);
                b_v0 = logspace(-1,2,10);
                
                for i = 1:length(b_w0)
                    a_w0(i) = (k_w-1)/b_w0(i);
                    a_v0(i) = (k_v-1)/b_v0(i);
                    prs0 = [a_w0(i) b_w0(i) a_v0(i) b_v0(i)];
                    LB = [-1 1 -1 1]; UB = [1 inf 1 inf];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    [models(k).optprs,models(k).fval,models(k).exitflag,...
                        models(k).output,models(k).lambda] = fmincon(@errfun_bayesianmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
               figure;plot(trls);
                end
                
            case 'leaky'
                k=3; models(k).name = modelname;
                b_w0 = 10; b_v0 = 10;
                tau_r0 = 10; tau_theta0 = 10;
                prs0 = [b_w0 b_v0 tau_r0 tau_theta0];
                LB = [0 0 0 0]; UB = [inf inf inf inf];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                [models(k).optprs,models(k).fval,models(k).exitflag,...
                    models(k).output,models(k).lambda] = fmincon(@errfun_temporalleakymodel,prs0,[],[],[],[],LB,UB,[],opts);
            case 'prior_leaky'
                k=4; models(k).name = modelname;
                k_w0 = 1; k_v0 = 1;
                prs0 = [k_w0 k_v0];
                LB = [-2 -2]; UB = [2 2];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
                k_w = optprs(1);
                k_v = optprs(2);
                tau_r0 = 100; tau_theta0 = 100; b_w0 = 1000; b_v0 = 1000; a_w0 = (k_w-1)/b_w0; a_v0 = (k_v-1)/b_v0;
                prs0 = [a_w0 b_w0 a_v0 b_v0 tau_r0 tau_theta0];
                LB = [-1 1 -1 1 0 0]; UB = [0 inf 0 inf inf inf];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                [models(k).optprs,models(k).fval,models(k).exitflag,...
                    models(k).output,models(k).lambda] = fmincon(@errfun_combinationmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);               
        end
    case 2
        switch modelname
            case 'fixedprior'
                k=1; models(k).name = 'prior';
                %fit prior to all trials
                k_w0 = 1; k_v0 = 1;
                prs0 = [k_w0 k_v0];
                LB = [-2 -2]; UB = [2 2];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
                k_w = optprs(1);
                k_v = optprs(2);
                b_w0 = 1000; a_w0 = (k_w-1)/b_w0; b_v0 = 1000; a_v0 = (k_v-1)/b_v0;
                prs0 = [a_w0 b_w0 a_v0 b_v0];
                LB = [-1 1 -1 1]; UB = [0 inf 0 inf];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                [models(k).optprs,models(k).fval,models(k).exitflag,...
                    models(k).output,models(k).lambda] = fmincon(@errfun_bayesianmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
                % calculate variance for each condition
                a_w = models(k).optprs(1); a_v = models(k).optprs(3);
                for i=length(conditions)
                    ntrls = length(conditions(i).trls);
                    x_f = NaN(ntrls,1); y_f = NaN(ntrls,1);
                    x_m0 = NaN(ntrls,1); y_m0 = NaN(ntrls,1);
                    speed = struct('w',cell(ntrls,1),'v',cell(ntrls,1));
                    for j=1:ntrls
                        x_f(j) = mean(conditions(i).trls(j).ch(1).wf);
                        y_f(j) = mean(conditions(i).trls(j).ch(2).wf);
                        x_m0(j) = conditions(i).trls(j).ch(3).wf(1);
                        y_m0(j) = conditions(i).trls(j).ch(4).wf(1);
                        speed(j).w = conditions(i).trls(j).ch(5).wf;
                        speed(j).v = conditions(i).trls(j).ch(6).wf;
                    end
                    conditions(i).models(k).name = modelname;
                    k_w0 = 1; k_v0 = 1;
                    prs0 = [k_w0 k_v0];
                    LB = [-2 -2]; UB = [2 2];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
                    k_w = optprs(1); k_v = optprs(2);
                    a_w0 = a_w; b_w0 = (k_w-1)/a_w0; a_v0 = a_v; b_v0 = (k_v-1)/a_v0;
                    prs0 = [a_w0 b_w0 a_v0 b_v0];
                    LB = [a_w 1 a_v 1]; UB = [a_w inf a_v inf];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    [conditions(i).models(k).optprs,conditions(i).models(k).fval,conditions(i).models(k).exitflag,...
                        conditions(i).models(k).output,conditions(i).models(k).lambda] = fmincon(@errfun_bayesianmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
                end
            case 'freeprior'
                k=2;
                for i=1:length(conditions)
                    ntrls = length(conditions(i).trls);
                    x_f = NaN(ntrls,1); y_f = NaN(ntrls,1);
                    x_m0 = NaN(ntrls,1); y_m0 = NaN(ntrls,1);
                    speed = struct('w',cell(ntrls,1),'v',cell(ntrls,1));
                    for j=1:ntrls
                        x_f(j) = mean(conditions(i).trls(j).ch(1).wf);
                        y_f(j) = mean(conditions(i).trls(j).ch(2).wf);
                        x_m0(j) = conditions(i).trls(j).ch(3).wf(1);
                        y_m0(j) = conditions(i).trls(j).ch(4).wf(1);
                        speed(j).w = conditions(i).trls(j).ch(5).wf;
                        speed(j).v = conditions(i).trls(j).ch(6).wf;
                    end
                    conditions(i).models(k).name = modelname;
                    k_w0 = 1; k_v0 = 1;
                    prs0 = [k_w0 k_v0];
                    LB = [-2 -2]; UB = [2 2];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    optprs = fmincon(@errfun_scalingmodel,prs0,[],[],[],[],LB,UB,[],opts);
                    k_w = optprs(1);
                    k_v = optprs(2);
                    b_w0 = 1000; a_w0 = (k_w-1)/b_w0; b_v0 = 1000; a_v0 = (k_v-1)/b_v0;
                    prs0 = [a_w0 b_w0 a_v0 b_v0];
                    LB = [-1 1 -1 1]; UB = [1 inf 1 inf];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    [conditions(i).models(k).optprs,conditions(i).models(k).fval,conditions(i).models(k).exitflag,...
                        conditions(i).models(k).output,conditions(i).models(k).lambda] = fmincon(@errfun_bayesianmodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
                end
            case 'fixedleak'
                k=3; models(k).name = 'leaky';
                %fit timeconstant to all trials
                b_w0 = 1000; b_v0 = 1000;
                tau_r0 = 100; tau_theta0 = 100;
                prs0 = [b_w0 b_v0 tau_r0 tau_theta0];
                LB = [0 0 0 0]; UB = [inf inf inf inf];
                opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                [models(k).optprs,models(k).fval,models(k).exitflag,...
                    models(k).output,models(k).lambda] = fmincon(@errfun_temporalleakymodel,prs0,[],[],[],[],LB,UB,[],opts);
                % calculate variance for each condition
                tau_r = models(k).optprs(3); tau_theta = models(k).optprs(4);
                for i=1:length(conditions)
                    ntrls = length(conditions(i).trls);
                    x_f = NaN(ntrls,1); y_f = NaN(ntrls,1);
                    x_m0 = NaN(ntrls,1); y_m0 = NaN(ntrls,1);
                    speed = struct('w',cell(ntrls,1),'v',cell(ntrls,1));
                    for j=1:ntrls
                        x_f(j) = mean(conditions(i).trls(j).ch(1).wf);
                        y_f(j) = mean(conditions(i).trls(j).ch(2).wf);
                        x_m0(j) = conditions(i).trls(j).ch(3).wf(1);
                        y_m0(j) = conditions(i).trls(j).ch(4).wf(1);
                        speed(j).w = conditions(i).trls(j).ch(5).wf;
                        speed(j).v = conditions(i).trls(j).ch(6).wf;
                    end
                    conditions(i).models(k).name = modelname;
                    b_w0 = 1000; b_v0 = 1000; tau_r0 = tau_r; tau_theta0 = tau_theta;
                    prs0 = [b_w0 b_v0 tau_r0 tau_theta0];
                    LB = [0 0 tau_r tau_theta]; UB = [inf inf tau_r tau_theta];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    [conditions(i).models(k).optprs,conditions(i).models(k).fval,conditions(i).models(k).exitflag,...
                        conditions(i).models(k).output,conditions(i).models(k).lambda] = fmincon(@errfun_temporalleakymodel,prs0,[],[],[],[],LB,UB,@nonlincon,opts);
                end
            case 'freeleak'
                k=4;
                for i=1:length(conditions)
                    ntrls = length(conditions(i).trls);
                    x_f = NaN(ntrls,1); y_f = NaN(ntrls,1);
                    x_m0 = NaN(ntrls,1); y_m0 = NaN(ntrls,1);
                    speed = struct('w',cell(ntrls,1),'v',cell(ntrls,1));
                    for j=1:ntrls
                        x_f(j) = mean(conditions(i).trls(j).ch(1).wf);
                        y_f(j) = mean(conditions(i).trls(j).ch(2).wf);
                        x_m0(j) = conditions(i).trls(j).ch(3).wf(1);
                        y_m0(j) = conditions(i).trls(j).ch(4).wf(1);
                        speed(j).w = conditions(i).trls(j).ch(5).wf;
                        speed(j).v = conditions(i).trls(j).ch(6).wf;
                    end
                    conditions(i).models(k).name = modelname;
                    b_w0 = 1000; b_v0 = 1000;
                    tau_r0 = 100; tau_theta0 = 100;
                    prs0 = [b_w0 b_v0 tau_r0 tau_theta0];
                    LB = [0 0 0 0]; UB = [inf inf inf inf];
                    opts = optimoptions(@fmincon,'maxiter',500,'Display','iter','Algorithm','interior-point');
                    [conditions(i).models(k).optprs,conditions(i).models(k).fval,conditions(i).models(k).exitflag,...
                        conditions(i).models(k).output,conditions(i).models(k).lambda] = fmincon(@errfun_temporalleakymodel,prs0,[],[],[],[],LB,UB,[],opts);
                end
        end
end

    function [C,Ceq] = nonlincon(x)
        C(1) = x(1)*x(2);
        C(2) = -x(1)*x(2)-1;
        C(3) = x(3)*x(4);
        C(4) = -x(3)*x(4)-1;
        Ceq = [];