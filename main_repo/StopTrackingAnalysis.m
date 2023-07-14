%%
%%
%% ADDITIONAL ANALYSIS OF SACCADE RELATIONSHIP WITH END OF STEERING
%%
%%
%% Distribution of saccades around end of steering/trial/button push
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
twin = 4;
condition = {'vestibular','visual','combined'};

tsac_rel = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlind = poolindx{i,s};
        
        tmp1 = []; tmp2 = []; tmp3 = []; tmp4 = [];
        for j = 1:numel(trlind)
            
        % trial end
        try tmp1 = [tmp1 subject(i).trials(trlind(j)).events.t_sac' - subject(i).trials(trlind(j)).events.t_end];
        catch; end
        % successful button push
        try tmp2 = [tmp2 subject(i).trials(trlind(j)).events.t_sac' - subject(i).trials(trlind(j)).events.button];
        catch; end
        % first button push (remember, sometimes 1st button push is the successful push)
        try push1st = subject(i).trials(trlind(j)).events.push_on(1);
            tmp3 = [tmp3  subject(i).trials(trlind(j)).events.t_sac' - (push1st + subject(i).trials(trlind(j)).events.t_beg)];
        catch; end
        % steering end (joystick let-go)
        jsy = -subject(i).trials(trlind(j)).mc.JS_X_Raw;           jsx = subject(i).trials(trlind(j)).mc.JS_Yaw_Raw;
        vmax = subject(i).trials(trlind(j)).prs.vmax;            wmax = subject(i).trials(trlind(j)).prs.wmax;
        try letgo_t = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
            tmp4 = [tmp4 subject(i).trials(trlind(j)).events.t_sac' - (subject(i).trials(trlind(j)).continuous.ts(letgo_t)+subject(i).trials(trlind(j)).events.t_beg)];
        catch; end
        end
        tsac_rel.t_end{i,s} = tmp1(abs(tmp1) < twin);
        tsac_rel.button_success{i,s} = tmp2(abs(tmp2) < twin);
        tsac_rel.button_first{i,s} = tmp3(abs(tmp3) < twin);
        tsac_rel.js_end{i,s} = tmp4(abs(tmp4) < twin);

        if 1
        % sanity check
        nx = 40;
        a = figure(111); clf; a.Position = [95 64 404 931];
        subplot(4,1,1); hold on;
        histogram(tsac_rel.t_end{i,s},nx); vline(0,'r'); xlabel('time from trial end [s]'); title('trial end');
        subplot(4,1,2); hold on;
        histogram(tsac_rel.button_success{i,s},nx); vline(0,'r'); xlabel('time from successful button push [s]'); title('successful button push');
        subplot(4,1,3); hold on;
        histogram(tsac_rel.button_first{i,s},nx); vline(0,'r'); xlabel('time from 1st button push [s]'); title('1st button push');
        subplot(4,1,4); hold on;
        histogram(tsac_rel.js_end{i,s},nx); vline(0,'r'); xlabel('time from JS let-go [s]'); title('Joystick let-go');
        suptitle([subject(i).name ' - ' condition{s}])
        end
    end
end


%% Difference between JS let-go and 1st push

default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
twin = 4;
condition = {'vestibular','visual','combined'};

t_rel = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlind = poolindx{i,s};
        
        tmp1 = []; tmp2 = []; tmp3 = []; tmp11 = []; tmp22 = []; tmp33 = []; tau = [];
        for j = 1:numel(trlind)
        % compute let-go
        jsy = -subject(i).trials(trlind(j)).mc.JS_X_Raw;           jsx = subject(i).trials(trlind(j)).mc.JS_Yaw_Raw;
        vmax = subject(i).trials(trlind(j)).prs.vmax;            wmax = subject(i).trials(trlind(j)).prs.wmax;
        try letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
            letgo_t = subject(i).trials(trlind(j)).continuous.ts(letgo_i) + subject(i).trials(trlind(j)).events.t_beg;
        catch; end
        % Let-go 2 trial end
        try tmp1 = [tmp1 letgo_t - subject(i).trials(trlind(j)).events.t_end];
        catch; end
        % Let-go 2 successful button push
        try tmp2 = [tmp2 letgo_t - subject(i).trials(trlind(j)).events.button];
        catch; end
        % Let-go 2 first button push
        try push1st = subject(i).trials(trlind(j)).events.push_on(1) + subject(i).trials(trlind(j)).events.t_beg;
            tmp3 = [tmp3  letgo_t - push1st];
            tau = [tau subject(i).trials(trlind(j)).prs.tau];
        catch; end
        
        % 1st push 2 trial end
        try tmp11 = [tmp11 push1st - subject(i).trials(trlind(j)).events.t_end];
        catch; end
        % 1st push 2 successful button push
        try tmp22 = [tmp22 push1st - subject(i).trials(trlind(j)).events.button];
        catch; end
        % 1st push 2 let-go
        try tmp33 = [tmp33  push1st - letgo_t];
        catch; end

        end
        % let-go
        t_rel.js2t_end{i,s} = tmp1(abs(tmp1) < twin);
        t_rel.js2button_success{i,s} = tmp2(abs(tmp2) < twin);
        t_rel.js2button_first{i,s} = tmp3(abs(tmp3) < twin); tau = tau(abs(tmp3) < twin);
        
        % 1st push
        t_rel.push2t_end{i,s} = tmp11(abs(tmp11) < twin);
        t_rel.push2button_success{i,s} = tmp22(abs(tmp22) < twin);
        t_rel.push2js{i,s} = tmp33(abs(tmp33) < twin);
        

        if 1
        % sanity check
        nx = 20;
        a = figure(222); clf; a.Position = [11 474 510 519];
        subplot(2,2,1); hold on;
        histogram(t_rel.js2button_success{i,s},nx); vline(0,'r'); xlabel('time from successful button push [s]'); title('Let-go to successful button push');
        subplot(2,2,3); hold on;
        histogram(t_rel.js2button_first{i,s},nx); vline(0,'r'); xlabel('time from 1st button push [s]'); title('Let-go to 1st button push');
        
        
        subplot(2,2,2); hold on;
        histogram(t_rel.push2button_success{i,s},nx); vline(0,'r'); xlabel('time from successful button push [s]'); title('1st push to successful button push');
        subplot(2,2,4); hold on;
        plot(tau,t_rel.js2button_first{i,s},'.'); xlabel('tau [s]'); ylabel('let-go 2 1st push [s]'); title('let-go to 1st push');
        
        suptitle([subject(i).name ' - ' condition{s}])
        end
    end
end


%% Relationship of final saccades and 1st push/let-go timing
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
sacwin = 6;
sacthresh = 10;
twin = 4;
wtf_time = 20;
condition = {'vestibular','visual','combined'};

data = extractEndOfSteeringVars(subject,params,jsthresh,sacwin,sacthresh);
push_success = data.push_success;
push_1st = data.push_1st;
tsac = data.tsac;
letgo = data.letgo;


rho = []; p = [];
for i = 1:Nsubs
    for s = 1:Nstim
               
        if 1
        % sanity check
        nx = 40;
        a = figure(111); clf; a.Position = [6 694 929 285];

        subplot(1,3,1); hold on;
        [rho.push_success(i,s),p.push_success(i,s)] = nancorr(push_success{i,s}(:),tsac{i,s}(:));
        plot(push_success{i,s},tsac{i,s},'.'); xlabel('successful push [s]'); ylabel('1st large saccade time [s]'); title('Successful push');
        
        subplot(1,3,2); hold on;
        [rho.push_1st(i,s),p.push_1st(i,s)] = nancorr(push_1st{i,s}(:),tsac{i,s}(:));
        plot(push_1st{i,s},tsac{i,s},'.'); xlabel('1st push [s]'); ylabel('1st large saccade time [s]'); title('1st push');
        
        subplot(1,3,3); hold on;
        [rho.letgo(i,s),p.letgo(i,s)] = nancorr(letgo{i,s}(:),tsac{i,s}(:));
        plot(letgo{i,s},tsac{i,s},'.'); xlabel('JS let-go [s]'); ylabel('1st large saccade time [s]'); title('JS let-go');
        
        suptitle([subject(i).name ' - ' condition{s}])
        end
    end
end

figure;
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
for s = 1:Nstim
subplot(1,3,1); hold on;
bar(s,nanmean(rho.push_success(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(rho.push_success(:,s)),nanstd(rho.push_success(:,s))./sqrt(Nsubs),'k','capsize',0);
xticks(1:Nstim); xticklabels(condition); ylabel('saccade time and successful push corr.'); title('Successful push'); ylim([0 1]);

subplot(1,3,2); hold on;
bar(s,nanmean(rho.push_1st(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(rho.push_1st(:,s)),nanstd(rho.push_1st(:,s))./sqrt(Nsubs),'k','capsize',0);
xticks(1:Nstim); xticklabels(condition); ylabel('saccade time and 1st push corr.'); title('1st push'); ylim([0 1]);

subplot(1,3,3); hold on;
bar(s,nanmean(rho.letgo(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(rho.letgo(:,s)),nanstd(rho.letgo(:,s))./sqrt(Nsubs),'k','capsize',0);
xticks(1:Nstim); xticklabels(condition); ylabel('saccade time and JS let-go corr.'); title('JS let-go'); ylim([0 1]);

end

%% Distribution of final saccades and 1st push/let-go timing
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
sacwin = 6;
sacthresh = 10;
wtf_time = 20;
condition = {'vestibular','visual','combined'};

dt_rel = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlind = poolindx{i,s};
        
        tau = []; tsac_mu = []; t_beg = []; t_end = []; push_success = []; push_1st = []; letgo_t = [];
        for j = 1:numel(trlind)
        
        t_beg(j) = subject(i).trials(trlind(j)).events.t_beg;
        t_end(j) = subject(i).trials(trlind(j)).events.t_end - t_beg(j); 
        % succesful push
        try push_success(j) = subject(i).trials(trlind(j)).events.button - t_beg(j); catch; push_success(j) = nan; end
        % first push
        try push_1st(j) = subject(i).trials(trlind(j)).events.push_on(1); catch; push_1st(j) = nan; end
        % steering end (joystick let-go)
        jsy = -subject(i).trials(trlind(j)).mc.JS_X_Raw;           jsx = subject(i).trials(trlind(j)).mc.JS_Yaw_Raw;
        vmax = subject(i).trials(trlind(j)).prs.vmax;            wmax = subject(i).trials(trlind(j)).prs.wmax;
        try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                letgo_t(j) = subject(i).trials(trlind(j)).continuous.ts(letgo_i);
        catch;  letgo_t(j) = nan; end
        % get taus
        tau(j) = subject(i).trials(trlind(j)).prs.tau;
        % saccade parameters
        t_sac = subject(i).trials(trlind(j)).events.t_sac - t_beg(j);
        sac_mag = subject(i).trials(trlind(j)).events.sac_mag;
        % max of final saccades [0-sacwin] sec before end of trial
        try     t_indx = t_sac > t_end(j) - sacwin; 
                sac_mag_tmp = sac_mag(t_indx);
                bigsac = sac_mag_tmp(find(sac_mag_tmp > sacthresh,1)); % first big saccade within [sacwin] from trial end
%                 bigsac = max(sac_mag_tmp); % biggest saccade within [sacwin] from trial end
                tsac_tmp = t_sac(sac_mag==bigsac);
                tsac_mu(j) = nanmean(tsac_tmp(tsac_tmp > t_end(j) - sacwin)); 
        catch;  tsac_mu(j) = nan;
        end
                
        end
        t_end(t_end>wtf_time)=nan;
        push_success(push_success>wtf_time)=nan;
        push_1st(push_1st>wtf_time)=nan;
        letgo_t(letgo_t>wtf_time)=nan;
        tsac_mu(tsac_mu>wtf_time)=nan;
        
        dt_rel.t_end{i,s} = tsac_mu - t_end;
        dt_rel.push_success{i,s} = tsac_mu - push_success;
        dt_rel.push_1st{i,s} = tsac_mu - push_1st;
        dt_rel.letgo{i,s} = tsac_mu - letgo_t;
        
        if 1
        % sanity check
        nx = 30;
        a = figure(111); clf; a.Position = [95 64 404 931];

        subplot(4,1,1); hold on;
        histogram(dt_rel.push_success{i,s},nx); vline(0,'r'); xlabel('time to successful push [s]'); ylabel('No. trials (or saccades)'); title('Successful push');
        
        subplot(4,1,2); hold on;
        histogram(dt_rel.push_1st{i,s},nx); vline(0,'r'); xlabel('time to 1st push [s]'); ylabel('No. trials (or saccades)'); title('1st push');
        
        subplot(4,1,3); hold on;
        histogram(dt_rel.letgo{i,s},nx); vline(0,'r'); xlabel('time to JS let-go [s]'); ylabel('No. trials (or saccades)'); title('JS let-go');
        
        suptitle([subject(i).name ' - ' condition{s}])
        end
    end
end

figure;
for s = 1:Nstim
subplot(2,1,1); hold on;
cellfun(@(x) plot(s+randn*0.1,nanmean(x),'o','markerfacecolor',colr(s,:),'markeredgecolor','none'),dt_rel.push_1st(:,s)); 
xticks(1:Nstim); xticklabels(condition); ylabel('saccade from 1st push [s]'); axis([0 4 -5 5]); hline(0,'k'); 

subplot(2,1,2); hold on;
cellfun(@(x) plot(s+randn*0.1,nanmean(x),'o','markerfacecolor',colr(s,:),'markeredgecolor','none'),dt_rel.letgo(:,s)); 
xticks(1:Nstim); xticklabels(condition); ylabel('saccade from JS let-go [s]'); axis([0 4 -5 5]); hline(0,'k'); 

end

%% Relationship of tau and final saccades timing
subject = subject_backup(keepindx);
subject = subject_backup;
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
sacwin = 6;
sacthresh = 10;
twin = 4;
condition = {'vestibular','visual','combined'};

tsac_rel = []; rho = []; p = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlind = poolindx{i,s};
        
        tau = []; tsac_mu = [];
        for j = 1:numel(trlind)
        
        t_end = subject(i).trials(trlind(j)).events.t_end;
        t_sac = subject(i).trials(trlind(j)).events.t_sac;
        sac_mag = subject(i).trials(trlind(j)).events.sac_mag;
        % get taus
        tau(j) = subject(i).trials(trlind(j)).prs.tau;
        % max of final saccades [0-4] sec before end of trial
        try     t_indx = t_sac > t_end - sacwin; 
                sac_mag_tmp = sac_mag(t_indx);
                bigsac = sac_mag_tmp(find(sac_mag_tmp > sacthresh,1)); % first big saccade within [sacwin] from trial end
%                 bigsac = max(sac_mag_tmp); % biggest saccade within [sacwin] from trial end
                tsac_tmp = t_sac(sac_mag==bigsac);
                tsac_mu(j) = nanmean(tsac_tmp(tsac_tmp > t_end - sacwin)); 
        catch;  tsac_mu(j) = nan;
        end
        % trial end
        try     tsac_rel.t_end{i,s}(j) = tsac_mu(j) - t_end;
        catch;  tsac_rel.t_end{i,s}(j) = nan;
        end
        % successful button push
        try     tsac_rel.button_success{i,s}(j) = tsac_mu(j) - subject(i).trials(trlind(j)).events.button;
        catch;  tsac_rel.button_success{i,s}(j) = nan;
        end
        % first button push (remember, sometimes 1st button push is the successful push)
        try push1st = subject(i).trials(trlind(j)).events.push_on(1) + subject(i).trials(trlind(j)).events.t_beg;
                tsac_rel.button_first{i,s}(j) = tsac_mu(j) - push1st;
        catch;  tsac_rel.button_first{i,s}(j) = nan;
        end
        % steering end (joystick let-go)
        jsy = -subject(i).trials(trlind(j)).mc.JS_X_Raw;           jsx = subject(i).trials(trlind(j)).mc.JS_Yaw_Raw;
        vmax = subject(i).trials(trlind(j)).prs.vmax;            wmax = subject(i).trials(trlind(j)).prs.wmax;
        try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                letgo_t = subject(i).trials(trlind(j)).continuous.ts(letgo_i) + subject(i).trials(trlind(j)).events.t_beg;
                tsac_rel.js_end{i,s}(j) = tsac_mu(j) - letgo_t;
        catch;  tsac_rel.js_end{i,s}(j) = nan;
        end
        
        end
        
        if 1
        % sanity check
        nx = 40;
        a = figure(111); clf; a.Position = [477 474 404 505];
        subplot(2,1,1); hold on;
        [rho.push(i,s),p.push(i,s)] = nancorr(tau(:),tsac_rel.button_first{i,s}(:));
        plot(tau,tsac_rel.button_first{i,s},'.'); xlabel('tau [s]'); ylabel('saccades from 1st button push [s]'); title('1st button push');ylim([-10 10]);

        subplot(2,1,2); hold on;
        [rho.letgo(i,s),p.letgo(i,s)] = nancorr(tau(:),tsac_rel.js_end{i,s}(:));
        plot(tau,tsac_rel.js_end{i,s},'.'); xlabel('tau [s]'); ylabel('saccades from JS let-go [s]'); title('Joystick let-go');ylim([-10 10]);
        suptitle([subject(i).name ' - ' condition{s}])
        end
    end
end

figure;
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
for s = 1:Nstim
subplot(2,1,1); hold on;
bar(s,nanmean(rho.push(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(rho.push(:,s)),nanstd(rho.push(:,s))./sqrt(Nsubs),'k','capsize',0);
xticks(1:Nstim); xticklabels(condition); ylabel('saccade timing and tau corr.'); title('1st button push'); ylim([0 1]);

subplot(2,1,2); hold on;
bar(s,nanmean(rho.letgo(:,s)),'facecolor',colr(s,:)); errorbar(s,nanmean(rho.letgo(:,s)),nanstd(rho.letgo(:,s))./sqrt(Nsubs),'k','capsize',0);
xticks(1:Nstim); xticklabels(condition); ylabel('saccade timing and tau corr.'); title('Joystick let-go'); ylim([0 1]);
end


%% Regression analysis of final saccade and JS let-go, 1st push, tau, and interaction terms
subject = subject_backup(keepindx);
subject = subject_backup(setdiff([1:numel(subject_backup)],[3 5 13]));
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
sacwin = 1/2;
sacthresh = 10;
condition = {'vestibular','visual','combined'};

data = extractEndOfSteeringVars(subject,params,jsthresh,sacwin,sacthresh);
push_success = data.push_success;
push_1st = data.push_1st;
tsac = data.tsac;
letgo = data.letgo;
tau = data.tau;

b = []; b_int = []; stats = []; pval = []; finalmodel = []; history = []; 
for i = 1:Nsubs
    for s = 1:Nstim
        Ntrls = numel(tsac{i,s});
        
        t_letGo = letgo{i,s}(:);
        t_1stPush = push_1st{i,s}(:);
        taus = tau{i,s}(:);
        t_letGoXtau = t_letGo.*taus;
        t_1stPushXtau = t_1stPush.*taus;
        
        X = [t_1stPush t_1stPushXtau]; % [t_letGo t_1stPush taus t_letGoXtau t_1stPushXtau];
        y = tsac{i,s}(:);
        
        [b{i,s},~,pval{i,s},finalmodel{i,s},stats{i,s},nextstep,history{i,s}] = stepwisefit(X,y,'PRemove',0.05); % yhat = stats.intercept + X(:,finalmodel)*b(finalmodel)
        
    end
end

model_select = [];
for s = 1:Nstim
    model_select(s,:) = sum(cell2mat(cellfun(@(x) x, finalmodel(:,s),'un',0)))/Nsubs;
end

model_select = table(model_select);
model_select = splitvars(model_select);
model_select.Properties.VariableNames = {'Button','ButtonXtau'} % {'letGo','Button','tau','letGoXtau','ButtonXtau'}

% seems like the model with just [Button, ButtonXtau] is adequate, based on the F-statistic
fstat_mu = nanmean(cellfun(@(x) x.fstat,stats));    fstat_sd = nanstd(cellfun(@(x) x.fstat,stats));
fstat_pval = sum(cellfun(@(x) x.pval<0.05,stats));
disp(['F-stat mean: ' num2str(fstat_mu)])
disp(['F-stat std: ' num2str(fstat_sd)]); disp('===')

USSR_mu = nanmean(cellfun(@(x) x.SSresid,stats));    USSR_sd = nanstd(cellfun(@(x) x.SSresid,stats));
disp(['SSResid mean: ' num2str(USSR_mu)])
disp(['SSResid std: ' num2str(USSR_sd)])

%% Sanity check of multiple regression: Individual trials 
sacwin = 6;
sacthresh = 10;
jsthresh = 0.1;

regrcoefs = cellfun(@(x1,x2) [x1.intercept ; x2],stats,b,'un',0);

% check whether button was recorded
subject_backup = buttonRecFlag(subject_backup);


% Individual trials
a = figure(345); a.Position = [12 327 320 575];
for i = 1:Nsubs
    for s = 1:Nstim
        trlindx = poolindx{i,s};
        for k = 1:numel(trlindx)
            ts = subject(i).trials(trlindx(k)).continuous.ts;
            tau = subject(i).trials(trlindx(k)).prs.tau;
            
            try     t_push = subject(i).trials(trlindx(k)).events.push_on(1);
                    t_end = subject(i).trials(trlindx(k)).events.t_end;
                    t_beg = subject(i).trials(trlindx(k)).events.t_beg;
                    t_sac = subject(i).trials(trlindx(k)).events.t_sac;
                    sac_mag = subject(i).trials(trlindx(k)).events.sac_mag;
                    
                    try     t_indx = t_sac > t_end - sacwin; % t_indx = t_sac > t_end(j) - sacwin;
                            sac_mag_tmp = sac_mag(t_indx);
                            bigsacmag = sac_mag_tmp(find(sac_mag_tmp > sacthresh,1));
                            bigsac = t_sac(sac_mag==bigsacmag) - t_beg;
                            novline = 0;
                    catch;  novline = 1; end
                    
                    % steering end (joystick let-go)
                    jsy = -subject(i).trials(trlindx(k)).mc.JS_X_Raw;           jsx = subject(i).trials(trlindx(k)).mc.JS_Yaw_Raw;
                    vmax = subject(i).trials(trlindx(k)).prs.vmax;            wmax = subject(i).trials(trlindx(k)).prs.wmax;
                    try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                        t_letgo = subject(i).trials(trlindx(k)).continuous.ts(letgo_i);
                    catch;  t_letgo = nan; end
                    
                    % get eye position
                    ye = 0.5*(subject(i).trials(trlindx(k)).continuous.yle + subject(i).trials(trlindx(k)).continuous.yre);
                    ze = 0.5*(subject(i).trials(trlindx(k)).continuous.zle + subject(i).trials(trlindx(k)).continuous.zre);
                    
                    % get target position
                    phi = subject(i).trials(trlindx(k)).continuous.phi;
                    x_monk = subject(i).trials(trlindx(k)).continuous.xmp;	y_monk = subject(i).trials(trlindx(k)).continuous.ymp;
                    fireflyposx = subject(i).trials(trlindx(k)).prs.fireflyposx;   fireflyposy = subject(i).trials(trlindx(k)).prs.fireflyposy;
                    x_fly = fireflyposx - x_monk;     y_fly = fireflyposy - y_monk;
                    R = @(phi) [cosd(phi)   sind(phi) ; -sind(phi)  cosd(phi)];
                    XY_fly = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi, x_fly, y_fly,'uniformoutput',false));
                    x_fly = XY_fly(:,1);    y_fly = XY_fly(:,2);
                    
                    y = y_fly; y(y < 0) = nan;
                    x = x_fly;
                    z = -prs.height;
                    delta = prs.interoculardist;
                    [yle_targ, zle_targ, yre_targ, zre_targ] = world2eye(x,y,z,delta);
                    yt = 0.5*(yle_targ + yre_targ);	zt = 0.5*(zle_targ + zre_targ);
                    
                    if numel(unique(ye)) > 1
                    % plot
                    X = [1 t_push t_push*tau]; % 
                    X = [1 t_letgo t_push tau t_letgo*tau t_push*tau];
                    
                    t_break = X*regrcoefs{i,s};
                    clf;
                    subplot(2,1,1); hold on;
                    plot(ts,ye,'k'); plot(ts,yt,'--','color',[.5 .5 .5]); vline(t_break,'r'); vline(t_push,'b'); if ~novline; vline(bigsac,'g');end
                    xlabel('time [s]'); ylabel('HOR [deg]'); title([condition{s} ': Trial = ' num2str(trlindx(k))]);
                    
                    subplot(2,1,2); hold on;
                    plot(ts,ze,'k'); plot(ts,zt,'--','color',[.5 .5 .5]); vline(t_break,'r'); vline(t_push,'b'); if ~novline; vline(bigsac,'g');end
                    plot(ts,movmean(ze,80),'r')
                    xlabel('time [s]'); ylabel('VER [deg]');
                    suptitle(subject(i).name); 
                    end
            catch
            end
            
        end
    end
end

%% Sanity check of multiple regression: average vertical eye position
maxlength = 1000;
pre_break = 500;
post_break = pre_break;
jsthresh = 0.1;
sacwin = 1/2;
velthresh = 10;
avgwin = 80;

regrcoefs = cellfun(@(x1,x2) [x1.intercept ; x2],stats,b,'un',0);

ze_break = []; ze_sacbreak = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlindx = poolindx{i,s};
        for k = 1:numel(trlindx)
            ts = subject(i).trials(trlindx(k)).continuous.ts;
            tau = subject(i).trials(trlindx(k)).prs.tau;
            t_end = subject(i).trials(trlindx(k)).events.t_end - subject(i).trials(trlindx(k)).events.t_beg;
            
            try     t_push = subject(i).trials(trlindx(k)).events.push_on(1);
            catch   button_rec = subject(i).trials(trlindx(k)).prs.button_rec;
                if button_rec
                    t_push = subject(i).trials(trlindx(k)).events.button;
                    if isempty(t_push); t_push = nan; end
                else
                    t_push = nan;
                end
            end
            % steering end (joystick let-go)
            jsy = -subject(i).trials(trlindx(k)).mc.JS_X_Raw;           jsx = subject(i).trials(trlindx(k)).mc.JS_Yaw_Raw;
            vmax = subject(i).trials(trlindx(k)).prs.vmax;            wmax = subject(i).trials(trlindx(k)).prs.wmax;
            try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                t_letgo = subject(i).trials(trlindx(k)).continuous.ts(letgo_i);
                if isempty(t_letgo); t_letgo = nan; end
            catch;  t_letgo = nan; end
            
            % get eye position
            ze = 0.5*(subject(i).trials(trlindx(k)).continuous.zle + subject(i).trials(trlindx(k)).continuous.zre);
            
            % get tracking-break time from actual saccade
            t_indx = diff(movmean(ze,avgwin))/dt > velthresh;
            tsac_tmp = ts(t_indx);
            try
                if sacwin <= 1; tsac_break = tsac_tmp(find(tsac_tmp > t_end - sacwin*t_end,1));
                else; tsac_break = tsac_tmp(find(tsac_tmp > t_end - sacwin,1));
                end
                sacbreakindx = find(ts <= tsac_break,1,'last');
                pre = sacbreakindx-pre_break; if pre<1; pre = 1; end
                post = sacbreakindx+post_break; if post>numel(ze); post=numel(ze); end
                ze_sacbreak{i,s}(k,:) = [nan(maxlength-(sacbreakindx-pre),1) ; ze(pre:post) ; nan(maxlength-(post-sacbreakindx),1) ];
            catch
                ze_sacbreak{i,s}(k,:) = nan(2*maxlength+1,1);
            end

            % get tracking-break time from model
            X = [1 t_push t_push*tau]; % 
            X = [1 t_letgo t_push tau t_letgo*tau t_push*tau];
            t_break = X*regrcoefs{i,s};
            
            breakindx = find(ts <= t_break,1,'last');
            if ~isempty(breakindx)
            pre = breakindx-pre_break; if pre<1; pre = 1; end
            post = breakindx+post_break; if post>numel(ze); post=numel(ze); end
            ze_break{i,s}(k,:) = [nan(maxlength-(breakindx-pre),1) ; ze(pre:post) ; nan(maxlength-(post-breakindx),1) ];
            else
               ze_break{i,s}(k,:) = nan(2*maxlength+1,1);
            end
        end
    end
end
%
dt = diff(subject(1).trials(1).continuous.ts(1:2));
ts = (-maxlength:maxlength)*dt;
aa = figure(345); aa.Position = [12 327 408 575];
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
for i = 1:Nsubs
        clf;
    for s = 1:Nstim
        % plot
        subplot(Nstim,2,s+(s-1)); hold on;
        plot(ts,nanmean(ze_sacbreak{i,s}),'color',colr(s,:)); axis([-3 3 -45 0]); vline(0,'k:');
        xlabel('time from tracking-break [s]'); ylabel('VER eye pos [deg]'); title(['NO model - ' condition{s}]); 
        suptitle(subject(i).name);
            
        subplot(Nstim,2,2*s); hold on;
        plot(ts,nanmean(ze_break{i,s}),'color',colr(s,:)); axis([-3 3 -45 0]); vline(0,'k:');
        xlabel('time from tracking-break [s]'); ylabel('VER eye pos [deg]'); title(['model - ' condition{s}]); 
        suptitle(subject(i).name);

    end
end
%% Sanity check of multiple regression: TTE
maxlength = 1000;
pre_break = 500;
post_break = pre_break;
jsthresh = 0.1;

regrcoefs = cellfun(@(x1,x2) [x1.intercept ; x2],stats,b,'un',0);

tte_model = []; tte_true = [];
for i = 1:Nsubs
    for s = 1:Nstim
        trlindx = poolindx{i,s};
        for k = 1:numel(trlindx)
            ts = subject(i).trials(trlindx(k)).continuous.ts;
            tau = subject(i).trials(trlindx(k)).prs.tau;
            t_end = subject(i).trials(trlindx(k)).events.t_end - subject(i).trials(trlindx(k)).events.t_beg;
                        
            try     t_push = subject(i).trials(trlindx(k)).events.push_on(1);
            catch   button_rec = subject(i).trials(trlindx(k)).prs.button_rec;
                if button_rec
                    t_push = subject(i).trials(trlindx(k)).events.button;
                    if isempty(t_push); t_push = nan; end
                else
                    t_push = nan;
                end
            end
            % steering end (joystick let-go)
            jsy = -subject(i).trials(trlindx(k)).mc.JS_X_Raw;           jsx = subject(i).trials(trlindx(k)).mc.JS_Yaw_Raw;
            vmax = subject(i).trials(trlindx(k)).prs.vmax;            wmax = subject(i).trials(trlindx(k)).prs.wmax;
            try     letgo_i = find(abs(jsy./vmax) > jsthresh & abs(jsx./wmax) > jsthresh,1,'last');
                t_letgo = subject(i).trials(trlindx(k)).continuous.ts(letgo_i);
                if isempty(t_letgo); t_letgo = nan; end
            catch;  t_letgo = nan; end
            
            % get eye position
            ye = 0.5*(subject(i).trials(trlindx(k)).continuous.yle + subject(i).trials(trlindx(k)).continuous.yre);
            ze = 0.5*(subject(i).trials(trlindx(k)).continuous.zle + subject(i).trials(trlindx(k)).continuous.zre);
            
            % get target position
            phi = subject(i).trials(trlindx(k)).continuous.phi;
            x_monk = subject(i).trials(trlindx(k)).continuous.xmp;	y_monk = subject(i).trials(trlindx(k)).continuous.ymp;
            fireflyposx = subject(i).trials(trlindx(k)).prs.fireflyposx;   fireflyposy = subject(i).trials(trlindx(k)).prs.fireflyposy;
            x_fly = fireflyposx - x_monk;     y_fly = fireflyposy - y_monk;
            R = @(phi) [cosd(phi)   sind(phi) ; -sind(phi)  cosd(phi)];
            XY_fly = cell2mat(arrayfun(@(phi,x,y) [x y]*R(phi), phi, x_fly, y_fly,'uniformoutput',false));
            x_fly = XY_fly(:,1);    y_fly = XY_fly(:,2);
            
            y = y_fly; y(y < 0) = nan;
            x = x_fly;
            z = -prs.height;
            delta = prs.interoculardist;
            [yle_targ, zle_targ, yre_targ, zre_targ] = world2eye(x,y,z,delta);
            yt = 0.5*(yle_targ + yre_targ);	zt = 0.5*(zle_targ + zre_targ);
            
            % get tracking-break time from actual saccade
            t_indx = diff(movmean(ze,avgwin))/dt > velthresh;
            tsac_tmp = ts(t_indx);
            try
                if sacwin <= 1; tsac_break = tsac_tmp(find(tsac_tmp > t_end - sacwin*t_end,1));
                else; tsac_break = tsac_tmp(find(tsac_tmp > t_end - sacwin,1));
                end
                sacbreakindx = find(ts <= tsac_break,1,'last');
                pre = sacbreakindx-pre_break; if pre<1; pre = 1; end
                post = sacbreakindx+post_break; if post>numel(ze); post=numel(ze); end
                tte_true{i,s}(k,:) = [nan(maxlength-(sacbreakindx-pre),1) ; sqrt((yt(pre:post)-ye(pre:post)).^2 + (zt(pre:post)-ze(pre:post)).^2) ; nan(maxlength-(post-sacbreakindx),1) ];
            catch
                tte_true{i,s}(k,:) = nan(2*maxlength+1,1);
            end

            % get tracking-break time from model
            X = [1 t_push t_push*tau]; % 
            X = [1 t_letgo t_push tau t_letgo*tau t_push*tau];
            t_break = X*regrcoefs{i,s};
            
            breakindx = find(ts <= t_break,1,'last');
            if ~isempty(breakindx)
            pre = breakindx-pre_break; if pre<1; pre = 1; end
            post = breakindx+post_break; if post>numel(yt); post=numel(yt); end
            tte_model{i,s}(k,:) = [nan(maxlength-(breakindx-pre),1) ; sqrt((yt(pre:post)-ye(pre:post)).^2 + (zt(pre:post)-ze(pre:post)).^2) ; nan(maxlength-(post-breakindx),1) ];
            else
               tte_model{i,s}(k,:) = nan(2*maxlength+1,1); 
            end
            
        end
    end
end

dt = diff(subject(1).trials(1).continuous.ts(1:2));
ts = (-maxlength:maxlength)*dt;
aa = figure(345); aa.Position = [12 327 408 575];
colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
for i = 1:Nsubs
        clf;
    for s = 1:Nstim
        % plot
        subplot(Nstim,2,s+(s-1)); hold on;
        plot(ts,nanmean(tte_true{i,s}),'color',colr(s,:)); axis([-3 3 0 50]); vline(0,'k:');
        xlabel('time from tracking-break [s]'); ylabel('TTE [deg]'); title(['NO model - ' condition{s}]); 
        suptitle(subject(i).name);
            
        subplot(Nstim,2,2*s); hold on;
        plot(ts,nanmean(tte_model{i,s}),'color',colr(s,:)); axis([-3 3 0 50]); vline(0,'k:');
        xlabel('time from tracking-break [s]'); ylabel('TTE [deg]'); title(['model - ' condition{s}]); 
        suptitle(subject(i).name);
    end
end
%% Regression between let-go and 1st button push from trial end
subject = subject_backup(keepindx);
subject = subject_backup(setdiff([1:numel(subject_backup)],[3 5 13]));
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);
jsthresh = 0.1;
sacwin = 6;
sacthresh = 10;
condition = {'vestibular','visual','combined'};

data = extractEndOfSteeringVars(subject,params,jsthresh,sacwin,sacthresh);
push_success = data.push_success;
push_1st = data.push_1st;
tsac = data.tsac;
letgo = data.letgo;
tau = data.tau;
t_end = data.tend;

b = cellfun(@(x,y,t) regress(y(:)-t(:),[ones(numel(x),1) x(:)-t(:)]),letgo,push_1st,t_end,'un',0);
[rho,pval] = cellfun(@(x,y,t) nancorr(x(:)-t(:),y(:)-t(:)),letgo,push_1st,t_end);

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
a = figure(888); a.Position = [187 554 665 205];
for i = 1:Nsubs
    clf;
    for s = 1:Nstim
        subplot(1,Nstim,s);hold on;
        plot(letgo{i,s}-t_end{i,s},push_1st{i,s}-t_end{i,s},'.','color',colr(s,:)); plot(-20:20,-20:20,'k--');
        xlabel('JS let-go from trial end');ylabel('1st push from trial end');axis equal;axis([-20 20 -20 20]);
        title(condition{s});
    end
    suptitle(subject(i).name);
    pause(1)
end
