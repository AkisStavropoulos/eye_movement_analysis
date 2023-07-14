
%% TTI for different error magnitudes
% Subjects SEPARATELY
%% Modalities
default_prs;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
upperbnd = 0;
BELIEF = 1;
[eye_fixation,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);
%
for i = 1:length(subject)
    figure; hold on;    whitebg([1 1 1]);    cnt = 0;
    for s = 1:size(poolindx,2)
        varexp_start{i,s} = (eye_movement{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.startaligned);  varexp_start{i,s}(varexp_start{i,s}<0) = 0;
        varexp_stop{i,s} = (eye_movement{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.stopaligned);    varexp_stop{i,s}(varexp_stop{i,s}<0) = 0;
        colr = brewermap(size(varexp_start{i,s},1)*size(poolindx,2),'Paired');
        
        T_start = 10; T_stop = 6; T_break = 2;
        dt = 1/60; % NOT SMR SAMPLING RATE, it's is downsampled to 60Hz
        nt_start = round(T_start/dt);
        nt_stop = round(T_stop/dt);
        nt_break = round(T_break/dt);
        if strcmp(params,'stimtau')
        if s<=size(poolindx,2)/3; fig=1; elseif s<=2*size(poolindx,2)/3; fig=2; elseif s<=3*size(poolindx,2)/3; fig=3; end
        else; fig=1; end
        subplot(1,size(poolindx,2)/3,fig); hold on;
        for n = 1:size(varexp_start{i,s},1)
        plot(dt:dt:nt_start*dt,sqrt(varexp_start{i,s}(n,1:nt_start)),'Color',colr(cnt+s+n-1,:));
        plot(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt),flip(sqrt(varexp_stop{i,s}(n,1:nt_stop))),'Color',colr(cnt+s+n-1,:),'HandleVisibility','off');
        end
        xlabel('Time (s)'); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
        cnt = cnt + 1;
    end
    title('modalities');    suptitle(subject(i).name);
end

% All trials
params = []; upperbnd = 0; BELIEF = 1;
[poolindx,~] = get_poolindx(subject,params);
prs.ngroups = 5; prs.boots = 0;
[eye_fixationALL,eye_movementALL] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);

for i = 1:length(subject)
    figure; hold on;    whitebg([1 1 1]);
    for s = 1:size(poolindx,2)
        varexp_start{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.startaligned);  varexp_start{i,s}(varexp_start{i,s}<0) = 0;
        varexp_stop{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.stopaligned);    varexp_stop{i,s}(varexp_stop{i,s}<0) = 0;
        colr = jet(size(varexp_start{i,s},1));
        
        T_start = 10; T_stop = 6; T_break = 2;
        dt = 1/60; % NOT SMR SAMPLING RATE, it's is downsampled to 60Hz
        nt_start = round(T_start/dt);
        nt_stop = round(T_stop/dt);
        nt_break = round(T_break/dt);
        for n = 1:size(varexp_start{i,s},1)
        plot(dt:dt:nt_start*dt,sqrt(varexp_start{i,s}(n,1:nt_start)),'Color',colr(n,:));
        plot(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt),flip(sqrt(varexp_stop{i,s}(n,1:nt_stop))),'Color',colr(n,:),'HandleVisibility','off');
        end
        xlabel('Time (s)'); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
        
    end
    suptitle(subject(i).name);
end

% All subjects
% Modalities
% nothing here (not even for modalities or tau)
default_prs;
params = 'tau';
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
upperbnd = 0; BELIEF = 1; prs.boots = 0;
[eye_fixation,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);

varexp_start = [];  varexp_stop = [];   startLNG = [];  stopLNG = [];
for i = 1:length(subject)
    for s = 1:size(poolindx,2)
        varexp_start{i,s} = (eye_movement{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.startaligned);  varexp_start{i,s}(varexp_start{i,s}<0) = 0;
        startLNG(i,s) = length(varexp_start{i,s});
        varexp_stop{i,s} = (eye_movement{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.stopaligned);    varexp_stop{i,s}(varexp_stop{i,s}<0) = 0;
        stopLNG(i,s) = length(varexp_stop{i,s});
    end
end
% Remove bad eye-tracking subjects
i = [ 5 ]; % bad eye-tracking subjects
varexp_start(i,:) = []; varexp_stop(i,:) = [];
startLNG(i,:) = [];     stopLNG(i,:) = [];

for s = 1:size(poolindx,2)
    figure; hold on;    whitebg([1 1 1]);         
    MINindxSTART = min(startLNG(:,s));
    MINindxSTOP = min(stopLNG(:,s));
    tempSTART = []; tempSTOP = [];      
    for i = 1:size(startLNG,1)
        for n = 1:size(varexp_start{i,s},1)
            tempSTART(i,n,:) = varexp_start{i,s}(n,1:MINindxSTART);
            tempSTOP(i,n,:) = varexp_stop{i,s}(n,1:MINindxSTOP);
        end
    end
    colr = jet(size(varexp_start{i,s},1));
    startAVG{s} = nanmean(tempSTART,1);     startSEM{s} = nanstd(tempSTART,1)./sqrt(length(subject));
    stopAVG{s} = nanmean(tempSTOP,1);       stopSEM{s} = nanstd(tempSTOP,1)./sqrt(length(subject));
    
    for n = 1:size(varexp_start{i,s},1)
    % plot
    T_start = 6; T_stop = 6; T_break = 2;
    dt = 1/60; % NOT SMR SAMPLING RATE, it's is downsampled to 60Hz
    nt_start = round(T_start/dt);
    nt_stop = round(T_stop/dt);
    nt_break = round(T_break/dt);
    
    shadedErrorBar(dt:dt:nt_start*dt,sqrt(startAVG{s}(1,n,1:nt_start)),startSEM{s}(1,n,1:nt_start),'lineprops',{'Color',colr(n,:)});
    shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt),flip(sqrt(stopAVG{s}(1,n,1:nt_stop))),flip(stopSEM{s}(1,n,1:nt_stop)),'lineprops',{'Color',colr(n,:)});
    xlabel('Time (s)'); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
    end
    legend(legend_input{s});
end

% All trials
% nothing here (not even for modalities or tau)
default_prs;
params = [];
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
upperbnd = 0; BELIEF = 1; prs.boots = 0; prs.ngroups = 5;
[eye_fixationALL,eye_movementALL] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF);

varexp_start = [];  varexp_stop = [];   startLNG = [];  stopLNG = [];
for i = 1:length(subject)
    for s = 1:size(poolindx,2)
        varexp_start{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.startaligned);  varexp_start{i,s}(varexp_start{i,s}<0) = 0;
        startLNG(i,s) = length(varexp_start{i,s});
        varexp_stop{i,s} = (eye_movementALL{i,s}.eyepos.pred_vs_true.var_explained_vs_error.mu.stopaligned);    varexp_stop{i,s}(varexp_stop{i,s}<0) = 0;
        stopLNG(i,s) = length(varexp_stop{i,s});
    end
end
% Remove bad eye-tracking subjects
i = [ 5 ]; % bad eye-tracking subjects
varexp_start(i,:) = []; varexp_stop(i,:) = [];
startLNG(i,:) = [];     stopLNG(i,:) = [];

for s = 1:size(poolindx,2)
    figure; hold on;    whitebg([1 1 1]);         
    MINindxSTART = min(startLNG(:,s));
    MINindxSTOP = min(stopLNG(:,s));
    tempSTART = []; tempSTOP = [];      
    for i = 1:size(startLNG,1)
        for n = 1:size(varexp_start{i,s},1)
            tempSTART(i,n,:) = varexp_start{i,s}(n,1:MINindxSTART);
            tempSTOP(i,n,:) = varexp_stop{i,s}(n,1:MINindxSTOP);
        end
    end
    colr = jet(size(varexp_start{i,s},1));
    startAVG{s} = nanmean(tempSTART,1);     startSEM{s} = nanstd(tempSTART,1)./sqrt(length(subject));
    stopAVG{s} = nanmean(tempSTOP,1);       stopSEM{s} = nanstd(tempSTOP,1)./sqrt(length(subject));
    
    for n = 1:size(varexp_start{i,s},1)
    % plot
    T_start = 6; T_stop = 6; T_break = 2;
    dt = 1/60; % NOT SMR SAMPLING RATE, it's is downsampled to 60Hz
    nt_start = round(T_start/dt);
    nt_stop = round(T_stop/dt);
    nt_break = round(T_break/dt);
    
    shadedErrorBar(dt:dt:nt_start*dt,sqrt(startAVG{s}(1,n,1:nt_start)),startSEM{s}(1,n,1:nt_start),'lineprops',{'Color',colr(n,:)});
    shadedErrorBar(nt_start*dt+ nt_break*dt + (dt:dt:nt_stop*dt),flip(sqrt(stopAVG{s}(1,n,1:nt_stop))),flip(stopSEM{s}(1,n,1:nt_stop)),'lineprops',{'Color',colr(n,:)});
    xlabel('Time (s)'); ylabel('Target-tracking index'); axis([0 T_start+T_stop+T_break 0 1]);
    end
    legend(legend_input{s});
end
