function data = extractEndOfSteeringVars(subject,params,jsthresh,sacwin,sacthresh)

%% extract end-of-trial events and behaviour
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = numel(subject);
Nstim = size(poolindx,2);

wtf_time = 20;
presac_t = 0.2;
postsac_t = presac_t;
velthresh = 10;
avgwin = 80;

dt = diff(subject(1).trials(1).continuous.ts(1:2)); 

for i = 1:Nsubs
    
    % find blocks that button was not recorded at all
    Ntrls = numel(subject(i).trials);
    norec_trls = nan(Ntrls,1);
    filename = unique(arrayfun(@(x) x.prs.filename,subject(i).trials,'un',0));
    for k = 1:numel(filename)
        fileindx = arrayfun(@(x) strcmp(x.prs.filename,filename{k}), subject(i).trials);
        push_trls = arrayfun(@(x) numel(x.events.push_on),subject(i).trials(fileindx)) > 1;
        norec = sum(push_trls) == 0;
        if norec; norec_trls(fileindx) = true; else; norec_trls(fileindx) = false; end
    end
    
    for s = 1:Nstim
        trlind = poolindx{i,s};
        
        tau = []; tsac_mu = []; t_beg = []; t_end = []; push_success = []; push_1st = []; letgo_t = [];
        for j = 1:numel(trlind)
        
        ts = subject(i).trials(trlind(j)).continuous.ts;
        ze = 0.5*(subject(i).trials(trlind(j)).continuous.zle + subject(i).trials(trlind(j)).continuous.zre);
        t_beg(j) = subject(i).trials(trlind(j)).events.t_beg;
        t_end(j) = subject(i).trials(trlind(j)).events.t_end - t_beg(j); 
        % succesful push
        try push_success(j) = subject(i).trials(trlind(j)).events.button - t_beg(j); catch; push_success(j) = nan; end
        % first push
        try     push_1st(j) = subject(i).trials(trlind(j)).events.push_on(1); 
        catch;  if 1%norec_trls(trlind(j)) 
                   push_1st(j) = nan; 
                else
                   push_1st(j) = push_success(j); % add successful push if it was the first push
                end
        end        
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
        rmindx = t_sac > ts(end); 
        t_sac(rmindx)=[]; sac_mag(rmindx)=[];
        ver_eye = 0.5*(subject(i).trials(trlind(j)).continuous.zle + subject(i).trials(trlind(j)).continuous.zre);
        pre_sac = arrayfun(@(x) find(ts >= (x - presac_t),1,'first'),t_sac);
        post_sac = arrayfun(@(x) find(ts <= (x + postsac_t),1,'last'),t_sac);
        sac_ver_dir = arrayfun(@(pre,post) ver_eye(post)-ver_eye(pre),pre_sac,post_sac);
        % max of upward final saccades [0-sacwin] sec before end of trial
        try     
%                 t_indx = t_sac > t_end(j) - sacwin; % t_indx = t_sac > t_end(j) - sacwin; 
%                 if sacwin <= 1; t_indx = t_sac > t_end(j) - sacwin*t_end(j); end % input is percentage of trial
%                 sac_mag_tmp = sac_mag(t_indx);
%                 sac_dir_tmp = sac_ver_dir(t_indx);
%                 bigsac = sac_mag_tmp(find(sac_mag_tmp > sacthresh & sac_dir_tmp > 0, 1)); % first big upward saccade within [sacwin] from trial end
% %                 bigsac = max(sac_mag_tmp); % biggest saccade within [sacwin] from trial end
%                 tsac_tmp = t_sac(sac_mag==bigsac);
%                 if sacwin <= 1; tsac_mu(j) = nanmean(tsac_tmp(tsac_tmp > t_end(j) - sacwin*t_end(j)));
%                 else; tsac_mu(j) = nanmean(tsac_tmp(tsac_tmp > t_end(j) - sacwin));
%                 end
                
                % alternative detection using vertical velocity threshold
                t_indx = diff(movmean(ze,avgwin))/dt > velthresh;
                tsac_tmp = ts(t_indx);
                if sacwin <= 1; tsac_mu(j) = tsac_tmp(find(tsac_tmp > t_end(j) - sacwin*t_end(j),1)); 
                else; tsac_mu(j) = tsac_tmp(find(tsac_tmp > t_end(j) - sacwin,1)); 
                end
        catch;  tsac_mu(j) = nan;
        end
                
        end
        t_end(t_end>wtf_time)=nan;
        push_success(push_success>wtf_time)=nan;
        push_1st(push_1st>wtf_time)=nan;
        letgo_t(letgo_t>wtf_time)=nan;
        tsac_mu(tsac_mu>wtf_time)=nan;
        
        % assign variables
        data.tau{i,s} = tau;
        data.tsac{i,s} = tsac_mu;
        data.push_success{i,s} = push_success;
        data.push_1st{i,s} = push_1st;
        data.letgo{i,s} = letgo_t;
        data.tend{i,s} = t_end;
    end
end