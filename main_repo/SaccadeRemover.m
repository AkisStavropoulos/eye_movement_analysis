function subject = SaccadeRemover(subject,prs)
%% Remove saccades from all trials and keep only smooth component

Nsubs = numel(subject);
dt = subject(1).trials(1).prs.dt;
targetoff = 60;
velthresh = prs.velthresh;

for i = 1:Nsubs
    ntrials = numel(subject(i).trials);
    for j = 1:ntrials
        zle = subject(i).trials(j).continuous.zle;
        yle = subject(i).trials(j).continuous.yle;
        zre = subject(i).trials(j).continuous.zre;
        yre = subject(i).trials(j).continuous.yre;
        
        ver_mean = nanmean([zle , zre],2); % mean vertical eye position (of the two eyes)
        hor_mean = nanmean([yle , yre],2); % mean horizontal eye position
        
        % remove saccades and concatenate eye position
        % both eyes
        startpos.ver = ver_mean(targetoff);
        startpos.hor = hor_mean(targetoff);        
        eyevel.ver = diff(ver_mean)/dt;
        eyevel.hor = diff(hor_mean)/dt;
        eyevel.mag = sqrt(eyevel.ver.^2 + eyevel.hor.^2);
        
        nanindx = unique([find(abs(eyevel.ver) > velthresh) ; find(abs(eyevel.hor) > velthresh) ; find(abs(eyevel.mag) > velthresh)]);
        eyevel.ver(nanindx) = nan;
        eyevel.hor(nanindx) = nan;
        pre_ver_mean = ver_mean;
        pre_hor_mean = hor_mean;
        ver_mean = [0 ; startpos.ver + cumsum(eyevel.ver,'omitnan')*dt];
        hor_mean = [0 ; startpos.hor + cumsum(eyevel.hor,'omitnan')*dt];
        
        % each eye separately
        startpos.zle = zle(targetoff);        eyevel.zle = diff(zle)/dt;
        startpos.zre = zre(targetoff);        eyevel.zre = diff(zre)/dt;
        startpos.yle = yle(targetoff);        eyevel.yle = diff(yle)/dt;
        startpos.yre = yre(targetoff);        eyevel.yre = diff(yre)/dt;
        
        eyevel.zle(nanindx) = nan;
        eyevel.zre(nanindx) = nan;
        eyevel.yle(nanindx) = nan;
        eyevel.yre(nanindx) = nan;
        pre_zle = zle;
        pre_zre = zre;
        pre_yle = yle;
        pre_yre = yre;
        zle = [0 ; startpos.zle + cumsum(eyevel.zle,'omitnan')*dt];
        zre = [0 ; startpos.zre + cumsum(eyevel.zre,'omitnan')*dt];
        yle = [0 ; startpos.yle + cumsum(eyevel.yle,'omitnan')*dt];
        yre = [0 ; startpos.yre + cumsum(eyevel.yre,'omitnan')*dt];
        
        % replace eye positions
        subject(i).trials(j).continuous.zle = zle;
        subject(i).trials(j).continuous.yle = yle;
        subject(i).trials(j).continuous.zre = zre;
        subject(i).trials(j).continuous.yre = yre;

    end
end
