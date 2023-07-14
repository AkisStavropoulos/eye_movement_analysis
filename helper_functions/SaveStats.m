function [subject,button] = SaveStats(subject,data_folder,savenow)

%% SaveStats
% save basic computed quantities for each subject
% any quantities you want to save, add them here

%% TAU binning properties
limtype = 'hard';
bins = input('Number of TAU bins? ');
button = input('Use position at 1st button push (1) or at end of trial (0)? ');
%% check distribution of taus for proper binning
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        taus(count) = subject(i).trials(j).prs.tau;
    end
end
[tau_cdf,taus_x] = ecdf(taus);
% figure;ecdf(taus);title('CDF of all taus from all trials');xlabel('\tau [s]');
%%
for i = 1:length(subject)
        % fix sign of w
        for j = 1:length(subject(i).trials)
            if sign(subject(i).trials(j).continuous.phi(end)) ~= sign(sum(subject(i).trials(j).continuous.w))
           subject(i).trials(j).continuous.w =  -subject(i).trials(j).continuous.w;
            end
        end
        % sort trials
        [subject(i).stimtype,subject(i).tau,subject(i).stimtau] = ...
            condsorting_bintau(subject(i).trials,bins,limtype,tau_cdf,taus_x);
        
        % save distance and angular error for every trial
        [r_tar,r_sub,theta_tar,theta_sub] = scatterDistAng(subject(i).trials,button,0,0);
        subject(i).trials = save_errors(subject(i).trials,r_tar,r_sub,theta_tar,theta_sub);
        clear('r_tar','r_sub','theta_tar','theta_sub')
        % save euclideian error
        eucl_err = compute_eucl_err(subject(i).trials);
        subject(i).trials = save_eucl_errors(subject(i).trials,eucl_err);
        clear('eucl_err')
        % save cum squared jerk for every trial
%         [xjerk,yjerk,jerk,cumjerk] = compute_jerk(subject(i).trials);
%         subject(i).trials = save_jerk(subject(i).trials,xjerk,yjerk,jerk,cumjerk);
        clear('xjerk','yjerk','cumjerk')
        % save cum sq acceleration
%         [xacc,yacc,acc,cumacc] = compute_acc(subject(i).trials);
%         subject(i).trials = save_acc(subject(i).trials,xacc,yacc,acc,cumacc);
        clear('xacc','yacc','cumacc')
        % Save SPC
        STOPVEL = 3;
        [~,~,stop_phase,pred_coast] = compare2predcoast(subject(i),STOPVEL,0);
        SPC = stop_phase{:}./pred_coast{:};
%         for j = 1:length(subject(i).trials)
%             subject(i).trials(j).stats.SPC = SPC(j);
%         end        
        % save stopping behavior
        STOPVEL = 3;
        [brake,fwdcorr,overbrake] = extract_stop_behavior(subject(i),STOPVEL,stop_phase,pred_coast,0);
        for j = 1:length(subject(i).trials)
            subject(i).trials(j).stats.brake = brake(j);
            subject(i).trials(j).stats.fwdcorr = fwdcorr(j);
            subject(i).trials(j).stats.overbrake = overbrake(j);
        end        
        clear('stop_phase','pred_coast','SPC','brake','overbrake','fwdcorr')
        if 0
        % save accuracy info (estimated normalized distance collapsed over all distance bins)
        binspace = 50;
        [dist_bins] = bin_dist(subject(i),binspace);
        params = [];
        [poolindx,~] = get_poolindx(subject(i),params);
        if i == 13 
            keyboard
        end
        [x_hat_avg,x_hat_std,x_avg,x_std] = var_given_dist(dist_bins,poolindx);
        
        for n = 1:length(dist_bins.indx)
            for j = 1:length(dist_bins.indx{n})
                subject(i).trials(dist_bins.indx{n}(j)).stats.d_hat_avg = x_hat_avg{n};
                subject(i).trials(dist_bins.indx{n}(j)).stats.d_hat_std = x_hat_std{n};
                subject(i).trials(dist_bins.indx{n}(j)).stats.d_avg = x_avg{n};
            end
        end
        clear('poolindx','dist_bins','x_hat_avg','x_hat_std','x_avg','x_std')
        end
        %% Condition sorting fix
%         taunames = fieldnames(subject(i).tau);
%         lims = [];
%         bins = [];
%         for n = 1:length(taunames)
%             if strcmp(taunames{n}(1),'l')
%                 subject(i).tau.lims.(taunames{n}) = subject(i).tau.(taunames{n});
%                 subject(i).tau = rmfield(subject(i).tau,taunames{n});
%                 
%                 lims = [lims n];
%             elseif strcmp(taunames{n}(1),'b')
%                 subject(i).tau.bins.(taunames{n}) = subject(i).tau.(taunames{n});
%                 subject(i).tau = rmfield(subject(i).tau,taunames{n});
%                 
%                 bins = [bins n];
%             end
%         end
        if savenow
        %% Save permanently on .mat files
        trials = [subject(i).trials];
        stimtype = [subject(i).stimtype];
        tau = [subject(i).tau];
        stimtau = [subject(i).stimtau];
        name = [subject(i).name];
        
        cd([data_folder subject(i).name]);
        save(subject(i).name,'trials','stimtype','tau','stimtau','name');
%         save(subject(i).name,'trials');
        disp(['Data saved: ' subject(i).name]);
        clear('trials','stimtype','tau','stimtau','name')
        end
end