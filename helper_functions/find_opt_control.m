function [T,s,T_trial,s_trial] = find_opt_control(trials,plt)
%% check with optimal trajectory
% plt = 0;
condition = {'vestibular','visual','combined','all trials'};
for i = 1:length(trials)%neg_brake(1:10)%[8:10,408:410]
    if ~ isempty(find(trials(i).mc.JS_X_Raw,1))
        % optimal
        ff_dist(i) = sqrt((trials(i).prs.fireflyposx - trials(i).continuous.xmp(1)).^2 + (trials(i).prs.fireflyposy - trials(i).continuous.ymp(1)).^2);
        T(i) = traveltime(trials(i).prs.tau,ff_dist(i),trials(i).prs.vmax); % optimal trial duration
        s(i) = switchtime(trials(i).prs.tau,T(i)); % optimal switch time
        
        [~,~,~,~,~,dist(i),vel{i},tt(i),~] = ...
            tau2vmax(trials(i).prs.tau, s(i), T(i)-s(i), [], trials(i).prs.x, trials(i).prs.T, plt);
        % actual
        move_indMC(i) = find(trials(i).mc.JS_X_Raw,1); % movement onset indx
        t_offsetMC(i) = trials(i).mc.timestamp(move_indMC(i));

        move_indSMR(i) = find(trials(i).continuous.ts >= t_offsetMC(i),1);
        t_offsetSMR(i) = trials(i).continuous.ts(move_indSMR(i));
        T_trial(i) = trials(i).continuous.ts(end-move_indSMR(i)); % actual trial duration

        [~,s_ind(i)] = max(trials(i).continuous.v);
        s_trial(i) = trials(i).continuous.ts(s_ind(i)-move_indSMR(i)); % actual switchtime
        
        if plt
            hold on;
            plot(trials(i).mc.timestamp(move_indMC(i):end) - t_offsetMC(i),trials(i).mc.JS_X_Raw(move_indMC(i):end));
            plot(trials(i).continuous.ts(move_indSMR(i):end) - t_offsetSMR(i),trials(i).continuous.v(move_indSMR(i):end));
            legend('optimal JS input','optimal velocity','actual JS input','actual velocity')
            suptitle([trials(1).prs.subject ' - ' condition{trials(i).prs.stimtype}]);
        end
    end
end
