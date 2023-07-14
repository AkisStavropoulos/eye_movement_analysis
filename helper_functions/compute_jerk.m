function [Jx,Jy,J,cum_sq_jerk] = compute_jerk(trials)
%% Calculate SQRT cumulative squared jerk for each trial
% calculate X,Y components and magnitude of jerk (J)
% WATCH OUT!!!! X IN MC FILE IS FORWARD, BUT IN SPIKE2 IS ANGULAR.
% HERE I SWITCH THE NOTATION TO MATCH THE SPIKE2 DATA WHERE IT IS GOING TO BE SAVED
dt = 1/60;
cum_sq_jerk = [];
Jx = [];
Jy = [];
J = [];
mc = [trials.mc];
% figure;hold on;
for k = 1:length(mc)
    if ~isnan(mc(k).JS_X_Raw)
        % divide by vmax,wmax to normalize JS input to range between [-1,1]        
        JS_X{k} = mc(k).JS_X_Raw/trials(k).prs.vmax;
        JS_Yaw{k} = mc(k).JS_Yaw_Raw/trials(k).prs.wmax;
            
        Jy{k} = [0 ; diff(JS_X{k})./dt];
        Jx{k} = [0 ; diff(JS_Yaw{k})./dt];
        
        J{k} = sqrt(Jx{k}.^2 + Jy{k}.^2);
        cum_sq_jerk{k} = cumsum(sqrt(Jx{k}.^2 + Jy{k}.^2)); % Xaq said: cumsum(sqrt(Jx{k}.^2 + Jy{k}.^2)) to avoid jumps
        
        % avoid jump at first timepoint
        indx = find(cum_sq_jerk{k} > 0 ,1);
        if isempty(indx)
            cum_sq_jerk{k} = nan;
            J{k} = nan;
            Jx{k} = nan;
            Jy{k} = nan;
        else
           if 0; cum_sq_jerk{k} = [cum_sq_jerk{k}(1:indx-1) ; (cum_sq_jerk{k}(indx:end) - cum_sq_jerk{k}(indx))]; end
%             Jx{k} = [Jx{k}(1:indx-1) ; (Jx{k}(indx:end) - Jx{k}(indx))];
%             Jy{k} = [Jy{k}(1:indx-1) ; (Jy{k}(indx:end) - Jy{k}(indx))];
        end
        
    else
        cum_sq_jerk{k} = nan;
        J{k} = nan;
        Jx{k} = nan;
        Jy{k} = nan;
    end
    
    
    %     plot(cum_sq_jerk{k});
% % check jerk profiles    
%     rmend = flip(cum_sq_jerk{k});
%     indx = find(abs(diff(rmend)) > 0,1);
%     rmend(1:indx) = nan;
%     plot(flip(rmend));
% give it the same length as the rest of continuous variables
%     cum_sq_jerk{k} = interp1(1:length(mc(k).JS_X_Raw),cum_sq_jerk{k},1:length(trials(k).continuous.ts));
end
% xlim([0 2500]);ylim([0 1200]);hline(0);
% title(trials(1).prs.subject);
