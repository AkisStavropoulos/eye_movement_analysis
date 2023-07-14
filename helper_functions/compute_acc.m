function [xaac,yacc,A,cum_sq_acc] = compute_acc(trials)
%% Calculate SQRT cumulative squared acceleration for each trial  
% calculate X,Y components and magnitude of  acc (A)
% WATCH OUT!!!! X IN MC FILE IS FORWARD, BUT IN SPIKE2 IS ANGULAR.
% HERE I SWITCH THE NOTATION TO MATCH THE SPIKE2 DATA WHERE IT IS GOING TO BE SAVED

dt = 1/60;
cum_sq_acc = [];
mc = [trials.mc];
for k = 1:length(mc)
    if ~isnan(mc(k).JS_X_Raw)
        % divide by vmax,wmax to normalize JS input to range between [-1,1]
        JS_X{k} = mc(k).JS_X_Raw/trials(k).prs.vmax;
        JS_Yaw{k} = mc(k).JS_Yaw_Raw/trials(k).prs.wmax;
        A{k} = sqrt(JS_X{k}.^2 + JS_Yaw{k}.^2);

        cum_sq_acc{k} = (cumsum(abs(JS_X{k}).^2 + abs(JS_Yaw{k}).^2));
    else
        A{k} = nan;
        cum_sq_acc{k} = nan;
    end
        % give it the same length as the rest of continuous variables
%     cum_sq_acc{k} = interp1(1:length(mc(k).JS_X_Raw),cum_sq_acc{k},1:length(trials(k).continuous.ts));

end

xaac = JS_Yaw;
yacc = JS_X;