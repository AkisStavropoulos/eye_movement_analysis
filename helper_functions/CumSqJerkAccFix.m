%% Calculate and save cumulative squared jerk and acceleration for each trial
dt = 1/60;
for i = 1:length(subject)
    mc = [];
    jerk = [];
    acc = [];
    mc = [subject(i).trials.mc];
    for k = 1:length(mc)
        if ~isnan(mc(k).JS_X_Raw)
        jerk{k} = cumsum((diff(mc(k).JS_X_Raw)./dt).^2 + (diff(mc(k).JS_Yaw_Raw)./dt).^2);
        acc{k} = cumsum((mc(k).JS_X_Raw./dt).^2 + (mc(k).JS_Yaw_Raw./dt).^2);
        else
            jerk{k} = nan;
            acc{k} = nan;
        end
        subject(i).trials(k).continuous.jerk = [0 ; jerk{k}];
        subject(i).trials(k).continuous.acc = acc{k};
    end   
end

