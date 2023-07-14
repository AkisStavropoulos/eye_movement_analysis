%% mean velocity and duration of trials
%% mean velocities
velocity = [trials.continuous];
fields = fieldnames(velocity);
velocity = rmfield(velocity,fields(~strcmp(fields,'v')));

mean_v = [];
for i = 1:length(velocity)
mean_v(i) = mean(velocity(i).v);
end
total_mean_v = mean(mean_v);

% across joystick controls
clear mean_v;clear std_v;
for n = 1:length(coefnames)-1
    indx = jscoef.(coefnames{n});
    for i = 1:length(indx)
        mean_v.(coefnames{n})(i) = mean(velocity(indx(i)).v);
    end
    clear indx;
    std_v.(coefnames{n}) = std(mean_v.(coefnames{n}));
    mean_v.(coefnames{n}) = mean(mean_v.(coefnames{n}));
end

% across stimulus types
for n = 1:length(stnames)-1
    indx = stimtype.(stnames{n});
    for i = 1:length(indx)
        mean_v.(stnames{n})(i) = mean(velocity(indx(i)).v);
    end
    clear indx;
    std_v.(stnames{n}) = std(mean_v.(stnames{n}));
    mean_v.(stnames{n}) = mean(mean_v.(stnames{n}));
end

% across all conditions
for n = 1:length(condnames)-1
    indx = stimjs.(condnames{n});
    for i = 1:length(indx)
        mean_v.(condnames{n})(i) = mean(velocity(indx(i)).v);
    end
    clear indx;
    std_v.(condnames{n}) = std(mean_v.(condnames{n}));
    mean_v.(condnames{n}) = mean(mean_v.(condnames{n}));
end

mean_v = orderfields(mean_v,std_v);
%% mean duration

mean_dur = [];
for i = 1:length(velocity)
    dur(i) = length(velocity(i).v)*prs.dt;
end
total_mean_dur = mean(dur);

% across joystick controls
clear mean_dur;clear std_dur;
for n = 1:length(coefnames)-1
    indx = jscoef.(coefnames{n});
    for i = 1:length(indx)
        mean_dur.(coefnames{n})(i) = length(velocity(indx(i)).v)*prs.dt;
    end
    clear indx;
    mean_dur.(coefnames{n}) = mean(mean_dur.(coefnames{n}));
end

% across stimulus types
for n = 1:length(stnames)-1
    indx = stimtype.(stnames{n});
    for i = 1:length(indx)
        mean_dur.(stnames{n})(i) = length(velocity(indx(i)).v)*prs.dt;
    end
    clear indx;
    mean_dur.(stnames{n}) = mean(mean_dur.(stnames{n}));
end

% across all conditions
for n = 1:length(condnames)-1
    indx = stimjs.(condnames{n});
    for i = 1:length(indx)
        mean_dur.(condnames{n})(i) = length(velocity(indx(i)).v)*prs.dt;
    end
    clear indx;
    mean_dur.(condnames{n}) = mean(mean_dur.(condnames{n}));
end

%% mean duration normalized by distance
        
norm_dur = [];
for i = 1:length(velocity)
    dur(i) = length(velocity(i).v)*prs.dt;
end
total_mean_dur = mean(dur);

% across joystick controls
clear norm_dur;
for n = 1:length(coefnames)-1
    indx = jscoef.(coefnames{n});
    for i = 1:length(indx)
        dist = target_dist(trials,indx(i));
        norm_dur.(coefnames{n})(i) = length(velocity(indx(i)).v)*prs.dt/dist;
    end
    clear indx;
    norm_dur.(coefnames{n}) = 100*mean(norm_dur.(coefnames{n}));
end

% across stimulus types
for n = 1:length(stnames)-1
    indx = stimtype.(stnames{n});
    for i = 1:length(indx)
        dist = target_dist(trials,indx(i));
        norm_dur.(stnames{n})(i) = length(velocity(indx(i)).v)*prs.dt/dist;
    end
    clear indx;
    norm_dur.(stnames{n}) = 100*mean(norm_dur.(stnames{n}));
end

% across all conditions
for n = 1:length(condnames)-1
    indx = stimjs.(condnames{n});
    for i = 1:length(indx)
        dist = target_dist(trials,indx(i));
        norm_dur.(condnames{n})(i) = length(velocity(indx(i)).v)*prs.dt/dist;
    end
    clear indx;
    norm_dur.(condnames{n}) = 100*mean(norm_dur.(condnames{n}));
end
