function trials = position_theta_fix(trials)
%% if Spike 2 position is faulty, generate trajectory from velocity and angle
for i = 1:length(trials)
[trials(i).continuous.xmp, trials(i).continuous.ymp,~, ~, ~, trials(i).continuous.phi] = ...
    gen_traj(trials(i).continuous.w, trials(i).continuous.v, trials(i).continuous.xmp(1), trials(i).continuous.ymp(1));
end
