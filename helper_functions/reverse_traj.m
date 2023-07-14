function subject = reverse_traj(subject,trialtaus,plts)
%% for every trial, find tangent of last part of trajectory to align based on target
% Refer to TestRotation script for the method
% trialtaus: index based on a trial_t/tau ratio (how many taus from the end of the trial)
%% Alternatively
% traj_perc: choose how much of the last part of the trajectory you want to get the average tangent
% suggest traj_perc = .65
% e.g.:
% indx = ceil(trialtaus*numel(subject(i).trials(j).continuous.ts)); % tangent of last .65 of trajectory
%%
% plts: 1 for plots 0 for no plots

% you can index based on the ratio of distance covered to target distance instead of time

for i = 1:length(subject)
    reverse_traj_x = [];
    reverse_traj_y = [];
    avg_phi = [];
    
    if plts
        figure;
        xlim([-400 400]);ylim([-700 100]);hline(0);vline(0);grid on;
        title([subject(i).name]);
    end
    for j = 1:length(subject(i).trials)
        phi = [];
        R = [];
        S1 = [];
        S2 = [];
        tau = subject(i).trials(j).prs.tau;
        t_trial = subject(i).trials(j).continuous.ts(end);
        indx = find(subject(i).trials(j).continuous.ts >= t_trial - tau*trialtaus, 1);

        % subtract endpoint from trajectory
        r_in = sqrt((subject(i).trials(j).continuous.xmp(end) - subject(i).trials(j).continuous.xmp(1)).^2 + (subject(i).trials(j).continuous.ymp(end) - subject(i).trials(j).continuous.ymp(1)).^2);
        
        reverse_traj_x{j} = subject(i).trials(j).continuous.xmp - subject(i).trials(j).continuous.xmp(end);
        reverse_traj_y{j} = subject(i).trials(j).continuous.ymp - subject(i).trials(j).continuous.ymp(end);
        
        r_sub = sqrt((reverse_traj_x{j}(end) - reverse_traj_x{j}(1)).^2 + (reverse_traj_y{j}(end) - reverse_traj_y{j}(1)).^2);
        
        phi = subject(i).trials(j).continuous.phi(indx:end);
        avg_phi(j) = mean(phi);
        % rotation matrix
        R = [cosd(avg_phi(j)) sind(avg_phi(j)) ;...
            -sind(avg_phi(j)) cosd(avg_phi(j))];
        S1 = [reverse_traj_x{j} reverse_traj_y{j}];
        
        S2 = S1*R;
        reverse_traj_x{j} = S2(:,1);
        reverse_traj_y{j} = S2(:,2);
        
        r_rev = sqrt((reverse_traj_x{j}(end) - reverse_traj_x{j}(1)).^2 + (reverse_traj_y{j}(end) - reverse_traj_y{j}(1)).^2);
        % save reversed trajectories in a struct
        subject(i).reversed_traj(j).xmp = flip(reverse_traj_x{j});
        subject(i).reversed_traj(j).ymp = flip(reverse_traj_y{j});
        
        if plts
            hold on;plot(reverse_traj_x{j},reverse_traj_y{j},'k');
        end
    end
    
end
