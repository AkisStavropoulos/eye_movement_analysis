function trials = save_errors(trials,r_tar,r_sub,theta_tar,theta_sub)
%     save target distance/angle, errors in trials struct
% signed error
for i = 1:length(trials)
    trials(i).prs.r_tar = r_tar(i);
    trials(i).prs.r_sub = r_sub(i);
    trials(i).prs.th_tar = theta_tar(i);
    trials(i).prs.th_sub = theta_sub(i);
    
   trials(i).stats.d_err = r_sub(i) - r_tar(i);
   
   % make sure you get under-/over-shooting correct
   if theta_tar(i) >= 0
   trials(i).stats.th_err = theta_sub(i) - theta_tar(i);
   elseif theta_tar(i) < 0
   trials(i).stats.th_err = -(theta_sub(i) - theta_tar(i));
   end
   
end
