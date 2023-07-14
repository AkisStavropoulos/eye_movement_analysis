function trials = save_errors2(trials,r_sub2,th_sub2)
%     save target distance/angle, errors in trials struct
% signed error
for i = 1:length(trials)
    trials(i).prs.r_sub2 = r_sub2(i);
    trials(i).prs.th_sub2 = th_sub2(i);
    
   trials(i).stats.d_err2 = r_sub2(i) - trials(i).prs.r_tar;
   
   if trials(i).prs.th_tar >= 0
   trials(i).stats.th_err2 = th_sub2(i) - trials(i).prs.th_tar;
   elseif trials(i).prs.th_tar < 0
   trials(i).stats.th_err2 = -(th_sub2(i) - trials(i).prs.th_tar);
   end
   
end
