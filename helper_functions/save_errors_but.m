function trials = save_errors_but(trials,r_sub_but,th_sub_but,v_but)
%     save target distance/angle, errors in trials struct
% signed error
for i = 1:length(trials)
    trials(i).prs.r_sub_but = r_sub_but(i);
    trials(i).prs.th_sub_but = th_sub_but(i);
    trials(i).prs.v_but = v_but(i);
        
   trials(i).stats.d_err_but = r_sub_but(i) - trials(i).prs.r_tar;
   
   if trials(i).prs.th_tar >= 0
   trials(i).stats.th_err_but = th_sub_but(i) - trials(i).prs.th_tar;
   elseif trials(i).prs.th_tar < 0
   trials(i).stats.th_err_but = -(th_sub_but(i) - trials(i).prs.th_tar);
   end
   
end
