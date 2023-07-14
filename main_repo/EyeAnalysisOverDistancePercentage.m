function tracking_distperc = EyeAnalysisOverDistancePercentage(tracking_regular,dd,align_flag,prs)
%% Re-sample quantities based on percentage of distance traveled
% align_flag: either 'align2targ', or 'align2stop', or 'alignall'.

disp('Re-sample Eye Data based on % of distance traveled..........')

Nsubs = size(tracking_regular,1);
Nstim = size(tracking_regular,2);

for i = 1:Nsubs
    for s = 1:Nstim
        tracking_distperc{i,s} = EyeDataOverDistancePerc(tracking_regular{i,s},dd,prs,align_flag);
    end
    disp(['..........Subject = ' num2str(i)])
end
