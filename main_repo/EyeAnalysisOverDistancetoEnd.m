function tracking_dist2end = EyeAnalysisOverDistancetoEnd(tracking_regular,dd,align_flag,prs)

%% Re-sample quantities based on distance to the end (end aligned)

disp('Re-sample Eye Data based on distance to the end..........')

Nsubs = size(tracking_regular,1);
Nstim = size(tracking_regular,2);

for i = 1:Nsubs
    for s = 1:Nstim
        tracking_dist2end{i,s} = EyeDataOverDistanceEndAligned(tracking_regular{i,s},dd,prs,align_flag);
    end
    disp(['..........Subject = ' num2str(i)])
end
