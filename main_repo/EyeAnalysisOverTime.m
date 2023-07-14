function tracking_regular = EyeAnalysisOverTime(subject,params,rmsac,prs)

%% Compare tracking of target position versus stopping position
[poolindx,legend_input] = get_poolindx(subject,params);
Nsubs = length(subject);
Nstim = size(poolindx,2);

disp('Eye Analysis over time..........')
for i = 1:Nsubs
    for s = 1:Nstim
        indx = poolindx{i,s};
        tracking_regular{i,s} = TargetVsStopPosition_tracking(subject(i).trials(indx),rmsac,prs);
    end
    disp(['........Subject = ' num2str(i)])
end
