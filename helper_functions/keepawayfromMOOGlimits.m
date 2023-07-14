function [indx_allow,Moog_X_Pos] = keepawayfromMOOGlimits(Moog_X_Pos,MaxPos,thresh,plots)

%% keep simulations that don't hit the Moog limits
% MaxPos must be in meters (m)
% plots: 0 no, 1 yes

indx_allow = []; % in terms of inputed struct!!!! (e.g. infostruct, vardist)
for i = 1:length(Moog_X_Pos)
    limit_hit = find(abs(MaxPos - Moog_X_Pos{i}) < thresh,1);
    if isempty(limit_hit)
        indx_allow = [indx_allow i];
    end
    limit_hit = [];
end


% plot
% give all vectors same length to plot
for i = 1:length(Moog_X_Pos)
    maxlength(i) = max(length(Moog_X_Pos{i}));
end
maxlength = max(maxlength);
for i = 1:length(Moog_X_Pos)
    dif = maxlength - length(Moog_X_Pos{i});
    if  dif
        Moog_X_Pos{i} = [Moog_X_Pos{i} ; zeros(dif,1)];
    end
end

dt = 1/60;
ts = dt*(1:maxlength);
numsims = 1:length(indx_allow);
if plots
    h(1) = figure;
    figure(h(1));surfc(numsims,ts,100*[Moog_X_Pos{indx_allow}]);xlabel('simulations of increasing \tau');title(['Simulations that are > ' num2str(thresh*100) ' cm away from Moog limits']);
    zlabel('Moog Position (cm)');ylabel('time (s)');

end