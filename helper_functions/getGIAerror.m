function [GIAerror,GIA,Moog_X_Pos,MaxPos] = getGIAerror(v_input)
%% simulate trajectories through the MC algorithm
% use the GIA error to restrict unwanted parameters for experiment
% simulates straight line trajectories

% MaxPos is in meters (m)

for i = 1:length(v_input)
    v = []; w = [];
    v = v_input{i};
    w = zeros(1,length(v_input{i}));
    v = v/100; % must be in m/s
    time = 1:length(v);
    [~,GIA,GIAerror{i},Moog_X_Pos{i},MaxPos] = Motion_Cueing(time,v,w);
end
