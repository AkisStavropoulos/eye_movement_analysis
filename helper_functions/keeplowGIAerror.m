function [indx_allow,GIAerror] = keeplowGIAerror(GIAerror,thresh,plots)

%% keep simulations with low GIA error
% threshold must be in m/s^2
% input threshold in cm/s^2 and transform
% plots: 0 no, 1 yes
thresh = thresh/100;
indx_allow = []; % in terms of inputed struct!!!! (e.g. infostruct, vardist)
for i = 1:length(GIAerror)
    big_err = find(abs(GIAerror{i}) > thresh);
    if isempty(big_err)
        indx_allow = [indx_allow i];
    end
end


% plot
% give all vectors same length to plot
for i = 1:length(GIAerror)
    maxlength(i) = max(length(GIAerror{i}));
end
maxlength = max(maxlength);
for i = 1:length(GIAerror)
    dif = maxlength - length(GIAerror{i});
    if  dif
        GIAerror{i} = [GIAerror{i} ; zeros(dif,1)];
    end
end

dt = 1/60;
ts = dt*(1:maxlength);
numsims = 1:length(indx_allow);
if plots
    h(1) = figure;
    figure(h(1));surfc(numsims,ts,[GIAerror{indx_allow}]);xlabel('simulations of increasing \tau');title(['Simulations with GIA errror < ' num2str(thresh) ' m/s^2']);
    zlabel('GIA error (m/s^2)');ylabel('time (s)');

end