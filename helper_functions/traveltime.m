function travel_T = traveltime(tau,x1,vmax)
%% calculate travel time for given distance and tau
% vmax, mean x, mean T are known
if tau > 0
travel_T = 2.*tau.*acosh(exp(x1./(2.*tau.*vmax)));
elseif tau == 0
    travel_T = x1./vmax;
end