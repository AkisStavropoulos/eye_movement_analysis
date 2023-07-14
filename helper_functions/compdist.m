function [err] = compdist(a,vmax,Tpulse,x)
% this function is used to help find an appropriate Tpulse input for a
% given target distance
if a > 0
    [~,~,~,~,~, dist_tr,~] = vmax2acc(a, vmax, Tpulse, [],[],[],0);
    err  = abs(dist_tr - x);
    
elseif a == 0
    [~,~,~,~,~,~, dist_tr] = vmax2acc(a, vmax, Tpulse, [],[],[],0);
    err  = abs(dist_tr - x);
end