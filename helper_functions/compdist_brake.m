function [err] = compdist_brake(a,vmax,sw,x,T)
% this function is used to help find an appropriate Tpulse input for a
% given target distance
% if a > 0
    [tt,~,~,~,~,~, dist_tr,~] = vmax2acc_brake(a, vmax, sw, [],x,T,0);
%    err  = (dist_tr - x).^2;
    err_x = (dist_tr - x).^2;
    err_t = (T - tt).^2;
    err = err_x + err_t;
% elseif a == 0
%     [~,~,~,~,~, dist_tr,~,tt] = vmax2acc_brake(a, vmax, Tpulse, Tbrake, [],x,T,0);
% %    err  = (dist_tr - x).^2;
%     err = (T - tt).^2;
% end