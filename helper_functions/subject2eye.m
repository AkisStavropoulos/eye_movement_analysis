function [yle,zle,yre,zre] = subject2eye(xt,yt,zt,delta)
% Transform world coordinates to eye coordinates

yle = atan2d(xt + delta, sqrt(yt.^2 + zt^2));
yre = atan2d(xt - delta, sqrt(yt.^2 + zt^2));
zle = atan2d(zt , sqrt(yt.^2 + (xt + delta).^2));
zre = atan2d(zt , sqrt(yt.^2 + (xt - delta).^2));

ver_mean = nanmean([zle , zre],2); % mean vertical eye position (of the two eyes)
hor_mean = nanmean([yle , yre],2); % mean horizontal eye position
