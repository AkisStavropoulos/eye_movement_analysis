function [r,theta,x,y] = eye2world(zle,zre,yle,yre,z)

% zle: left eye elevation
% zre: right eye elevation
% yle: left eye version
% yre: right eye version

thresh = 0; % ZE = -5 deg <=> Y = 800 m, ZE = -4 deg <=> Y = 1000 m, ...
zle(zle > thresh) = nan;
zre(zre > thresh) = nan;


ver_mean = 0.5*(zle + zre);
hor_mean = 0.5*(yle + yre);
x_squared = z^2*[(tand(hor_mean).^2)./(tand(ver_mean).^2)].*[(1 + tand(ver_mean).^2)./(1 + tand(hor_mean).^2)];
y_squared = (z^2)./(tand(ver_mean).^2) - x_squared;
r = sqrt(x_squared + y_squared);
% r = sqrt(z.^2./tand(beta).^2);
theta = hor_mean;


x = r.*sind(theta);
y = r.*cosd(theta);
