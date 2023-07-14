function [r,theta,x,y] = eyepos2flypos(beta_l,beta_r,alpha_l,alpha_r,z)

% beta_l: left eye elevation
% beta_r: right eye elevation
% alpha_l: left eye version
% alpha_r: right eye version

thresh = 0; % ZE = -5 deg <=> Y = 800 m, ZE = -4 deg <=> Y = 1000 m, ...
beta_l(beta_l > thresh) = nan;
beta_r(beta_r > thresh) = nan;


beta = 0.5*(beta_l + beta_r);
alpha = 0.5*(alpha_l + alpha_r);
x_squared = z^2*[(tand(alpha).^2)./(tand(beta).^2)].*[(1 + tand(beta).^2)./(1 + tand(alpha).^2)];
y_squared = (z^2)./(tand(beta).^2) - x_squared;
r = sqrt(x_squared + y_squared);
% r = sqrt(z.^2./tand(beta).^2);
theta = alpha;


x = r.*sind(theta);
y = r.*cosd(theta);
