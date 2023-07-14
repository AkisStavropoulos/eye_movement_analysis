function [mu_xt, mu_yt, mu_x, mu_y, var_x, var_y, covar_xy, theta] = gen_traj2(mu_w, mu_v, var_w, var_v, x0, y0)

% generates trajectory given 
% linear speed mean and variance (mu_v/var_v)
% angular speed mean and variance (mu_w/var_w)
% initial positions (x0/y0)
% outputs: mu_xt, mu_yt, mu_x, mu_y, var_x, var_y, covar_xy 

% sampling rate
fsamp = 5000/6; % needs to match downsampling rate, sampling rate of smr file
dt = (1/fsamp);

% select first dimension
sz = size(mu_v);

% initialize
mu_xt = zeros(sz);
var_xt = zeros(sz);
mu_yt = zeros(sz);
var_yt = zeros(sz);
mu_thetat = zeros(sz);
var_thetat = zeros(sz);
covar_xyt = zeros(sz);
mu_xt(1) = x0;
mu_yt(1) = y0;

% construct trajectory
for j=1:length(mu_xt)
    vt_x = mu_v(j).*sin(mu_thetat(j));
    vt_y = mu_v(j).*cos(mu_thetat(j));
    
    mu_xt(j+1) = mu_xt(j) + vt_x*dt;
    var_xt(j+1) = var_xt(j) + ...
        ((sin(mu_thetat(j)).^2).*var_v(j) + ...
        ((mu_v(j).*cos(mu_thetat(j))).^2).*var_thetat(j))*(dt);
    
    mu_yt(j+1) = mu_yt(j) + vt_y*dt;
    var_yt(j+1) = var_yt(j) + ...
        ((cos(mu_thetat(j))).^2.*var_v(j) + ...
        ((mu_v(j).*sin(mu_thetat(j))).^2).*var_thetat(j))*(dt);
    
    mu_thetat(j+1) = mu_thetat(j) - (mu_w(j)*pi/180)*dt;
    mu_thetat(j+1) = (mu_thetat(j+1)>-pi & mu_thetat(j+1)<=pi).*mu_thetat(j+1) + ...
        (mu_thetat(j+1)>pi).*(mu_thetat(j+1) - 2*pi) + ...
            (mu_thetat(j+1)<=-pi).*(mu_thetat(j+1) + 2*pi);
    var_thetat(j+1) = var_thetat(j) + var_w(j)*((pi/180)^2)*(dt);
    
    covar_xyt(j+1) = covar_xyt(j) + ...
        var_v(j).*sin(mu_thetat(j)).*cos(mu_thetat(j)).*(dt) ...
        + (mu_v(j)^2)*cos(2*mu_thetat(j))*(dt^2);
end

mu_x = mu_xt(end); mu_y = mu_yt(end);
var_x = var_xt(end); var_y = var_yt(end);
covar_xy = covar_xyt(end);
theta = mu_thetat(end);