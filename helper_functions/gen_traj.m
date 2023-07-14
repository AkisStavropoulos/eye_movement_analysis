function [mu_xt, mu_yt, mu_x, mu_y, mu_theta, mu_thetat] = gen_traj(mu_w, mu_v, x0, y0, dt)

% generates trajectory given w,v,x,y
% linear speed mean (mu_v)
% angular speed mean (mu_w)
% initial positions (x0/y0)
% outputs: mu_xt, mu_yt, mu_x, mu_y 

% sampling rate
default_prs;
factor_downsample = prs.factor_downsample; % from default_prs
fsamp = (5000/6)/factor_downsample; % needs to match downsampling rate
if isempty(dt);   dt = (1/60);    end

% select first dimension
sz = size(mu_v);

% initialize
mu_xt = zeros(sz);
mu_yt = zeros(sz);
mu_thetat = zeros(sz);
mu_xt(1) = x0;
mu_yt(1) = y0;

% construct trajectory
for j=1:length(mu_xt)
    vt_x = mu_v(j).*sin(mu_thetat(j)); % input rads
    vt_y = mu_v(j).*cos(mu_thetat(j));
    mu_xt(j+1) = mu_xt(j) + vt_x*dt;
    mu_yt(j+1) = mu_yt(j) + vt_y*dt;
    mu_thetat(j+1) = mu_thetat(j) + (mu_w(j)*pi/180)*dt; % converts degrees to rads
    mu_thetat(j+1) = (mu_thetat(j+1)>-pi & mu_thetat(j+1)<=pi).*mu_thetat(j+1) + ...
        (mu_thetat(j+1)>pi).*(mu_thetat(j+1) - 2*pi) + (mu_thetat(j+1)<=-pi).*(mu_thetat(j+1) + 2*pi);
end

mu_x = mu_xt(end); mu_y = mu_yt(end);
% re-convert rads to degrees
mu_thetat = mu_thetat*180/pi;
mu_theta = mu_thetat(end);