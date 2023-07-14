function [mu_xt, mu_yt, mu_x, mu_y, mu_theta, mu_thetat] = gen_traj_HP(mu_w, mu_v, x0, y0, theta0)

% generates trajectory given w,v,x,y
% linear speed mean (mu_v)
% angular speed mean (mu_w)
% initial positions (x0/y0)
% outputs: mu_xt, mu_yt, mu_x, mu_y 

% sampling rate
fsamp = 5000/48; % needs to match downsampling rate
dt = (1/fsamp);

% select first dimension
sz = size(mu_v);

% initializeclo
mu_xt = zeros(sz);
mu_yt = zeros(sz);
mu_thetat = zeros(sz);
mu_xt(1) = x0;
mu_yt(1) = y0;
mu_thetat(1) = theta0;

% construct trajectory
for j=1:sz(2)
    vt_x = mu_v(j).*sin(mu_thetat(j));
    vt_y = mu_v(j).*cos(mu_thetat(j));
    mu_xt(j+1) = mu_xt(j) + vt_x*dt;
    mu_yt(j+1) = mu_yt(j) + vt_y*dt;
    mu_thetat(j+1) = mu_thetat(j) + (mu_w(j)*pi/180)*dt;
    mu_thetat(j+1) = (mu_thetat(j+1)>-pi & mu_thetat(j+1)<=pi).*mu_thetat(j+1) + ...
        (mu_thetat(j+1)>pi).*(mu_thetat(j+1) - 2*pi) + ...
            (mu_thetat(j+1)<=-pi).*(mu_thetat(j+1) + 2*pi);
end

mu_yt = mu_yt;
mu_x = mu_xt(end); mu_y = mu_yt(end);
mu_theta = mu_thetat(end);