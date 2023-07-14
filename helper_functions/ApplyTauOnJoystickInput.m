%% Control dynamics manipulation
%% Overview
% j_x: X joystick position (angular component of motion)
% j_y: Y joystick position (linear component of motion)

% parameters
tau = 0.9; % choose value
x = 300; % average desired travel distance (cm)
theta = 80; % maximum target angle x2 (deg): we want to be able to turn the max angle even for closer targets
T = 3; % average desired travel time (s)
dt = 1/60; % same as refresh rate

% calculate max velocities
if tau > 0
vmax = (x/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
wmax = (theta/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
else
vmax = x/T;
wmax = theta/T;
end

% calculate velocity output
a = exp(-dt/tau);
beta = (1 - a);
v(t) = a*v(t-1) + vmax*beta*j_y(t);
w(t) = a*w(t-1) + wmax*beta*j_x(t);

%% Example (get parameters and max velocities from above section)
bang_bang = 1;
Tpulse = 2; % seconds
Tbrake = 0; % seconds
Tsim = 12; % seconds

% switchtime for bang-bang control (for a target T seconds away)
if bang_bang
    if tau > 0;    sw = tau.*log((1 + exp(T./tau))./2);
    else           sw = T;  end
    Tpulse = sw;
    Tbrake = T - sw;
end

% Create Pulse input
Tpulse = round(Tpulse/dt);
Tbrake = round(Tbrake/dt);
Tsim = round(Tsim/dt);
j_y = ([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tpulse-Tbrake)]); % joystick translational input
j_x = ([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tpulse-Tbrake)]); % joystick translational input
ts = (1:Tsim)*dt;

% Generate velocity output
v = zeros(1,Tsim);  w = zeros(1,Tsim);
for t = 2:Tsim
    beta = (1 - a);
    v(t) = a*v(t-1) + vmax*beta*j_y(t);
    w(t) = a*w(t-1) + wmax*beta*j_x(t);
end

