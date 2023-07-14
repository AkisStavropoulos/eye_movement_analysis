%% Joystick coefficient function (a = exp(-dt/tau)  or  tau = -dt/log(a))

% a is the input coefficient [0,1)
% tau is the timescale on which control affects the velocity v
% vmax is the maximum velocity (which will depend on tau)
% u is the control (max of 1)

% equation for joystick function: v(t) = v(t-1) + dt/tau * (-v(t-1) + vmax*u(t))  ,   u(t) ranges from 0 to 1
% vmax*u(t) = joystick_input (j_x, j_w, ...)
% equation for determining vmax based on tau: vmax(tau) = x/T * (1 + b*log(1 + exp((2*tau - T)/T*b))
% equation for tau given a: a = exp(-dt/tau)    ,   -log(a) = dt/tau

% b = 0: vmax = X/T , tau < T/2
%        vmax = (X/T)*(2*tau - T + 1) , tau > T/2
%        vmax = 
a = 0.9; % coefficient input
b = 0.5; 
dt = 1/60; % time step, 
tau = -dt/log(a); % averaging time in units of dt

vmax_control = 90; % 90 cm/s: the value that I'll choose for the PURE velocity control condition in the moog.
wmax_control = 90; % 90 deg/s: the value that I'll choose for the PURE velocity control condition in the moog.
x = 400; % mean distance in cm to targets
T = x/vmax_control; % mean duration in sec to reach targets
Nt  = round(T/dt); % number of timesteps

vmax = (x/T)*(1 + b*log(1 + exp((2*tau - T)/T*b))); % vmax as function of tau
ratio = vmax/vmax_control;
wmax = ratio*wmax_control; 

%% Option 1: Choose square wave input
j_x = vmax*([zeros(1,round(Nt/5)) ones(1,round(Nt*3/5)) zeros(1,round(Nt/5))]); % joystick translational input
j_w = wmax*([zeros(1,round(Nt/5)) ones(1,round(Nt*3/5)) zeros(1,round(Nt/5))]); % joystick rotational input

%% Option 2: Choose random input in [0,200]
j_x = vmax*rand(1,Nt); % joystick translational input
j_w = wmax*rand(1,Nt); % joystick rotational input

%% Apply the joystick coefficient function
v = zeros(1,Nt); w = zeros(1,Nt);
for t = 2:Nt
    if a == 0
        v(t) = j_x(t);
        w(t) = j_w(t);
    elseif a == 1
        v(t) = v(t-1) + j_x(t)*dt;
        w(t) = w(t-1) + j_w(t)*dt;
    else
        v(t) = v(t-1) - log(a) * (-v(t-1) + j_x(t));
        w(t) = w(t-1) - log(a) * (-w(t-1) + j_w(t));
    end
end

figure;
subplot(2,1,1);
plot([1:Nt]*dt,j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
plot([1:Nt]*dt,v);ylabel('translational input (cm/s)');xlabel('timesteps');
subplot(2,1,2);
plot([1:Nt]*dt,j_w);ylabel('rotational input (deg/s)');xlabel('timesteps');hold on;
plot([1:Nt]*dt,w);hold off;


