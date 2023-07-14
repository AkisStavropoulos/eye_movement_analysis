function [vmax, amax_filt, amax_unfilt, tau, tau_exp, dist_ac, dist_vc, vel] = coef2vmax(a, Tpulse, Tsim, x, T, prints)
%% Simulation
% full inputs: coef2vmax(a, Tpulse, Tsim, prints)
% provide at least a, b, Tpulse

% a: joystick coefficient
% tau = -dt/log(a)
% b: constant that determines how fast the desired vmax increases as a function of tau 
% Tpulse: duration of square maximum joystick input
% Tsim: length of simulation, just big enough for the plot to include the whole simulation
% x: distance we want to travel, we use 600 cm (max target distance)
% T: time duration of trial, we use max desired time duration 5 sec
% prints: print results and plots, 1 or 0

% insert parameters and see the results, to choose which fit 
% a = .99; % choose from a = .9; a = .95; a = .99; a = .995; a = .999;
% b = .25; % choose from b = .1; b = .25; b = .5; b = .75; b = .9; b = 1;
% x = 600; % cm
% T = 5; % sec
% 
% example runs: coef2vmax(.99,.5,10,[],[],[]), coef2vmax(.99,.5,10,[],[],[],1), coef2vmax(.99,.5,10,[],[],[],0)
%
% amax_filt must be less than 130 cm/s^2 (which is the value for velocity
% control acceleration after the low pass filter)
% for tau=0: vmax should be 114cm/s, which means acc_max=167cm/s^2
% actual tau (including hardcoded tau of motion cuing): a=.975: tau=1.25
%                                                       a=.99:  tau=2

% input variables sequence: a, b, Tpulse, Tsim, x, T, prints

% prints: 0 or 1 , display messages or not
minargs = 3;maxargs = 7;
narginchk(minargs,maxargs);
  if isempty(Tsim)
      Tsim = Tpulse*4;
  end
if isempty(x)
    x = 600; % cm
end
if isempty(T)
    T = 6; % sec
end
if nargin ~= 6
    if nargin == 5
        prints = 1;
    else
        error(['Wrong input. Provide one of the following: ' ...
            'a) a, b, Tpulse, Tsim,  x, T, prints ' ...
            'b) a, b, Tpulse, Tsim,  x, T ' ...
            'c) a, b, Tpulse, [], [], [] ' ...
            'c) a, b, Tpulse, [], [], [], prints '])
    end
else 
end

Tsim_in = Tsim;
Tpulse_in = Tpulse;
%% Generate Tau and vmax
dt = 1/60;
tau = -dt/log(a);
% vmax = (x/T)*(1 + b*log(1 + exp((2*tau-T)/(T*b)))) ;
vmax = (x/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
%% Create Pulse input
Tpulse = round(Tpulse/dt);
Tsim = round(Tsim/dt);
Rt = .2; % reaction time (200ms)
Rt = round(Rt/dt);
j_x = vmax*([zeros(1,Rt) ones(1,Tpulse) zeros(1,Tsim-Tpulse-Rt)]); % joystick translational input
ts = (1:Tsim)*dt;
%% Apply the joystick coefficient function
if a > 0
    
    v = zeros(1,Tsim);
    for t = 2:Tsim
        beta = (1 - a);
        v(t) = a*v(t-1) + beta*j_x(t);
    end
        
    acc = diff(v)/dt;
    acc = [0 acc];
    % plot(ts,acc,'.');
    amax_unfilt = max(acc);
    
    %% Low pass filter / Motion cuing algorithm
    %% Apply low pass filter
    FA = [1, -1.9260, 0.9286];
    FB = [0.0007, 0.0013, 0.0007];
    % velocity control input (j_x) into motion cuing algorithm / for comparison plot
    v1 = zeros(1,Tsim); % y(n) = 0;
    for n = 3:Tsim
        
        v1(n) = FB(1)*j_x(n) + FB(2)*j_x(n-1) + FB(3)*j_x(n-2);
        
        v1(n) = v1(n) - (FA(2)*v1(n-1) + FA(3)*v1(n-2));
        v1(n) = v1(n)/FA(1);
    end
    
    % acceleration control input (v) into motion cuing algorithm
    v2 = zeros(1,Tsim); % y(n) = 0;
    for n = 3:Tsim
        
        v2(n) = FB(1)*v(n) + FB(2)*v(n-1) + FB(3)*v(n-2);
        
        v2(n) = v2(n) - (FA(2)*v2(n-1) + FA(3)*v2(n-2));
        v2(n) = v2(n)/FA(1);
    end
    
    if prints
        
        figure;
        plot(ts,j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
        plot(ts,v);ylabel('translational input (cm/s)');xlabel('time (s)');
        plot(ts,v1,'*k');
        plot(ts,v2,'*g');
        legend('joystick input','velocity input (unfiltered)','vel. control case after filter','acc. control case after filter');
        hold off;

    end

    
    amax_filt = max(diff(v2)/dt);
    % find tau for combined acceleration control and applied filter
    indx = find(v2 >= (.63)*vmax, 1);
    tau_exp = indx*dt;
    %
    dist_vc = cumsum(v1*dt);
    dist_vc = dist_vc(end);
    
    dist_ac = cumsum(v2*dt);
    dist_ac = dist_ac(end);
    vel = v2;
    %%
    %% amax for velocity control edge condition
    vmax0 = 90;
    % Tpulse = 4;
    % Tsim = 8;
    % Tpulse = round(Tpulse/dt);
    % Tsim = round(Tsim/dt);
    j_x = vmax0*([ones(1,Tpulse) zeros(1,Tsim-Tpulse)]); % joystick translational input
    % ts = (1:Tsim)*dt;
    %% Apply low pass filter
    FA = [1, -1.9260, 0.9286];
    FB = [0.0007, 0.0013, 0.0007];
    % velocity control input
    v0 = zeros(1,Tsim); % y(n) = 0;
    for n = 3:Tsim
        
        v0(n) = FB(1)*j_x(n) + FB(2)*j_x(n-1) + FB(3)*j_x(n-2);
        
        v0(n) = v0(n) - (FA(2)*v0(n-1) + FA(3)*v0(n-2));
        v0(n) = v0(n)/FA(1);
    end
    % figure;
    % plot(ts,j_x);title(['joystick coefficient = 0']);hold on;
    % plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
    % hold off;
    
    acc0 = diff(v0)/dt;
    acc0 = [0 acc0];
    % plot(ts,acc,'.');
    amax0 = max(acc0);
    
    if prints
        
        fprintf('for a = %4.3f, x = %4.2f, T = %4.2f, Tsim = %4.2f, Tpulse = %4.2f: \n', a,x,T,Tsim_in,Tpulse_in)
        fprintf('    vmax = %4.2f, amax_unfilt = %4.2f, amax_filt = %4.2f, tau = %4.2f, tau_exp = %4.3f, dist_vc = %4.2f, dist_ac = %4.2f.\n',vmax,amax_unfilt,amax_filt,tau,tau_exp,dist_vc,dist_ac)
        
        if amax_filt >  amax0
            disp('acceleration(amax_filt) is larger than the moog limit acceleration for vmax=90cm/s.')
            disp(['limit amax = ' num2str(amax0) ' cm/s^2'])
        end
    end
    %

else
    
        %% amax for velocity control edge condition
    vmax0 = x/T;
    % Tpulse = 4;
    % Tsim = 8;
    % Tpulse = round(Tpulse/dt);
    % Tsim = round(Tsim/dt);
    j_x = vmax0*([ones(1,Tpulse) zeros(1,Tsim-Tpulse)]); % joystick translational input
    % ts = (1:Tsim)*dt;
    %% Apply low pass filter
    FA = [1, -1.9260, 0.9286];
    FB = [0.0007, 0.0013, 0.0007];
    % velocity control input
    v0 = zeros(1,Tsim); % y(n) = 0;
    for n = 3:Tsim
        
        v0(n) = FB(1)*j_x(n) + FB(2)*j_x(n-1) + FB(3)*j_x(n-2);
        
        v0(n) = v0(n) - (FA(2)*v0(n-1) + FA(3)*v0(n-2));
        v0(n) = v0(n)/FA(1);
    end
    % figure;
    % plot(ts,j_x);title(['joystick coefficient = 0']);hold on;
    % plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
    % hold off;
    acc0 = diff(v0)/dt;
    acc0 = [0 acc0];
    % plot(ts,acc,'.');
    amax0 = max(acc0);
    amax_unfilt = amax0;
    indx = find(v0 >= (.63)*vmax0, 1);
    tau = 0;
    tau_exp = indx*dt;
    dist_vc = cumsum(v0*dt);
    dist_vc = dist_vc(end);
    dist_ac = [];
    amax_filt = amax_unfilt;
    vel = v0;
    vmax = vmax0;
    if prints
        
        figure;
        plot(ts,j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
        plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
        hold off; legend('velocity control after filter');
        
        fprintf('for a = %4.3f, x = %4.2f, T = %4.2f, Tsim = %4.2f, Tpulse = %4.2f: \n', a,x,T,Tsim_in,Tpulse_in)
        fprintf('    vmax = %4.2f, amax_lp_filt = %4.2f, tau = %4.2f, dist_vc = %4.2f.\n',vmax0,amax0,tau_exp,dist_vc)
    end
end
if size(vel,2)==1
    vel = vel';
end