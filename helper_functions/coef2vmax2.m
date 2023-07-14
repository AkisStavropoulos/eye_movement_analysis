function [vmax, amax_filt, amax_unfilt, tau, tau_exp, dist, vel, tt, Tbrake_in] = coef2vmax2(a, Tpulse, Tbrake, Tsim, x, T, prints)
%% Simulation
% NEW VERSION
% full inputs: coef2vmax(a, Tpulse, Tbrake,Tsim, x, T, prints)
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
    Tsim = Tpulse*6;
end
if isempty(Tbrake)
    Tbrake = 0;
end

if isempty(x)
    x = 400; % cm
end
if isempty(T)
    T = 7; % sec
end
if nargin ~= 7
    if nargin == 6
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
% 
% if Tbrake < 0
%     Tbrake = 4;
% end
Tsim_in = Tsim;
Tpulse_in = Tpulse;
Tbrake_in = Tbrake;

%% Generate Tau and vmax
dt = 1/60;
tau = -dt/log(a);
sw = switchtime(tau,T);
% vmax = (x/T)*(1 + b*log(1 + exp((2*tau-T)/(T*b)))) ;
vmax = (x/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
%% Create Pulse input
Tpulse = round(Tpulse/dt);
Tbrake = round(Tbrake/dt);
Tsim = round(Tsim/dt);
% Rt = .2; % reaction time (200ms)
% Rt = round(Rt/dt);
j_x = vmax*([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tpulse-Tbrake)]); % joystick translational input
ts = (1:Tsim)*dt;
%% Apply the joystick coefficient function
if a > 0
    
    v = zeros(1,Tsim);
    for t = 2:Tsim
        beta = (1 - a);
        v(t) = a*v(t-1) + beta*j_x(t);
%         
%         if Tpulse_in ~= sw
%         if t > Tpulse
%             if v(t) <= 0
%                 disp('went in the loop to terminate braking');
%                 v(t:Tsim) = 0;
%                 Tbrake = t - Tpulse + 1;
%                 j_x_new = vmax*([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tpulse-Tbrake)]);
%                 break;
%             end
%         end
%         end
        
    end
        
    acc = diff(v)/dt;
    acc = [0 acc];
    % plot(ts,acc,'.');
    amax_unfilt = max(acc);
    
    %% Low pass filter / Motion cueing algorithm
    %% Apply low pass filter
%     FA = [1, -1.9260, 0.9286];
%     FB = [0.0007, 0.0013, 0.0007];
%     % velocity control input (j_x) into motion cuing algorithm / for comparison plot
% %     v1 = zeros(1,Tsim); % y(n) = 0;
% %     for n = 3:Tsim
% %         
% %         v1(n) = FB(1)*j_x(n) + FB(2)*j_x(n-1) + FB(3)*j_x(n-2);
% %         
% %         v1(n) = v1(n) - (FA(2)*v1(n-1) + FA(3)*v1(n-2));
% %         v1(n) = v1(n)/FA(1);
% %     end
%     
%     % acceleration control input (v) into motion cuing algorithm
%     v2 = zeros(1,Tsim); % y(n) = 0;
%     for n = 3:Tsim
%         
%         v2(n) = FB(1)*v(n) + FB(2)*v(n-1) + FB(3)*v(n-2);
%         
%         v2(n) = v2(n) - (FA(2)*v2(n-1) + FA(3)*v2(n-2));
%         v2(n) = v2(n)/FA(1);
%     end
    % find end of trial, when velocity < 1 cm
    [last_peak,last_peak_ind] = max(abs(v(Tpulse+Tbrake+30:end))); % find when control was let go, 30 timesteps = .5 secs, filter delay
    indx0 = find(abs(v(Tpulse+Tbrake+30:end)) <= (.63)*last_peak, 1);
    
    indx = find(abs(v(Tpulse+Tbrake+30+indx0:end))<=1,1,'first'); % calculate stop after the control was let go
    indx = Tpulse + Tbrake + indx0 + 30 + indx;
    tt = indx*dt;
    if isempty(tt)
        indx = find(abs(v)>1,1,'last');
        tt = indx*dt;
    end

    
    vel = v;
    dist = cumsum(vel(1:indx))*dt;
    dist = dist(end);
    Tbrake_in = Tbrake*dt;
%     if Tpulse_in ~= sw
%     brakeind = find(j_x_new < 0);
%     Tbrake_in = length(brakeind)*dt;
%     end
    %     dist = cumsum(vel*dt);
    %     dist = dist(end);
    
    
    if prints
        
        figure;
        plot(ts,j_x);title(['\tau = ' num2str(tau) ', vmax = ' num2str(vmax) ', T = ' num2str(T) ', ff_dist = ' num2str(dist)]);hold on;
        plot(ts,v);ylabel('translational input (cm/s)');xlabel('time (s)');
%         plot(ts,v1,'*k');
%         if Tpulse_in ~= sw
%             plot((1:length(j_x_new))*dt,j_x_new,'m'); % new joystick input
%         end
        vline(tt); xlim([0 Tpulse_in*4]);
        legend('joystick input','velocity output');
        hold off;

    end

    
    amax_filt = max(abs(diff(v)/dt));
    % find tau for combined acceleration control and applied filter
    indx = find(v >= (.63)*vmax, 1);
    tau_exp = indx*dt;
    %
    % check obtained tau
    if isempty(tau_exp)
        Tsim_new = Tsim*10;
        j_x_tau = vmax*(ones(1,Tsim_new)); % joystick translational input
        
        v_tau = zeros(1,Tsim_new);
        for t = 2:Tsim_new
            beta = (1 - a);
            v_tau(t) = a*v_tau(t-1) + beta*j_x_tau(t);
        end
        
        v_tau_filt = zeros(1,Tsim_new); % y(n) = 0;
        for n = 3:Tsim_new
            
            v_tau_filt(n) = FB(1)*v_tau(n) + FB(2)*v_tau(n-1) + FB(3)*v_tau(n-2);
            
            v_tau_filt(n) = v_tau_filt(n) - (FA(2)*v_tau_filt(n-1) + FA(3)*v_tau_filt(n-2));
            v_tau_filt(n) = v_tau_filt(n)/FA(1);
        end
        indx = find(v_tau_filt >= (.63)*vmax, 1);
        tau_exp = indx*dt;
    end
    
    %%
    %% amax for velocity control edge condition
    vmax0 = 90;
    % Tpulse = 4;
    % Tsim = 8;
    % Tpulse = round(Tpulse/dt);
    % Tsim = round(Tsim/dt);
    j_x0 = vmax0*([ones(1,Tpulse) zeros(1,Tsim-Tpulse)]); % joystick translational input
    % ts = (1:Tsim)*dt;
    %
    %
%     if prints
%         
%         fprintf('for a = %4.3f, x = %4.2f, T = %4.2f, Tsim = %4.2f, Tpulse = %4.2f: \n', a,x,T,Tsim_in,Tpulse_in)
%         fprintf('    vmax = %4.2f, amax_unfilt = %4.2f, amax_filt = %4.2f, tau = %4.2f, tau_exp = %4.3f, dist = %4.2f.\n',vmax,amax_unfilt,amax_filt,tau,tau_exp,dist)
%         
%         if amax_filt >  amax0
%             disp('acceleration(amax_filt) is larger than the moog limit acceleration for vmax=90cm/s.')
%             disp(['limit amax = ' num2str(amax0) ' cm/s^2'])
%         end
%     end
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
    amax0 = max(abs(acc0));
    amax_unfilt = amax0;
    indx = find(v0 >= (.63)*vmax0, 1);
    tau = 0;
    tau_exp = indx*dt;
    dist = cumsum(v0*dt);
    dist = dist(end);
    amax_filt = amax_unfilt;
    vel = v0;
    vmax = vmax0;
    
    [last_peak,last_peak_ind] = max(abs(v0(Tpulse+Tbrake+30:end))); % find when control was let go, 30 timesteps = .5 secs, filter delay
    indx0 = find(abs(v0(Tpulse+Tbrake+30:end)) <= (.63)*last_peak, 1);
    
    indx = find(abs(v0(Tpulse+Tbrake+30+indx0:end))<=1,1,'first'); % calculate stop after the control was let go
    indx = Tpulse + Tbrake + indx0 + 30 + indx;
    tt = indx*dt;
    if isempty(tt)
        indx = find(abs(v0)>1,1,'last');
        tt = indx*dt;
    end
    %     j_x_new = j_x;
    %     brakeind = find(j_x_new < 0);
    %     Tbrake_in = length(brakeind)*dt;
        %
    % check obtained tau
    if isempty(tau_exp)
        Tsim_new = Tsim*10;
        j_x_tau = vmax*(ones(1,Tsim_new)); % joystick translational input
        
        v_tau_filt = zeros(1,Tsim_new); % y(n) = 0;
        for n = 3:Tsim_new
            
            v_tau_filt(n) = FB(1)*j_x_tau(n) + FB(2)*j_x_tau(n-1) + FB(3)*j_x_tau(n-2);
            
            v_tau_filt(n) = v_tau_filt(n) - (FA(2)*v_tau_filt(n-1) + FA(3)*v_tau_filt(n-2));
            v_tau_filt(n) = v_tau_filt(n)/FA(1);
        end
        indx = find(v_tau_filt >= (.63)*vmax, 1);
        tau_exp = indx*dt;
    end

    if prints
        
        figure;
        plot(ts,j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
        plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
        vline(tt);
        hold off; legend('velocity control after filter');
        
        fprintf('for a = %4.3f, x = %4.2f, T = %4.2f, Tsim = %4.2f, Tpulse = %4.2f: \n', a,x,T,Tsim_in,Tpulse_in)
        fprintf('    vmax = %4.2f, amax_lp_filt = %4.2f, tau = %4.2f, dist = %4.2f.\n',vmax0,amax0,tau_exp,dist)
    end
end
if isempty(tau_exp)
    keyboard;
end
if size(vel,2)==1
    vel = vel';
end