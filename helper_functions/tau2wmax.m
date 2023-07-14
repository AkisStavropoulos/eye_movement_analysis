function [wmax, awmax_filt, awmax_unfilt, tau, tau_exp, angl, angvel, tt, Tbrake_in] = tau2wmax(tau, Tpulse, Tbrake, Tsim, theta, T, prints)
%% Simulation
% NO BUTTERWORTH FILTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% NEW VERSION

% a: joystick coefficient
% tau = -dt/log(a)
% Tpulse: duration of square maximum joystick input
% Tsim: length of simulation, just big enough for the plot to include the whole simulation
% theta: angle we want to travel, we use 40 deg (max target angle)
% T: time duration of trial, we use max desired time duration 5 sec
% prints: print results and plots, 1 or 0

% prints: 0 or 1 , display messages or not
minargs = 3;maxargs = 7;
narginchk(minargs,maxargs);
if isempty(Tsim)
    Tsim = Tpulse*6;
end
if isempty(Tbrake)
    Tbrake = 0;
end

if isempty(theta)
    theta = 400; % cm
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
a = exp(-dt/tau);
sw = switchtime(tau,T);
wmax = findwmax(theta,T,tau);
%% Create Pulse input
Tpulse = round(Tpulse/dt);
Tbrake = round(Tbrake/dt);
Tsim = round(Tsim/dt);
% Rt = .2; % reaction time (200ms)
% Rt = round(Rt/dt);
j_w = wmax*([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tpulse-Tbrake)]); % joystick translational input
ts = (1:Tsim)*dt;
%% Apply the joystick coefficient function
if a > 0
    
    w = zeros(1,Tsim);
    for t = 2:Tsim
        beta = (1 - a);
        w(t) = a*w(t-1) + beta*j_w(t);
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
        
    accw = diff(w)/dt;
    accw = [0 accw];
    % plot(ts,acc,'.');
    awmax_unfilt = max(accw);
    %%
    % find end of trial, when velocity < 1 cm
    [last_peak,last_peak_ind] = max(abs(w(Tpulse+Tbrake+30:end))); % find when control was let go, 30 timesteps = .5 secs, filter delay
    indx0 = find(abs(w(Tpulse+Tbrake+30:end)) <= (.63)*last_peak, 1);
    
    indx = find(abs(w(Tpulse+Tbrake+30+indx0:end))<=1,1,'first'); % calculate stop after the control was let go
    indx = Tpulse + Tbrake + indx0 + 30 + indx;
    tt = indx*dt;
    if isempty(tt)
        indx = find(abs(w)>1,1,'last');
        tt = indx*dt;
    end

    
    angvel = w;
    angl = cumsum(angvel(1:indx))*dt;
    angl = angl(end);
    Tbrake_in = Tbrake*dt;
%     if Tpulse_in ~= sw
%     brakeind = find(j_x_new < 0);
%     Tbrake_in = length(brakeind)*dt;
%     end
    %     dist = cumsum(vel*dt);
    %     dist = dist(end);
    
    
    if prints
        
        figure;
        plot(ts,j_w);title(['joystick coefficient = ' num2str(a) ', wmax = ' num2str(wmax) ', T = ' num2str(T)]);hold on;
        plot(ts,w);ylabel('rotational input (deg/s)');xlabel('time (s)');
%         if Tpulse_in ~= sw
%             plot((1:length(j_x_new))*dt,j_x_new,'m'); % new joystick input
%         end
        vline(tt); xlim([0 Tpulse_in*4]);
        legend('joystick rotational input','ang. velocity output','new joystick input');
        hold off;

    end

    
    awmax_filt = max(abs(diff(w)/dt));
    % find tau for combined acceleration control and applied filter
    indx = find(w >= (.63)*wmax, 1);
    tau_exp = indx*dt;
        %
    % check obtained tau
    if isempty(tau_exp)
        Tsim_new = Tsim*10;
        j_w_tau = wmax*(ones(1,Tsim_new)); % joystick translational input
        
        w_tau = zeros(1,Tsim_new);
        for t = 2:Tsim_new
            beta = (1 - a);
            w_tau(t) = a*w_tau(t-1) + beta*j_w_tau(t);
        end
        
        w_tau_filt = zeros(1,Tsim_new); % y(n) = 0;
        for n = 3:Tsim_new
            
            w_tau_filt(n) = FB(1)*w_tau(n) + FB(2)*w_tau(n-1) + FB(3)*w_tau(n-2);
            
            w_tau_filt(n) = w_tau_filt(n) - (FA(2)*w_tau_filt(n-1) + FA(3)*w_tau_filt(n-2));
            w_tau_filt(n) = w_tau_filt(n)/FA(1);
        end
        indx = find(w_tau_filt >= (.63)*wmax, 1);
        tau_exp = indx*dt;
    end
    %
    %%
else
    
        %% amax for velocity control edge condition
    wmax0 = theta/T;
    % Tpulse = 4;
    % Tsim = 8;
    % Tpulse = round(Tpulse/dt);
    % Tsim = round(Tsim/dt);
    j_w = wmax0*([ones(1,Tpulse) zeros(1,Tsim-Tpulse)]); % joystick translational input
    % ts = (1:Tsim)*dt;
    w0 = j_w;
    % figure;
    % plot(ts,j_x);title(['joystick coefficient = 0']);hold on;
    % plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
    % hold off;
    accw0 = diff(w0)/dt;
    accw0 = [0 accw0];
    % plot(ts,acc,'.');
    awmax0 = max(abs(accw0));
    awmax_unfilt = awmax0;
    indx = find(w0 >= (.63)*wmax0, 1);
    tau_exp = indx*dt;
    angl = cumsum(w0*dt);
    angl = angl(end);
    awmax_filt = awmax_unfilt;
    angvel = w0;
    wmax = wmax0;

    [last_peak,last_peak_ind] = max(abs(w0(Tpulse+Tbrake+30:end))); % find when control was let go, 30 timesteps = .5 secs, filter delay
    indx0 = find(abs(w0(Tpulse+Tbrake+30:end)) <= (.63)*last_peak, 1);
    
    indx = find(abs(w0(Tpulse+Tbrake+30+indx0:end))<=1,1,'first'); % calculate stop after the control was let go
    indx = Tpulse + Tbrake + indx0 + 30 + indx;
    tt = indx*dt;
    if isempty(tt)
        indx = find(abs(w0)>1,1,'last');
        tt = indx*dt;
    end
%     j_x_new = j_x;
%     brakeind = find(j_x_new < 0);
%     Tbrake_in = length(brakeind)*dt;

    if prints
        
        figure;
        plot(ts,j_w);title(['joystick coefficient = ' num2str(a)]);hold on;
        plot(ts,w0,'k');ylabel('translational input (cm/s)');xlabel('time (s)');
        vline(tt);
        hold off; legend('velocity control after filter');
        
        fprintf('for a = %4.3f, x = %4.2f, T = %4.2f, Tsim = %4.2f, Tpulse = %4.2f: \n', a,theta,T,Tsim_in,Tpulse_in)
        fprintf('    vmax = %4.2f, amax_lp_filt = %4.2f, tau = %4.2f, dist = %4.2f.\n',wmax0,awmax0,tau_exp,angl)
    end
    
        % check obtained tau
    if isempty(tau_exp)
        Tsim_new = Tsim*10;
        j_w_tau = wmax*(ones(1,Tsim_new)); % joystick translational input
        
        w_tau = zeros(1,Tsim_new);
        for t = 2:Tsim_new
            beta = (1 - a);
            w_tau(t) = a*w_tau(t-1) + beta*j_w_tau(t);
        end
        
        w_tau_filt = zeros(1,Tsim_new); % y(n) = 0;
        for n = 3:Tsim_new
            
            w_tau_filt(n) = FB(1)*w_tau(n) + FB(2)*w_tau(n-1) + FB(3)*w_tau(n-2);
            
            w_tau_filt(n) = w_tau_filt(n) - (FA(2)*w_tau_filt(n-1) + FA(3)*w_tau_filt(n-2));
            w_tau_filt(n) = w_tau_filt(n)/FA(1);
        end
        indx = find(w_tau_filt >= (.63)*wmax, 1);
        tau_exp = indx*dt;
    end
    %

end
if size(angvel,2)==1
    angvel = angvel';
end