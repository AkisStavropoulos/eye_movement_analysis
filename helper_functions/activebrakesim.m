function [tt,sw,vel,dist] = activebrakesim(a,vmax,ff_dist,T,prints)

% calculates time traveled for a given target distance given joystick coefficient and vmax
% a: JS coefficient
% vmax: cm/s
% ff_dist: cm

if nargin < 5
    prints = 0;
end

x = ff_dist;
dt = 1/60;
%% Generate Pulse input
% find appropriate pulse input for the trial
% minfun = @(Tpulse,Tbrake)(T - vmax2acc_brake(a, vmax, Tpulse, Tbrake, [],x,T,0)).^2;
minfun  = @(sw)compdist_brake(a,vmax,sw,x,T);
sw0 = T*3/4;
params0 = [sw0];
% [params,fval] = fminsearch(@(prms) minfun(prms(1),prms(2)), params0);%,[],[],[],[],[0 0],[inf inf]);
[params,fval] = fminsearch(minfun, params0);
sw = params;
Tpulse = sw;
Tbrake = T - sw;
%
if Tpulse < 5
    Tsim = 20;
else
    Tsim = Tpulse*4;
end
Tpulse = round(Tpulse/dt); % in timebins
if ~exist('Tbrake')
    brakes = input('Want brakes or not? (1 for yes, 0 for no): ');
    if brakes
        Tbrake = floor(Tpulse/4);
    else
        Tbrake = 0;
    end
else Tbrake = round(Tbrake/dt);
end
Tsim = round(Tsim/dt);
j_x = ([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tbrake-Tpulse)]); % joystick translational input
ts = (1:Tsim)*dt;
%% Apply the joystick coefficient function
if a > 0
    
    v = zeros(1,Tsim);
    for t = 2:Tsim
        beta = vmax*(1 - a);
        v(t) = a*v(t-1) + beta*j_x(t);
    end
    
    %% Low pass filter
    %% Apply low pass filter
    FA = [1, -1.9260, 0.9286];
    FB = [0.0007, 0.0013, 0.0007];
    
    % acceleration control input
    v2 = zeros(1,Tsim); % y(n) = 0;
    for n = 3:Tsim
        
        v2(n) = FB(1)*v(n) + FB(2)*v(n-1) + FB(3)*v(n-2);
        
        v2(n) = v2(n) - (FA(2)*v2(n-1) + FA(3)*v2(n-2));
        v2(n) = v2(n)/FA(1);
    end
    % find end of trial, when velocity < 1 cm
    [last_peak,last_peak_ind] = max(abs(v2(Tpulse+Tbrake+30:end))); % find when control was let go, 30 timesteps = .5 secs, filter delay
    indx0 = find(abs(v2(Tpulse+Tbrake+30:end)) <= (.63)*last_peak, 1);
    
    indx = find(abs(v2(Tpulse+Tbrake+30+indx0:end))<=1,1,'first'); % calculate stop after the control was let go
    indx = Tpulse + Tbrake + indx0 + 30 + indx;
    tt = indx*dt;
   
    if prints
        
        figure;
        plot(ts,vmax*j_x);title(['joystick coefficient = ' num2str(a) ', vmax = ' num2str(vmax) ', ff_d = ' num2str(x)]);hold on;
        plot(ts,v);ylabel('translational input (cm/s)');xlabel('time (s)');
        plot(ts,v2,'*g');
        vline(tt);
        legend('joystick input','velocity input (unfiltered)','control after MC filter');
        hold off;
        
        fprintf('for a = %4.3f, vmax = %4.2f, firefly_dist = %4.2f \n', a,vmax,x)
        fprintf('    Time traveled is t = %4.2f seconds.\n',tt)

    end
    vel = v2;
    dist = cumsum(vel(1:indx))*dt;
    dist = dist(end);
    dist;
%     amax_filt = max(diff(v2)/dt);
%     % find tau for combined acceleration control and applied filter
%     indx = find(v2 >= (.63)*vmax, 1);
%     tau_exp = indx*dt;
%     %    
%     dist_ac = cumsum(v2*dt);
%     dist_ac = dist_ac(end);
elseif a == 0
    
        %% amax for velocity control edge condition
    vmax0 = vmax;
    % Tpulse = 4;
    % Tsim = 8;
    % Tpulse = round(Tpulse/dt);
    % Tsim = round(Tsim/dt);
    j_x = vmax0*([ones(1,Tpulse) -ones(1,Tbrake) zeros(1,Tsim-Tbrake-Tpulse)]); % joystick translational input
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
        
   indx = find(v0(Tpulse+Tbrake:end)<=1,1,'first');
   indx = Tpulse + Tbrake + indx;
   tt = indx*dt;

%     acc0 = diff(v0)/dt;
%     acc0 = [0 acc0];
%     amax0 = max(acc0);
%     amax_unfilt = amax0;
%     indx = find(v0 >= (.63)*vmax0, 1);
%     tau_exp = indx*dt;
%     dist_vc = cumsum(v0*dt);
%     dist_vc = dist_vc(end);
%     dist_ac = [];
%     amax_filt = [];
%     tau = [];
    
    if prints
        
        figure;
        plot(ts,j_x);title(['joystick coefficient = ' num2str(a) ', vmax = ' num2str(vmax) ', ff_d = ' num2str(x)]);hold on;
        plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
        vline(tt);
        legend('joystick/velocity input','control after filter');hold off;
        
        fprintf('for a = %4.3f, vmax = %4.2f, firefly_dist = %4.2f \n', a,vmax,x)
        fprintf('    Time traveled is t = %4.2f seconds.\n',tt)
    end
    vel = v0;
    dist = cumsum(vel(1:indx))*dt;
    dist = dist(end);
end

%% BONUS: command to generate velocity and acceleration plot
% 
% [tt,vel] = timetravelsim(0,40,500);
% acc = diff(vel)/dt;
% acc = [0 acc];figure;plot(vel);hold on;plot(acc);legend('velocity','acceleration');