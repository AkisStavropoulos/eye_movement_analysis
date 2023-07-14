%% Acceleration control parameters
a = .95; % choose from a = .9; a = .95; a = .99; a = .995; a = .999;
b = .25; % choose from b = .1; b = .25; b = .5; b = .75; b = .9; b = 1;
x = 600; % cm
T = 4; % sec
dt = 1/60;

tau = -dt/log(a)
vmax1 = (x/T)*(1 + b*log(1 + exp((2*tau-T)/(T*b)))) 
% 
% x = 400cm, T = 3s
% for a=.9, b=.1: tau=.16, vmax=133cm/s, amax=842cm/s^2, amax_filt = 175 cm/s^2
% for a=.95, b=.1: tau=.32, vmax=133cm/s, amax=410cm/s^2, amax_filt = 145 cm/s^2
% for a=.99, b=.1: tau=1.66, vmax=151cm/s, amax=91cm/s^2 , amax_filt = 68 cm/s^2   
% for a=.995, b=.1: tau=3.32, vmax=295cm/s, amax=89cm/s^2, amax_filt = 78 cm/s^2
% 
% for a=.9, b=.25: tau=.16, vmax=134cm/s, amax=849cm/s^2, amax_filt = 177 cm/s^2
% for a=.95, b=.25: tau=.32, vmax=135cm/s,amax=415cm/s^2, amax_filt =147 cm/s^2
% for a=.99, b=.25: tau=1.66, vmax=164cm/s, amax=99cm/s^2, amax_filt = 74 cm/s^2
% for a=.995, b=.25: tau=3.32, vmax=296cm/s, amax=89cm/s^2, amax_filt = 78 cm/s^2
%      
% for a=.9, b=.5: tau=.16, vmax=144cm/s,amax=908cm/s^2, amax_filt = 189 cm/s^2
% for a=.95, b=.5: tau=.32, vmax=146cm/s,amax=449cm/s^2 , amax_filt = 158 cm/s^2
% for a=.99, b=.5: tau=1.66, vmax=187cm/s ,amax=112cm/s^2, amax_filt = 84 cm/s^2
% for a=.995, b=.5: tau=3.32, vmax=301cm/s, amax=90cm/s^2 , amax_filt = 79 cm/s^2  
% 
% for a=.9, b=.75: tau=.16, vmax=160cm/s,amax=1010cm/s^2  , amax_filt = 210 cm/s^2
% for a=.95, b=.75: tau=.32, vmax=163cm/s,amax=503cm/s^2  , amax_filt = 178 cm/s^2
% for a=.99, b=.75: tau=1.66, vmax=210cm/s,amax=126cm/s^2 , amax_filt = 94 cm/s^2  
% for a=.995, b=.75: tau=3.32, vmax=314cm/s,amax=94cm/s^2 , amax_filt = 82 cm/s^2  
%     
% for a=.9, b=.9: tau=.16, vmax=171cm/s,amax=1081cm/s^2  , amax_filt = 225 cm/s^2 
% for a=.95, b=.9: tau=.32, vmax=175cm/s, amax=540cm/s^2 , amax_filt = 191 cm/s^2  
% for a=.99, b=.9: tau=1.66, vmax=224cm/s, amax=135cm/s^2, amax_filt = 101 cm/s^2   
% for a=.995, b=.9: tau=3.32, vmax=323cm/s, amax=97cm/s^2 , amax_filt = 85 cm/s^2   
%     
% for a=.9, b=1: tau=.16, vmax=179cm/s, amax=1131cm/s^2 , amax_filt = 235 cm/s^2  
% for a=.95, b=1: tau=.32, vmax=184cm/s, amax=564cm/s^2 , amax_filt = 200 cm/s^2  
% for a=.99, b=1: tau=1.66, vmax=232cm/s, amax=140cm/s^2 , amax_filt = 105 cm/s^2  
% for a=.995, b=1: tau=3.32, vmax=330cm/s, amax=99cm/s^2 , amax_filt = 87 cm/s^2     
%     
% for a=.9, b=1.25: tau=.16, vmax=200cm/s, amax=1262cm/s^2, amax_filt = 263 cm/s^2   
% for a=.95, b=1.25: tau=.32, vmax=204cm/s, amax=630cm/s^2, amax_filt = 223 cm/s^2   
% for a=.99, b=1.25: tau=1.66, vmax=256cm/s, amax=154cm/s^2 , amax_filt = 115 cm/s^2  
%     
% for a=.9, b=1.5: tau=.16, vmax=221cm/s, amax=1398cm/s^2, amax_filt = 291 cm/s^2   
% for a=.95, b=1.5: tau=.32, vmax=226cm/s, amax=697cm/s^2 , amax_filt = 246 cm/s^2  
% for a=.99, b=1.5: tau=1.66, vmax=279cm/s, amax=168cm/s^2, amax_filt = 125 cm/s^2      
%     
% for a=.9, b=1.75: tau=.16, vmax=242cm/s, amax=1535cm/s^2 , amax_filt = 319 cm/s^2  
% for a=.95, b=1.75: tau=.32, vmax=248cm/s, amax=765cm/s^2 , amax_filt = 270 cm/s^2  
% for a=.99, b=1.75: tau=1.66, vmax=302cm/s, amax=182cm/s^2, amax_filt = 136 cm/s^2   
%     
% for a=.9, b=2: tau=.16, vmax=265cm/s, amax=1676cm/s^2, amax_filt = 349 cm/s^2
% for a=.95, b=2: tau=.32, vmax=271cm/s, amax=834cm/s^2 , amax_filt = 295 cm/s^2  
% for a=.99, b=2: tau=1.66, vmax=325cm/s, amax=196cm/s^2 , amax_filt = 146 cm/s^2  
%     
%% Option 1: Choose square wave input
Tsim = 6; % Tsim = 2*Nt
Nt  = round(Tsim/dt); % number of timesteps
j_x = vmax1*([zeros(1,round(Nt/5)) ones(1,round(Nt*3/5)) zeros(1,round(Nt*6/5))]); % joystick translational input
%% Apply the joystick coefficient function
v = zeros(1,2*Nt);
for t = 2:2*Nt
    v(t) = v(t-1) - log(a) * (-v(t-1) + j_x(t));
end

figure;
plot([1:2*Nt]*dt,j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
plot([1:2*Nt]*dt,v);ylabel('translational input (cm/s)');xlabel('time (s)');

acc0 = diff(v)/dt;
acc0 = [0 acc0];
plot([1:2*Nt]*dt,acc0,'.');
amax = max(acc0)
%% Low pass filter
%% Apply low pass filter
FA = [1, -1.9260, 0.9286];
FB = [0.0007, 0.0013, 0.0007];
% velocity control input
v0 = zeros(1,2*Nt); % y(n) = 0;
for n = 3:2*Nt

    v0(n) = FB(1)*j_x(n) + FB(2)*j_x(n-1) + FB(3)*j_x(n-2);
    
    v0(n) = v0(n) - (FA(2)*v0(n-1) + FA(3)*v0(n-2));
    v0(n) = v0(n)/FA(1);
end

hold on;
plot([1:2*Nt]*dt,v0,'*');
hold off;
% acceleration control input
v2 = zeros(1,2*Nt); % y(n) = 0;
for n = 3:2*Nt

    v2(n) = FB(1)*v(n) + FB(2)*v(n-1) + FB(3)*v(n-2);
    
    v2(n) = v2(n) - (FA(2)*v2(n-1) + FA(3)*v2(n-2));
    v2(n) = v2(n)/FA(1);
end
hold on;
plot([1:2*Nt]*dt,v2,'*');
amax_filt = max(diff(v2)/dt)

dist = cumsum(v0*dt);
dist = dist(end)

dist = cumsum(v2*dt);
dist = dist(end)
%%
%% amax for velocity control edge condition
vmax1 = 114;
Tpulse = 4;
Tsim = 8;
Tpulse = round(Tpulse/dt);
Tsim = round(Tsim/dt);
j_x = vmax1*([ones(1,Tpulse) zeros(1,Tsim-Tpulse)]); % joystick translational input
ts = (1:Tsim)*dt;
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
figure;
plot(ts,j_x);title(['joystick coefficient = 0']);hold on;
plot(ts,v0,'*k');ylabel('translational input (cm/s)');xlabel('time (s)');
hold off;

acc0 = diff(v0)/dt;
acc0 = [0 acc0];
% plot(ts,acc,'.');
amax0 = max(acc0)




