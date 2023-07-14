%% Joystick coefficient function (tau = dt/1-a)

a = 0.01; % joystick coefficient, a ? [0,1]
v_max = 200; % cm/s
w_max = 90; % deg/s
dt = 1/60;
%% Option 1: Choose square wave input
j_x = v_max*([zeros(1,20) ones(1,60) zeros(1,20)]);
j_w = w_max*([zeros(1,20) ones(1,60) zeros(1,20)]);

%% Option 2: Choose random input in [0,200]
j_x = v_max*rand(1,100); % joystick translational input
j_w = w_max*rand(1,100); % joystick rotational input

%% Apply the joystick coefficient function
v = zeros(1,100); w = zeros(1,100);
for t = 2:100
    if a == 0
        v(t) = j_x(t);
        w(t) = j_w(t);
    elseif a == 1
        v(t) = v(t-1) + j_x(t)*dt;
        w(t) = w(t-1) + j_w(t)*dt;
    else
        tau = dt/1-a;
        v(t) = j_x(t) + tau*v(t-1);
        w(t) = j_w(t) + tau*w(t-1);
    end
end

figure;
subplot(2,1,1);
plot(j_x);title(['joystick coefficient = ' num2str(a)]);hold on;
plot(v);ylabel('translational input (cm/s)');xlabel('time (s)');
subplot(2,1,2);
plot(j_w);ylabel('rotational input (deg/s)');xlabel('time (s)');hold on;
plot(w);hold off;
