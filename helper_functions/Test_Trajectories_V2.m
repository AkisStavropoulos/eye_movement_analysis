% Sampling rate
dt = 1/60 ;

Parameters.headHeight = 0.7 ;

Parameters.t0 = 0.07 ; 
Parameters.t1 = 0.3 ; 
Parameters.t2 = 1 ;

Parameters.k0 = -0.5471; 
Parameters.k1= 1.8929; 
Parameters.k2=-0.3457;
Parameters.t3 = 3 ;

Parameters.maxBall_X_Vel = 1 ;
Parameters.maxBall_Y_Vel = 1 ;
Parameters.maxBall_Yaw_Vel = 90 ;

Parameters.Arena_Radius = 1 ; 
Parameters.Avoidance_Radius = 0.6 ; 
Parameters.MaxPos = 0.18 ;
Parameters.MaxVel = 0.3 ;
Parameters.MaxAcc = 4 ;
Parameters.MaxTiltPos = 5 ;
Parameters.MaxTiltVel = 25 ;
Parameters.MaxTiltAcc = 150 ;
Parameters.LPCutoff = 1 ;

%% Read simulation results and overlay
X = dlmread('filteredBall_vel') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,2); hold on; plot(time,X(:,1:2),':k','LineWidth',2)
subplot(4,4,2); hold on; plot(time,X(:,3)/180,':k','LineWidth',2)

X = dlmread('VR_XY_Vel') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,3); hold on; plot(time,X,':k','LineWidth',2)
X = dlmread('VR_Yaw_Vel') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,3); hold on; plot(time,X/180,':k','LineWidth',2)

X = dlmread('VR_XY_Pos') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,4); hold on; plot(time,X,':k','LineWidth',2)
X = dlmread('VR_Yaw_Pos') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,4); hold on; plot(time,(mod(X+180,360)-180)/180,':k','LineWidth',2)

X = dlmread('Rat_XY_Acc') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,6); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('desiredMoog_XY_Acc') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,7); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('MC_XY_Acc') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,8); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('MC_TiltXY_Pos') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,9); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('Moog_XY_Vel') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,11); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('Moog_XY_Pos') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,12); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('Moog_TiltXY_Acc') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,13); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('Moog_TiltXY_Vel') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,14); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('Moog_TiltXY_Pos') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,15); hold on; plot(time,X,':k','LineWidth',2)

X = dlmread('finalMoog_XY_GIA') ;time = ((1:size(X,1))/60)+t(1);
subplot(4,4,16); hold on; plot(time,X,':k','LineWidth',2)

%% Simple sinusoidal motion, should not exceed limits
f = 0.5 ; % Frequency
A = 0.05 ; % Amplitude in meters
t = (-20:dt:20)' ;
pos = A*sin(t*f*2*pi) ;
ramp = normpdf(t,0,5); ramp=ramp/max(ramp) ;
pos=pos.*ramp;

Ball_X_Vel = gradient(pos)/dt ;
Ball_Y_Vel = t*0 ; Ball_Yaw_Vel = t*0;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)

%% Velocity step, should exceed limits

t = (-5:dt:15)' ;
vel=t*0;vel(t>0&t<10)=1 ;vel(t>0&t<1)=t(t>0&t<1)  ;

Ball_X_Vel = vel*0.5 ;
Ball_Y_Vel = t*0 ; Ball_Yaw_Vel = t*0;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)
%% Circular trajectory

t = (-5:dt:35.02)' ;
vel=t*0;vel(t>0)=0.1 ;
omega = t*0;omega(t>2)=28.6479;

Ball_X_Vel = vel;
Ball_Y_Vel = t*0 ; 
Ball_Yaw_Vel = omega;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)

%% Test arena borders
t = (-5:dt:90)' ;
vel = t*0+0.5 ;
ramp = t*0;ramp(t>0&t<1)=t(t>0&t<1) ; ramp(t>=1)=1 ;
ramp(t>=60&t<62)=61-t(t>=60&t<62) ;ramp(t>=62)=-1;
ramp(t>80&t<=81)=-81+t(t>80&t<=81) ;
ramp(t>81)=0 ;
vel=vel.*ramp ;
% go forward for 10s: should get stuck on the wall
% next rotate 60 deg
omega = t*0; 
omega(t>10&t<=11)=60/2 ; % Should 'slide' and get stuck to a new position
omega(t>30&t<=40)=3 ; % Should go on sliding for 10 s


Ball_X_Vel = vel*0.5 ;
Ball_Y_Vel = t*0 ; 
Ball_Yaw_Vel = omega;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)

