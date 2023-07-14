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

Parameters.Arena_Radius = 1; 
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

%% Go full speed and see what happens

t = (-5:dt:15)' ;
vel=t*0;vel(t>0&t<3)=-0.1 ;  ;
vel(t>5)=10 ;
Ball_X_Vel = vel ;
Ball_Y_Vel = t*0 ; Ball_Yaw_Vel = t*0;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)

%% Velocity step, should exceed limits

t = (-5:dt:15)' ;
vel=t*0;vel(t>0&t<3)=1 ;vel(t>0&t<1)=t(t>0&t<1)  ;

Ball_X_Vel = vel*0.5 ;
Ball_Y_Vel = t*0 ; Ball_Yaw_Vel = t*0;

Motion_Cueing_Complete(t,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)

%% Rotation and velocity step

t = (-5:dt:15)' ;
omega = t*0;omega(t>0&t<10)=0 ;
vel=t*0;vel(t>12&t<13)=0.5  ;

Ball_X_Vel = vel ;
Ball_Y_Vel = t*0 ; Ball_Yaw_Vel = omega;

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

%% Load smr file

cd('C:\Users\Jean\Dropbox\Houston_Matlab_Analysis\Motion Cueing\MCtest') ;

CEDS64LoadLib('C:\Users\Jean\Dropbox\Houston_Matlab_Analysis\matson\CEDS64ML');

fhand1 = CEDS64Open('test.smr');
if fhand1 <0, error(['Could not open file ']);end
[~, BallX] = CEDS64ReadWaveF( fhand1,11, 1000000, 0 );
[~, BallY] = CEDS64ReadWaveF( fhand1,10, 1000000, 0 );
[~, Yaw] = CEDS64ReadWaveF( fhand1,7, 1000000, 0 );
CEDS64CloseAll();

Yaw=Yaw*90/2.25;
BallX=BallX(1:231890) ;
BallY=BallY(1:231890) ;
Yaw=Yaw(1:231890) ;
time = (1:size(Yaw,1))/833.33;%+Commands(99).parameters.end_time ;
[b,a] = butter(2,20/833.33,'low') ;
for i = 1:2
    BallX = filter(b,a,BallX(end:-1:1)) ;
    BallY = filter(b,a,BallY(end:-1:1)) ;
    Yaw = filter(b,a,Yaw(end:-1:1)) ;
end

t = (0:dt:time(end))' ;
BallX = interp1(time,BallX,t) ;
BallY = interp1(time,BallY,t) ;
Yaw = interp1(time,Yaw,t) ;
BallX(isnan(BallX))=0 ;
BallY(isnan(BallY))=0 ;
Yaw(isnan(Yaw))=0 ;
Motion_Cueing_Complete(t,BallX,BallY,Yaw,Parameters)
