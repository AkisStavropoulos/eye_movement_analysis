function [] = Motion_Cueing_Complete(time,Ball_X_Vel,Ball_Y_Vel,Ball_Yaw_Vel,Parameters)
%% Parameters
dt = 1/60 ;
% Parameters.headHeight = 0.7 ;

% Parameters.t0 = 0.07 ; 
% Parameters.t1 = 0.3 ; 
% Parameters.t2 = 1 ;
% Parameters.k0 = -0.5471; 
% Parameters.k1= 1.8929; 
% Parameters.k2=-0.3457;
% Parameters.t3 = 3 ;
% 
% Parameters.maxBall_X_Vel = 1 ;
% Parameters.maxBall_Y_Vel = 1 ;
% Parameters.maxBall_Yaw_Vel = 90 ;
% 
% Parameters.Arena_Radius = 1 ; 
% Parameters.Avoidance_Radius = 0.6 ; 
% Parameters.MaxPos = 0.2 ;
% Parameters.MaxVel = 0.4 ;
% Parameters.MaxAcc = 8 ;
% Parameters.MaxTiltPos = 5 ;
% Parameters.MaxTiltVel = 40 ;
% Parameters.MaxTiltAcc = 150 ;
% Parameters.LPCutoff = 1 ;

if Parameters.LPCutoff < 60
    [filter_b,filter_a]=butter(2,Parameters.LPCutoff*dt,'low');
else
    filter_b = [1 0 0]; filter_a = [1 0 0]; 
end

% Initialize
n = length(time) ;

filteredBall_X_Vel = zeros(n,1) ;
filteredBall_Y_Vel = zeros(n,1) ;
filteredBall_Yaw_Vel = zeros(n,1) ;

VR_Yaw_Pos = zeros(n,1) ;
Moog_Yaw_Pos = zeros(n,1) ;

VR_X_Vel = zeros(n,1) ;
VR_Y_Vel = zeros(n,1) ;
VR_Yaw_Vel = zeros(n,1) ;

VR_X_Pos = zeros(n,1) ;
VR_Y_Pos = zeros(n,1) ;
VR_Yaw_Pos = zeros(n,1) ;

VR_X_Acc = zeros(n,1) ;
VR_Y_Acc = zeros(n,1) ;

Rat_X_Acc = zeros(n,1) ;
Rat_Y_Acc = zeros(n,1) ;

desiredMoog_X_Acc = zeros(n,1) ;
desiredMoog_Y_Acc = zeros(n,1) ;

MC_InternalVariable0_X = 0 ;
MC_InternalVariable1_X = 0 ;
MC_InternalVariable2_X = 0 ;
MC_InternalVariable0_Y = 0 ;
MC_InternalVariable1_Y = 0 ;
MC_InternalVariable2_Y = 0 ;

MC_X_Acc = zeros(n,1) ;
MC_Y_Acc = zeros(n,1) ;

MC_TiltX_Pos = zeros(n,1) ;
MC_TiltY_Pos = zeros(n,1) ;
MC_TiltX_Vel = zeros(n,1) ;
MC_TiltY_Vel = zeros(n,1) ;
MC_TiltX_Acc = zeros(n,1) ;
MC_TiltY_Acc = zeros(n,1) ;

Moog_X_Pos = zeros(n,1) ;
Moog_Y_Pos = zeros(n,1) ;
Moog_X_Vel = zeros(n,1) ;
Moog_Y_Vel = zeros(n,1) ;
Moog_X_Acc = zeros(n,1) ;
Moog_Y_Acc = zeros(n,1) ;

Moog_TiltX_Pos = zeros(n,1) ;
Moog_TiltY_Pos = zeros(n,1) ;
Moog_TiltX_Vel = zeros(n,1) ;
Moog_TiltY_Vel = zeros(n,1) ;
Moog_TiltX_Acc = zeros(n,1) ;
Moog_TiltY_Acc = zeros(n,1) ;

finalMoog_X_Vel = zeros(n,1) ;
finalMoog_Y_Vel = zeros(n,1) ;
finalMoog_X_Acc = zeros(n,1) ;
finalMoog_Y_Acc = zeros(n,1) ;

finalMoog_X_GIA = zeros(n,1) ;
finalMoog_Y_GIA = zeros(n,1) ;

Rat_X_GIA = zeros(n,1) ;
Rat_Y_GIA = zeros(n,1) ;

Rat_X_GIAerror = zeros(n,1) ;
Rat_Y_GIAerror = zeros(n,1) ;

VR_X_GIAerror = zeros(n,1) ;
VR_Y_GIAerror = zeros(n,1) ;

finalVR_X_Vel = zeros(n,1) ;
finalVR_Y_Vel = zeros(n,1) ;

for i = 3:length(time)

    % ***** STEP 1 ***** 
    if Ball_X_Vel(i) > Parameters.maxBall_X_Vel, Ball_X_Vel(i) = Parameters.maxBall_X_Vel;end
    if Ball_X_Vel(i) < -Parameters.maxBall_X_Vel, Ball_X_Vel(i) = -Parameters.maxBall_X_Vel;end
    if Ball_Y_Vel(i) > Parameters.maxBall_Y_Vel, Ball_Y_Vel(i) = Parameters.maxBall_Y_Vel;end
    if Ball_Y_Vel(i) < -Parameters.maxBall_Y_Vel, Ball_Y_Vel(i) = -Parameters.maxBall_Y_Vel;end
    if Ball_Yaw_Vel(i) > Parameters.maxBall_Yaw_Vel, Ball_Yaw_Vel(i) = Parameters.maxBall_Yaw_Vel;end
    if Ball_Yaw_Vel(i) < -Parameters.maxBall_Yaw_Vel, Ball_Yaw_Vel(i) = -Parameters.maxBall_Yaw_Vel;end
    
    filteredBall_X_Vel(i) = (filter_b(1)*Ball_X_Vel(i)+filter_b(2)*Ball_X_Vel(i-1)+filter_b(3)*Ball_X_Vel(i-2)-filter_a(2)*filteredBall_X_Vel(i-1)-filter_a(3)*filteredBall_X_Vel(i-2))/filter_a(1) ;
    filteredBall_Y_Vel(i) = (filter_b(1)*Ball_Y_Vel(i)+filter_b(2)*Ball_Y_Vel(i-1)+filter_b(3)*Ball_Y_Vel(i-2)-filter_a(2)*filteredBall_Y_Vel(i-1)-filter_a(3)*filteredBall_Y_Vel(i-2))/filter_a(1) ;
    filteredBall_Yaw_Vel(i) = (filter_b(1)*Ball_Yaw_Vel(i)+filter_b(2)*Ball_Yaw_Vel(i-1)+filter_b(3)*Ball_Yaw_Vel(i-2)-filter_a(2)*filteredBall_Yaw_Vel(i-1)-filter_a(3)*filteredBall_Yaw_Vel(i-2))/filter_a(1) ;

    % -------- PART 1: Virtual Reality Computations
    
    % ***** STEP 2 ***** 
    % Convert velocity in VR coordinates            
    VR_X_Vel(i) = filteredBall_X_Vel(i)*cosd(VR_Yaw_Pos(i-1))-filteredBall_Y_Vel(i)*sind(VR_Yaw_Pos(i-1)) ;
    VR_Y_Vel(i) = filteredBall_X_Vel(i)*sind(VR_Yaw_Pos(i-1))+filteredBall_Y_Vel(i)*cosd(VR_Yaw_Pos(i-1)) ;
    VR_Yaw_Vel(i) = filteredBall_Yaw_Vel(i) ;
          
    % ***** STEP 3 ***** 
    % Slow down when you get close to the walls
    VR_radius = sqrt(VR_X_Pos(i-1)^2+VR_Y_Pos(i-1)^2) ;
    if VR_radius>Parameters.Avoidance_Radius
        
        VR_centrifuge_velocity = (VR_X_Vel(i)*VR_X_Pos(i-1)+VR_Y_Vel(i)*VR_Y_Pos(i-1))/VR_radius^2 ;
        VR_tangential_velocity = (VR_Y_Vel(i)*VR_X_Pos(i-1)-VR_X_Vel(i)*VR_Y_Pos(i-1))/VR_radius^2 ;
        
        if VR_centrifuge_velocity>0           
            r = (VR_radius-Parameters.Avoidance_Radius)/(Parameters.Arena_Radius-Parameters.Avoidance_Radius) ;
            g = 1-1.2*r.^2; g(r<0)=1;            
            VR_centrifuge_velocity = VR_centrifuge_velocity*g ;
        end
        
        VR_X_Vel(i) = VR_centrifuge_velocity*VR_X_Pos(i-1)-VR_tangential_velocity*VR_Y_Pos(i-1) ;
        VR_Y_Vel(i) = VR_centrifuge_velocity*VR_Y_Pos(i-1)+VR_tangential_velocity*VR_X_Pos(i-1) ;
    end
       
    % ***** STEP 4 *****
    % Acceleration in VR coordinates
    VR_X_Acc(i) = (VR_X_Vel(i)-VR_X_Vel(i-1))/dt ;
    VR_Y_Acc(i) = (VR_Y_Vel(i)-VR_Y_Vel(i-1))/dt  ;    
    
    % Project to Rat coordinates
    Rat_X_Acc(i) = VR_X_Acc(i)*cosd(VR_Yaw_Pos(i-1))+VR_Y_Acc(i)*sind(VR_Yaw_Pos(i-1)) ;
    Rat_Y_Acc(i) = VR_X_Acc(i)*sind(VR_Yaw_Pos(i-1))-VR_Y_Acc(i)*cosd(VR_Yaw_Pos(i-1)) ;
      
    
    % -------- PART 2: MOOG computations
    % ***** STEP 5 ***** 
    % Here Moog_Yaw_Pos encodes the position of the yaw axis (need to read it)
    % Transform acceleration to Moog coordinates based on Moog position (rotation)
    desiredMoog_X_Acc(i) = Rat_X_Acc(i)*cosd(Moog_Yaw_Pos(i-1))-Rat_Y_Acc(i)*sind(Moog_Yaw_Pos(i-1)) ;
    desiredMoog_Y_Acc(i) = Rat_X_Acc(i)*sind(Moog_Yaw_Pos(i-1))+Rat_Y_Acc(i)*cosd(Moog_Yaw_Pos(i-1)) ;
    
%      Rat_X_GIA(i) = finalMoog_X_GIA(i)*cosd(Moog_Yaw_Pos(i-1))+finalMoog_Y_GIA(i)*sind(Moog_Yaw_Pos(i-1)) ;
%     Rat_Y_GIA(i) = finalMoog_X_GIA(i)*sind(Moog_Yaw_Pos(i-1))-finalMoog_Y_GIA(i)*cosd(Moog_Yaw_Pos(i-1)) ;
   
    
    
    % ***** STEP 6 ***** 
    % Apply motion cueing algorithm
    MC_InternalVariable0_X = exp(-dt/Parameters.t0)*MC_InternalVariable0_X + desiredMoog_X_Acc(i)-desiredMoog_X_Acc(i-1) ;
    MC_InternalVariable1_X = exp(-dt/Parameters.t1)*MC_InternalVariable1_X + desiredMoog_X_Acc(i)-desiredMoog_X_Acc(i-1) ;
    MC_InternalVariable2_X = exp(-dt/Parameters.t2)*MC_InternalVariable2_X + desiredMoog_X_Acc(i)-desiredMoog_X_Acc(i-1) ;
    MC_X_Acc(i) = Parameters.k0*MC_InternalVariable0_X + Parameters.k1*MC_InternalVariable1_X + Parameters.k2*MC_InternalVariable2_X ;
       
    MC_InternalVariable0_Y = exp(-dt/Parameters.t0)*MC_InternalVariable0_Y + desiredMoog_Y_Acc(i)-desiredMoog_Y_Acc(i-1) ;
    MC_InternalVariable1_Y = exp(-dt/Parameters.t1)*MC_InternalVariable1_Y + desiredMoog_Y_Acc(i)-desiredMoog_Y_Acc(i-1) ;
    MC_InternalVariable2_Y = exp(-dt/Parameters.t2)*MC_InternalVariable2_Y + desiredMoog_Y_Acc(i)-desiredMoog_Y_Acc(i-1) ;
    MC_Y_Acc(i) = Parameters.k0*MC_InternalVariable0_Y + Parameters.k1*MC_InternalVariable1_Y + Parameters.k2*MC_InternalVariable2_Y ;
        
    MC_TiltX_Pos(i) = asind((desiredMoog_X_Acc(i) - MC_X_Acc(i))/9.81) ;
    MC_TiltY_Pos(i) = asind((desiredMoog_Y_Acc(i) - MC_Y_Acc(i))/9.81) ;
    
    % ***** STEP 7 ***** 
    % Tilt control: go as fast as possible towards MC
    MC_TiltX_Vel(i) = (MC_TiltX_Pos(i) - Moog_TiltX_Pos(i-1))/dt ; 
    MC_TiltY_Vel(i) = (MC_TiltY_Pos(i) - Moog_TiltY_Pos(i-1))/dt ; 
    
    MC_TiltX_Acc(i) = (MC_TiltX_Vel(i) - Moog_TiltX_Vel(i-1))/dt ; 
    MC_TiltY_Acc(i) = (MC_TiltY_Vel(i) - Moog_TiltY_Vel(i-1))/dt ; 
    
    % Compute Moog motor command, taking limits into account
    beta_function = @(x,xmax)(xmax*(-1+2*power(1+exp(-2*x/xmax),-1))) ;
    
    Moog_X_Acc(i) = beta_function(MC_X_Acc(i)+Parameters.headHeight*MC_TiltX_Acc(i)*pi/180,Parameters.MaxAcc) ;
    Moog_Y_Acc(i) = beta_function(MC_Y_Acc(i)+Parameters.headHeight*MC_TiltY_Acc(i)*pi/180,Parameters.MaxAcc) ;    
    Moog_X_Vel(i) = beta_function(Moog_X_Vel(i-1)+dt*Moog_X_Acc(i)-Moog_X_Pos(i-1)*dt/Parameters.t3,Parameters.MaxVel) ;
    Moog_Y_Vel(i) = beta_function(Moog_Y_Vel(i-1)+dt*Moog_Y_Acc(i)-Moog_Y_Pos(i-1)*dt/Parameters.t3,Parameters.MaxVel) ;
    Moog_X_Pos(i) = beta_function((1-dt/Parameters.t3)*Moog_X_Pos(i-1)+dt*Moog_X_Vel(i),Parameters.MaxPos) ;
    Moog_Y_Pos(i) = beta_function((1-dt/Parameters.t3)*Moog_Y_Pos(i-1)+dt*Moog_Y_Vel(i),Parameters.MaxPos) ;
    
    Moog_TiltX_Acc(i) = beta_function(MC_TiltX_Acc(i),Parameters.MaxTiltAcc) ;
    Moog_TiltY_Acc(i) = beta_function(MC_TiltY_Acc(i),Parameters.MaxTiltAcc) ;    
    Moog_TiltX_Vel(i) = beta_function(Moog_TiltX_Vel(i-1)+dt*Moog_TiltX_Acc(i),Parameters.MaxTiltVel) ;
    Moog_TiltY_Vel(i) = beta_function(Moog_TiltY_Vel(i-1)+dt*Moog_TiltY_Acc(i),Parameters.MaxTiltVel) ;
    Moog_TiltX_Pos(i) = beta_function(Moog_TiltX_Pos(i-1)+dt*Moog_TiltX_Vel(i),Parameters.MaxTiltPos) ;
    Moog_TiltY_Pos(i) = beta_function(Moog_TiltY_Pos(i-1)+dt*Moog_TiltY_Vel(i),Parameters.MaxTiltPos) ;
    
    % ***** STEP 8 ***** 
    % --- The motor commands to the Moog will be:
    % Moog_X_Pos, Moog_Y_Pos, Moog_TiltX_Pos, Moog_TiltY_Pos
    
    % ***** STEP 9 *****
    % --- The linear acceleration from the Moog will be:
    finalMoog_X_Vel(i) = (Moog_X_Pos(i) - Moog_X_Pos(i-1))/dt ; 
    finalMoog_Y_Vel(i) = (Moog_Y_Pos(i) - Moog_Y_Pos(i-1))/dt ;     
    finalMoog_X_Acc(i) = (finalMoog_X_Vel(i) - finalMoog_X_Vel(i-1))/dt ; 
    finalMoog_Y_Acc(i) = (finalMoog_Y_Vel(i) - finalMoog_Y_Vel(i-1))/dt ; 
    
    finalMoog_X_GIA(i) = finalMoog_X_Acc(i) + sind(Moog_TiltX_Pos(i))*9.81 - Parameters.headHeight*Moog_TiltX_Acc(i)*pi/180 ;
    finalMoog_Y_GIA(i) = finalMoog_Y_Acc(i) + sind(Moog_TiltY_Pos(i))*9.81 - Parameters.headHeight*Moog_TiltY_Acc(i)*pi/180 ; 

    % -------- PART 3: GIA error and feedback to the VR system
    % ***** STEP 10 *****
    % Project to Rat coordinates
    Rat_X_GIA(i) = finalMoog_X_GIA(i)*cosd(Moog_Yaw_Pos(i-1))+finalMoog_Y_GIA(i)*sind(Moog_Yaw_Pos(i-1)) ;
    Rat_Y_GIA(i) = finalMoog_X_GIA(i)*sind(Moog_Yaw_Pos(i-1))-finalMoog_Y_GIA(i)*cosd(Moog_Yaw_Pos(i-1)) ;
    
    Rat_X_GIAerror(i) = Rat_X_GIA(i)-Rat_X_Acc(i) ;
    Rat_Y_GIAerror(i) = Rat_Y_GIA(i)-Rat_Y_Acc(i) ;
    
    % ***** STEP 11 *****
    % Transform to VR coordinates
    VR_X_GIAerror(i) = Rat_X_GIAerror(i)*cosd(VR_Yaw_Pos(i-1))-Rat_Y_GIAerror(i)*sind(VR_Yaw_Pos(i-1)) ;
    VR_Y_GIAerror(i) = Rat_X_GIAerror(i)*sind(VR_Yaw_Pos(i-1))+Rat_Y_GIAerror(i)*cosd(VR_Yaw_Pos(i-1)) ;

    finalVR_X_Vel(i) = VR_X_Vel(i) - VR_X_GIAerror(i)*dt ;
    finalVR_Y_Vel(i) = VR_Y_Vel(i) - VR_Y_GIAerror(i)*dt ;
    
    % ***** STEP 12 *****    
    % Update position in VR coordinates
    VR_X_Pos(i) = VR_X_Pos(i-1) + dt*finalVR_X_Vel(i) ;
    VR_Y_Pos(i) = VR_Y_Pos(i-1) + dt*finalVR_Y_Vel(i) ;
    VR_Yaw_Pos(i) = VR_Yaw_Pos(i-1) + dt*VR_Yaw_Vel(i) ;
    
    % Moog Yaw position signal
    Moog_Yaw_Pos(i) = VR_Yaw_Pos(i);
    
end


%% Graphical output
clf
clear hsub
for i = 1:16
    hsub(i)=subplot(4,4,i) ;hold on
    grid on
end
linkaxes(hsub,'x') ;

% Motion command
subplot(hsub(1)) ;
plot(time,[Ball_X_Vel Ball_Y_Vel]) ;
plot(time,Ball_Yaw_Vel/180) ;
axis([time(1) time(end) -1 1]) ; 
title('Ball Velocity Input')

% Motion command - filtered
subplot(hsub(2)) ;
plot(time,[filteredBall_X_Vel filteredBall_Y_Vel]) ;
plot(time,filteredBall_Yaw_Vel/180) ;
axis([time(1) time(end) -1 1]) ; 
title('Ball Filtered Velocity Input')

% Velocity in VR coordinates
subplot(hsub(3)) ;
% plot(time,[VR_X_Vel VR_Y_Vel]) ;
plot(time,[finalVR_X_Vel finalVR_Y_Vel]) ;
plot(time,VR_Yaw_Vel/180) ;
axis([time(1) time(end) -1 1]) ; 
title('FINAL VR Velocity')

% Position in VR coordinates
subplot(hsub(4)) ;
plot(time,[VR_X_Pos VR_Y_Pos]) ;
plot(time,mod(VR_Yaw_Pos/180+1,2)-1) ;
axis([time(1) time(end) -1 1]) ; 
title('VR Position')

% Acceleration in VR coordinates
subplot(hsub(5)) ;
plot(time,[VR_X_Acc VR_Y_Acc]) ;
axis([time(1) time(end) -10 10]) ; 
title('VR Acceleration')

% Acceleration in Rat coordinates
subplot(hsub(6)) ;
plot(time,[Rat_X_Acc Rat_Y_Acc]) ;
axis([time(1) time(end) -10 10]) ; 
title('Rat Acceleration')

% Acceleration in Moog coordinates
subplot(hsub(7)) ;
plot(time,[desiredMoog_X_Acc desiredMoog_Y_Acc]) ;
axis([time(1) time(end) -10 10]) ; 
title('desired Moog Acceleration')

% Motion Cueing Acceleration Output
subplot(hsub(8)) ;
plot(time,[MC_X_Acc MC_Y_Acc]) ;
axis([time(1) time(end) -10 10]) ; 
title('Motion Cueing Acceleration Output')

% Motion Cueing Tilt Output
subplot(hsub(9)) ;
plot(time,[MC_TiltX_Pos MC_TiltY_Pos]) ;
axis([time(1) time(end) -5 5]) ; 
title('Motion Cueing Tilt Position Output')

% Motion Cueing Acceleration Output
subplot(hsub(10)) ;
% plot(time,[MC_X_Acc],'--','Color',[0 0 1]*0.5,'LineWidth',2) ;
% plot(time,[MC_Y_Acc],'--','Color',[0 0.5 0]*0.3,'LineWidth',2) ;
plot(time,[finalMoog_X_Acc],'-','Color',[0 0 1],'LineWidth',1) ;
plot(time,[finalMoog_Y_Acc],'-','Color',[0 0.5 0],'LineWidth',1) ;
axis([time(1) time(end) -10 10]) ; 
title('Moog final linear acceleration')

% Motion Cueing vel Output
subplot(hsub(11)) ;
plot(time,[Moog_X_Vel Moog_Y_Vel]) ;
axis([time(1) time(end) -0.5 0.5]) ; 
title('Moog linear velocity command')

% Motion Cueing pos Output
subplot(hsub(12)) ;
plot(time,[Moog_X_Pos Moog_Y_Pos]) ;
axis([time(1) time(end) -0.3 0.3]) ; 
title('Moog linear position command')

% Motion Cueing Tilt Output
subplot(hsub(13)) ;
plot(time,[Moog_TiltX_Acc Moog_TiltY_Acc]) ;
axis([time(1) time(end) -100 100]) ; 
title('Moog angular acc command')

% Motion Cueing Tilt Output
subplot(hsub(14)) ;
plot(time,[Moog_TiltX_Vel Moog_TiltY_Vel]) ;
axis([time(1) time(end) -30 30]) ; 
title('Moog angular vel command')

% Motion Cueing Tilt Output
subplot(hsub(15)) ;
% plot(time,[MC_TiltX_Pos],'--','Color',[0 0 1]*0.5,'LineWidth',2) ;
% plot(time,[MC_TiltY_Pos],'--','Color',[0 0.5 0]*0.3,'LineWidth',2) ;
plot(time,[Moog_TiltX_Pos],'-','Color',[0 0 1],'LineWidth',1) ;
plot(time,[Moog_TiltY_Pos],'-','Color',[0 0.5 0],'LineWidth',1) ;
axis([time(1) time(end) -5 5]) ; 
title('Moog angular position command')

% Final Moog GIA
subplot(hsub(16)) ;
% plot(time,[desiredMoog_X_Acc ],'--','Color',[0 0 1]*0.5,'LineWidth',2) ;
% plot(time,[desiredMoog_Y_Acc ],'--','Color',[0 0.5 0]*0.3,'LineWidth',2) ;
plot(time,[finalMoog_X_GIA],'-','Color',[0 0 1],'LineWidth',1) ;
plot(time,[finalMoog_Y_GIA],'-','Color',[0 0.5 0],'LineWidth',1) ;


axis([time(1) time(end) -10 10]) ; 
title('Final Moog GIA')

linkaxes(hsub([5 6 7 8 10 16]),'y')
