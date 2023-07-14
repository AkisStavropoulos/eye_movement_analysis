%% Test Moog and sensors
dt = 1/60; % 60 Hz
HeadHeight = 23.5; % don't forget to write code to extract head center information from .LOG file!!
% EVERYTHING IS IN METERS!!!!
%% Load MC variables
data = ImportSMR('kl767.smr');
[~,ch] = AddSMRData(data,prs);
f = fopen('MC_Variables767','r') ;
data = dlmread('MC_Variables767');
[mcvariables] = AddMCData(data,prs);

extractMCvariables(mcvariables)
n = length(Moog_TiltX_Pos);

% clear tline
% j = 1 ;
% finish = 0 ;
% while ~finish
%     for i = 1:58
%         str = fgets(f) ;
%         if str==-1, finish=1;end
%         [~,k] = find(str=='=') ;
%         str = [str(1:(k-2)) '[' num2str(j) ',1]' str(k-1:end-2) ';'] ;
%         str = eval(str) ;
%     end
%     fgets(f) ;
%     fgets(f) ;
%     j=j+1 ;
% end
% fclose(f)
% n=j-1 ;

timestamp = (timestamp - timestamp(1))*.001;
flag = flag.*timestamp;
indx = find(flag);
% Moog_Pitch = sind(Moog_TiltX_Pos(1:n).*cos(Moog_Yaw_Pos(1:n))+Moog_TiltY_Pos(1:n).*sin(Moog_Yaw_Pos(1:n)))*981 ;
%% Plots
% gravity components are in moog coordinates, should be transform for subject
% VR Y components are not really centrifugal because they're not in subject's coordinates

% forward
figure;suptitle('Forward');
% plot(gradient(gradient(VR_X_Pos(1:n)))/(dt^2));
% plot(VR_X_Acc(1:n)-gradient(gradient(VR_X_Pos(1:n)))/(dt^2));
subplot(2,1,1);hold on; % VR
plot(timestamp,VR_X_Acc(1:n));  % corrected through the feedback loop
plot(timestamp,gradient(finalVR_X_Vel)/dt);
plot(timestamp,VR_X_GIAerror(1:n));
vline(flag(indx),'k');set(gca,'YLim',[-1 1]);title('VR');xlabel('time (s)');
legend('VR Desired Forward Acceleration','VR Achieved Forward Acceleration','VR forward GIA Error')

subplot(2,1,2);hold on; % Subject
plot(timestamp,Rat_X_Acc); % subject, corrected through the feedback loop
plot(timestamp,Rat_X_GIA); % subject
plot(timestamp,Rat_X_GIAerror); % subject
% plot(sind(Moog_TiltX_Pos)*9.81) % gravity component in subject's forward axis because of tilt
% plot(sind(Moog_TiltX_Acc*HeadHeight*(pi/180))); % acceleration of tilting back
vline(flag(indx),'k');set(gca,'YLim',[-1 1]);title('Subject');xlabel('time (s)');
legend('Subj Desired Forward Acceleration','Subj Achieved Forward Acceleration','Subj forward GIA Error')
%    ...'Gravity Component in Forward Axis','Tilting Back Acceleration');

% centrifugal
figure;suptitle('Centrifugal');
subplot(2,1,1);hold on; % VR
plot(timestamp,VR_Y_Acc(1:n)); % corrected through the feedback loop
plot(timestamp,gradient(finalVR_Y_Vel(1:n))/dt);
plot(timestamp,VR_Y_GIAerror(1:n));
vline(flag(indx),'k');set(gca,'YLim',[-1 1]);title('VR');xlabel('time (s)');
legend('VR Desired Centrifugal Acceleration','VR Achieved Centrifugal Acceleration','VR Centrifugal GIA Error')

subplot(2,1,2);hold on; % Subject
plot(timestamp,Rat_Y_Acc); % subject, corrected through the feedback loop
plot(timestamp,Rat_Y_GIA); % subject
plot(timestamp,Rat_Y_GIAerror); % subject
% plot(sind(Moog_TiltY_Pos)*9.81) % gravity component in subject's forward axis because of tilt
% plot(sind(Moog_TiltY_Acc*HeadHeight*(pi/180))); % acceleration of tilting back
vline(flag(indx),'k');set(gca,'YLim',[-1 1]);title('Subject');xlabel('time (s)');
legend('Subj Desired Centrifugal Acceleration','Subj Achieved Centrifugal GIA',' Subj Centrifugal GIA Error')
%     ...'Gravity Component in Centrifugal Axis','Centrifugal Tilt Acceleration');

%%
X=[];Y=[];

for ii = 1:length(trials)
    try
        Omega = trials(ii).continuous.w ;
        vel = trials(ii).continuous.vyaw-Offset_Acc(3);
        
        X = cat(1,X,[Omega]) ;
        Y = cat(1,Y,[vel*80]) ;
        %
    end
end
trials(ii).prs

x = -50:1:50 ;y =  [] ;
for i = 1:length(x), y(i)=median(Y(abs(X-x(i))<2));end
plot(x,-y,'-o')

axis([-60 60 -60 60]);axis square;grid on
xlabel('Motor command (channel 5) (°/s)');ylabel('Motor velocity (read from sensor) (channel 15) (°/s)')
set(gca,'FontSize',10,'FontWeight','bold')

% title('From Janna data')
set(gcf,'PaperPositionMode','auto','Color','w','InvertHardCopy','off');
%%
JS=[] ;X = [] ; Y = [] ;Z=[] ;A=[];B=[] ;C=[];
figure;
% ii=ii+1
for ii = 1:5:length(trials)
    try
        VX = trials(ii).continuous.v ;
        Omega = trials(ii).continuous.w ;
        acc = [trials(ii).continuous.yac trials(ii).continuous.xac] ;
        vel = [trials(ii).continuous.vyaw trials(ii).continuous.vpit trials(ii).continuous.vrol] ;
        time = (1:length(VX))/833.333 ;
                vx = cumsum(vel) ;
        a = mean(vx(1:500,:));b=mean(vx(end-500:end,:)) ;
        offset = (b-a)/(length(time)-500) ;
        for i = 2:3, vel(:,i)=vel(:,i)-offset(i) ;end
        
        offset = mean([vel(1:500,:) ; vel(end-500:end,:)]) ;
        for i = 1, vel(:,i)=vel(:,i)-offset(i) ;end
        
        
        
        dt = 1/833.333 ;
        Yaw = cumsum(Omega)*dt ;
        V_Allo = [cosd(Yaw).*VX sind(Yaw).*VX] ;
        Pos_Allo = cumsum(V_Allo)*dt ;
        [~,Acc_Allo] = gradient(V_Allo,dt) ;
        Acc_X = cosd(Yaw).*Acc_Allo(:,1)+sind(Yaw).*Acc_Allo(:,2) ;
        Acc_Y = -sind(Yaw).*Acc_Allo(:,1)+cosd(Yaw).*Acc_Allo(:,2) ;
        
        GX = -cumsum(vel(:,2))*80*dt*981/57 ;
        GY = -cumsum(vel(:,3))*80*dt*981/57 ;
        
        acc(:,1)=acc(:,1)-Offset_Acc(1);acc(:,2)=acc(:,2)-Offset_Acc(2);

        subplot(2,4,1)
        plot(Pos_Allo(:,1),Pos_Allo(:,2),'k') ;
        hold on
        plot(Pos_Allo(1:83:end,1),Pos_Allo(1:83:end,2),'ok','MarkerFaceColor','k','MarkerSize',3) ;
        hold off
        set(gca,'CameraUpVector',[1 0 0]) ;
        
        axis equal
        grid on
        
        subplot(2,4,2) ;
        plot(time,[VX cumsum(acc(:,1))*dt*2000]);title('Forwards Vel')
        
        subplot(2,4,3)
        plot(time, [Omega -JLWF(vel(:,1)*80,50)]);title('Angular Vel')
        
        subplot(2,4,5)
        plot(time, [Acc_X JLWF(acc(:,1),300)*2000 GX]);legend('Desired','Actual','Tilt')
        set(gca,'YLim',[-100 100])
        title('Forward Acc')
        
        subplot(2,4,6)
        plot(time, [-Acc_Y JLWF(acc(:,2),300)*2000 GY]);legend('Desired','Actual','Tilt')
        set(gca,'YLim',[-100 100])
        title('Lateral Acc')
        
        X = cat(1,X,[Acc_X Acc_Y]) ;
        Y = cat(1,Y,[GX GY]) ;
        Z = cat(1,Z,acc) ;
        A = cat(1,A,VX) ;
        B = cat(1,B,cumsum(acc(:,1))*dt) ;
        C = cat(1,C,vel*80) ;
        JS = cat(1,JS,VX*0+trials(ii).prs.js_coef) ;
        
        %
        
        Forward(ii,:)=[sum(cumsum(Acc_X)) sum(cumsum(acc(:,1)))*1473]*dt*dt;
    end
end
trials(ii).prs;
%% Regression (Akis' version)
% Z = sind(JLWF([trials.GyroVel]',100)) ; Z = gradient(Z)*833 ;
Y = JLWF(ch.yac,200) ;
X = JLWF(diff(ch.v)/ch.dt,200) ;
JS = [trials.JS]' ;
phase = [trials.phase]' ;

% X=[trials.AlignedAcc] ; X=X(:,1:2:end) ;X=X(:) ;
% Y=[trials.AlignedGyroAcc] ; Y=Y(:,1:2:end) ;Y=Y(:);

lag = 150 ;
Y = Y(lag:end,:);
% Z = Z(lag:end,:);
JS=JS(lag:end,:) ;
phase=phase(lag:end,:) ;

X=X(1:(end-lag+1),:) ;
clf
jsI = [0 0.975 0.99] ;
for j = 1:length(jsI)
    hold on
    x = (-40:2:40)' ;
    y = x*0 ;
    for i = 1:length(x)
        y(i)=nanmean(Y(abs(X(:,1)-x(i))<5 &JS==jsI(j)&~ismember(phase,[1 7]) ,1)) ;
%         I = abs(X(:,1)-x(i))<3 &phase==j ;
%         if sum(I)>1000
%         y(i)=mean(Y(I,1)) ;
%         else
%             y(i)=NaN;
%         end
    end
    b = regress(y,[x x.^0]) 
    
    plot(x,y-b(2),'o-')
end
axis equal
axis([-100 100 -100 100]*0.6);axis square;grid on
plot([-100 100],[-100 100],'k')

xlabel('Visual Acceleration (cm/s^2)')
ylabel('Vestibular Acceleration (cm/s^2)')
legend('Velocity control','Intermediate','Acceleration control')
% legend('Init','Acc','Stop Acc','Plateau','Dec','Stop Dec','Wait')
set(gcf,'Position',[782 734 780 780],'PaperPositionMode','auto','Color','w','InvertHardCopy','off');

%%
[v] = coef2vmax(0.99, 0.05, 15, 20, [], [], 0);v=v*1 ;
Acc = gradient(v',1/60) ;
time60 = (1:length(Acc))/60 ;

[R] = Motion_Cueing_VR(time60,Acc/100,Acc*0,Acc*0,Parameters) ;

figure(3) ;clf

subplot(311)
hold on
plot(time,[nanmean(X,2) nanmean(Y,2) nanmean(Z,2)])
axis([0 5 ylim])
xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
title('Measured')
legend('Desired GIA','Measured G','Measured GIA')
box on; grid on

subplot(312)
plot(time60,[Acc sind(R.Moog_TiltX_Pos)*981 R.Rat_X_GIA*100])
axis([0 5 ylim])
xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
legend('Desired GIA','G','Final GIA')

title('Simulated')
box on; grid on

subplot(313)
plot(time60,[JLWF([R.Moog_X_Acc*100 100*R.Moog_TiltX_Acc*pi/180*1.22 R.finalMoog_X_Acc*100-R.Moog_X_Acc*100],2)])
axis([0 5 ylim])
xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
title('GIA = A+G-(Tangential A)')
box on; grid on
legend('Desired Linear A','Tangential A','Linear A Error (Moog limits)')


%% Gyro...
figure(1);clf
dt = [-300:4000]/833.333 ;

O=[trials.AlignedGyroVel] ; O=O(:,1:3:end) ;
Om=[trials.AlignedOmega] ;
for i = 1:size(O,2), 
    O(:,i)=-O(:,i)*sign(Om(1500,i));
    Om(:,i)=Om(:,i)*sign(Om(1500,i));
end
O = JLWF(O,15) ;
plot([nanmean(O,2) nanmean(Om,2)])
%%
stimtype=gF2(trials,'prs.stimtype') ;

trials_vo = trials(stimtype==2) ;
trials = trials(stimtype~=2) ;
js = gF2(trials,'prs.js_coef') ;
trials = trials(js==0.99) ;

% 
% for ii = 1:length(trials)
%     
% try
% 370 - switch mode
ii=0 ;
%%
stimtype=gF2(trials,'prs.stimtype') ;

trials_vo = trials(stimtype==2) ;
trials = trials(stimtype~=2) ;
x = [] ;
for i = 1:length(trials_vo)
   x = cat(1,x, [trials_vo(i).continuous.yac trials_vo(i).continuous.xac trials_vo(i).continuous.vyaw] ) ;
end
Offset_Acc = median(x) ;

%% Import a smr file
rawSR=833.333 ;
% load('MC_Variables20180925114112.mat')
D = dlmread('kl111.txt','') ;
data = ImportSMR('kl111.smr'); 
default_prs;
[trials,ch] = AddSMRData(data,prs); 




R = struct('X',ch.ymp,'Y',ch.xmp,'TSI',ch.tsi,'V',ch.v,'W',ch.w,'XACC',ch.yac,'YACC',ch.xac) ;
ind = find(ch.tsi>0,1);

timeSMR = (1:length(R(1).X))*prs.dt ;
timeSMR = timeSMR(1:end-ind+1);
timeMD = (1:size(D,1))/58.3 ;
% R(1).X = R(1).X-median(R(1).X(timeSMR<10)) ;
cla; hold on
plot(timeSMR,R(1).X(ind:end)-.2) ;
plot(timeMD, D(:,54)*.18) ;

legend('SMR','MoogDots')
% 
% plot(trials.continuous.ts,trials.continuous.ymp - trials.continuous.ymp(1));hold on;
% plot(timeMD, D(:,54)*.2) ;
% legend('SMR','MoogDots')
