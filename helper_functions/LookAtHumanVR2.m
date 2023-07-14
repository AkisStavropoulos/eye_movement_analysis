default_prs;
trials = AddTrials2Behaviour(prs);
load('C:\Users\Jean\Dropbox\Houston_Matlab_Analysis\Motion Cueing\Human\Offset.mat')

%%
h = waitbar(0,'Calibrating....');
I = [] ;
for ii = 1:length(trials)
    try
        
        waitbar(ii/length(trials),h)
        %%
        VX = trials(ii).continuous.v ;
        Omega = trials(ii).continuous.w ;
        acc = [trials(ii).continuous.yac trials(ii).continuous.xac] ;
        vel = [trials(ii).continuous.vyaw trials(ii).continuous.vpit trials(ii).continuous.vrol] ;
        time = (1:length(VX))/833.333 ;
        Y = -trials(ii).continuous.xmp ;
        X = trials(ii).continuous.ymp ;
        X=X-median(X(1:300));Y=Y-median(Y(1:300));
        
        acc(:,1)=acc(:,1)-Offset_Acc(1);acc(:,2)=-acc(:,2)+Offset_Acc(2);
        
        [b,a]=butter(2,10/(833.333/2),'low') ;
        acc = filter(b,a,acc(end:-1:1,:)) ;
        acc = filter(b,a,acc(end:-1:1,:)) ;

%          X = filter(b,a,X(end:-1:1,:)) ;
%          X = filter(b,a,X(end:-1:1,:)) ;
%          Y = filter(b,a,Y(end:-1:1,:)) ;
%          Y = filter(b,a,Y(end:-1:1,:)) ;
         

          
        vx = cumsum(vel) ;
        a = mean(vx(1:500,:));b=mean(vx(end-500:end,:)) ;
        offset = (b-a)/(length(time)-500) ;
        %     offset=mean(vel(1:500,:)) ;
        for i = 2:3, vel(:,i)=vel(:,i)-offset(i) ;end
      
        [b,a]=butter(2,30/(833.333/2),'low') ;
        vel = filter(b,a,vel(end:-1:1,:)) ;
        vel = filter(b,a,vel(end:-1:1,:)) ;
     
        offset = mean([vel(1:500,:) ; vel(end-500:end,:)]) ;
        for i = 1, vel(:,i)=vel(:,i)-offset(i) ;end
        
        trials(ii).GyroVel = vel'*80 ;
        trials(ii).GyroAcc = acc'*1473 ;
        
        dt = 1/833.333 ;
        Yaw = cumsum(Omega)*dt ;
        V_Allo = [cosd(Yaw).*VX sind(Yaw).*VX] ;
        
        
        Pos_Allo = cumsum(V_Allo)*dt ;
        [~,Acc_Allo] = gradient(V_Allo,dt) ;
        Acc_X = cosd(Yaw).*Acc_Allo(:,1)+sind(Yaw).*Acc_Allo(:,2) ;
        Acc_Y = -sind(Yaw).*Acc_Allo(:,1)+cosd(Yaw).*Acc_Allo(:,2) ;        
        trials(ii).CommandedAcc = [Acc_X' ; Acc_Y'] ;

        Pos_Allo = [X Y] ;[~,V_Allo] = gradient(Pos_Allo,dt) ;V_Allo(end-300:end,:)=0 ;
        [~,Acc_Allo] = gradient(JLWF(V_Allo,300),dt) ;
        Acc_X = cosd(Yaw).*Acc_Allo(:,1)+sind(Yaw).*Acc_Allo(:,2) ;
        Acc_Y = -sind(Yaw).*Acc_Allo(:,1)+cosd(Yaw).*Acc_Allo(:,2) ;
        
        trials(ii).DesiredAcc = [Acc_X' ; Acc_Y'] ;
        
        trials(ii).DesiredAcc = JLWF(trials(ii).DesiredAcc(1,:)',20)'  ;
        
        GX = -sind(cumsum(vel(:,2))*80*dt)*981 ;
        GY = sind(cumsum(vel(:,3))*80*dt)*981 ;
        trials(ii).MoogTilt = [GX'; GY'] ;
        
        trials(ii).JS = zeros(1,length(trials(ii).GyroVel))+trials(ii).prs.js_coef ;
              
        trials(ii).phase = trials(ii).JS*0 ;
%         [~,k] = find(trials(ii).DesiredAcc(1,:)>10,1,'first') ;
%         trials(ii).phase(k:end)=1 ;
%         [~,k] = find(trials(ii).DesiredAcc(1,:)>10,1,'last') ;
%         trials(ii).phase(k:end)=2 ;
%         [~,k] = find(trials(ii).DesiredAcc(1,:)<-10,1,'first') ;
%         trials(ii).phase(k:end)=3 ;
%         [~,k] = find(trials(ii).DesiredAcc(1,:)<-10,1,'last') ;
%         trials(ii).phase(k:end)=4 ;
        n = (1:length(trials(ii).DesiredAcc))' ;T=[] ;
        tmp = JLWF(trials(ii).DesiredAcc(1,:)',100) ;
        k = find(trials(ii).DesiredAcc(1,:)'>10,1,'first') ;
        l = find(trials(ii).DesiredAcc(1,:)'>20,1,'first') ;
        b = regress((k:l)',[trials(ii).DesiredAcc(1,k:l)' trials(ii).DesiredAcc(1,k:l)'.^0]) ;
        T(1) = round(b(2)) ;
        k = find(trials(ii).DesiredAcc(1,:)'<20&n>l,1,'first') ;
        l = find(trials(ii).DesiredAcc(1,:)'<10&n>k,1,'first') ;
        b = regress((k:l)',[trials(ii).DesiredAcc(1,k:l)' trials(ii).DesiredAcc(1,k:l)'.^0]) ;
        T(3) = round(b(2)) ;
        [~,k]=max(tmp(T(1):T(3))) ;
        T(2)=T(1)+k-1 ;
        
         k = find(trials(ii).DesiredAcc(1,:)'<-10,1,'first') ;
        l = find(trials(ii).DesiredAcc(1,:)'<-20,1,'first') ;
        b = regress((k:l)',[trials(ii).DesiredAcc(1,k:l)' trials(ii).DesiredAcc(1,k:l)'.^0]) ;
        T(4) = round(b(2)) ;
        k = find(trials(ii).DesiredAcc(1,:)'>-20&n>l,1,'first') ;
        l = find(trials(ii).DesiredAcc(1,:)'>-10&n>k,1,'first') ;
        b = regress((k:l)',[trials(ii).DesiredAcc(1,k:l)' trials(ii).DesiredAcc(1,k:l)'.^0]) ;
        T(6) = round(b(2)) ;
        [~,k]=min(tmp(T(4):T(6))) ;
        T(5)=T(4)+k-1 ;
        T = [1 T n(end)] ;
        for i = 1:7, trials(ii).phase(T(i):T(i+1))=i;end
        trials(ii).T=T ;

%         trials(ii).DesiredAcc(end-500:end)=NaN ;
%         trials(ii).GyroAcc(:,1:500)=NaN ;
%         trials(ii).GyroAcc(:,end-500:end)=NaN ;
        
        time = (1:length(Yaw))*dt ;
        time60 = (1/60):(1/60):time(end) ;
        trials(ii).SimulateAcc = interp1(time,trials(ii).DesiredAcc',time60) ;
        trials(ii).SimulateYaw = interp1(time,Yaw,time60) ;
        
        [k] = find(trials(ii).DesiredAcc(1,:)'>5,1,'first') ;
        [l] = find(trials(ii).DesiredAcc(1,:)'>20,1,'first') ;
        b = regress((k:l)',[trials(ii).DesiredAcc(1,k:l)' trials(ii).DesiredAcc(1,k:l)'.^0]) ;
        k = round(b(2)) ;Ik=k+(-300:4000) ;
        trials(ii).time = time+trials(ii).events.t_beg ;
        try
            trials(ii).AlignedAcc = trials(ii).DesiredAcc(:,Ik)' ;
            trials(ii).AlignedCommandAcc = trials(ii).CommandedAcc(:,Ik)' ;
            trials(ii).AlignedMoogTilt = [GX(Ik)'; GY(Ik)']' ;
            trials(ii).AlignedGyroAcc = acc(Ik,:)*1473 ;
            trials(ii).AlignedGyroVel = vel(Ik,:)*80 ;
            trials(ii).AlignedOmega = Omega(Ik,:);
        catch
            trials(ii).AlignedAcc = nan(2,length(Ik))' ;
            trials(ii).AlignedCommandAcc = nan(2,length(Ik))' ;
            trials(ii).AlignedMoogTilt = nan(2,length(Ik))' ;
            trials(ii).AlignedGyroAcc = nan(2,length(Ik))' ;
            trials(ii).AlignedGyroVel = nan(3,length(Ik))' ;
            trials(ii).AlignedOmega = nan(1,length(Ik))';
 end
    catch
        I = cat(1,I ,ii) ; 
    end
end
close(h)
% trials(I)=[] ;

%%
cla
timeVR = (1:length(Rat_XY_Acc.x))/58.55 ;
plot([trials.time],[trials.DesiredAcc]')
plot(timeVR,Rat_XY_Acc.x*100)

%% Test for lags
js = gF2(trials,'prs.js_coef') ;

jsI = 0.975 ;
X = [trials(ismember(js,jsI)).DesiredAcc]' ;
Y = JLWF([trials(ismember(js,jsI)).GyroAcc]',300) ;
Z = [trials(ismember(js,jsI)).MoogTilt]' ; 
dt = [-3000:3000]/833.3333 ;
C = xcorr(Z(:,1),X(:,1),3000);
plot(C)
[~,k]=max(C) ;[k-3000 dt(k)]

%% Lag from aligned data
figure(1);clf
js = gF2(trials,'prs.js_coef');unique(js)
T = find(js==0.0) ;
dt = [-300:4000]/833.333 ;
X=[trials(T).AlignedAcc] ; X=X(:,1:2:end) ;
plot(X,'k'); hold on
C=[trials(T).AlignedCommandAcc] ; C=C(:,1:2:end) ;
Y=[trials(T).AlignedMoogTilt] ; Y=Y(:,1:2:end) ;
Z=[trials(T).AlignedGyroAcc] ; Z=Z(:,1:2:end)-25 ;
O=[trials(T).AlignedGyroVel] ; O=-O(:,2:3:end) ;

I = find(X(1000,:)>30&X(1000,:)<80&X(2000,:)>10&X(2000,:)<25) ;
I = find(X(1000,:)>30&X(1000,:)<60) ;
% I = find(X(800,:)>25&X(760,:)<45) ;
I = find(X(900,:)>0) ;
size(I)
% X=X(:,I);Y=Y(:,I);Z=Z(:,I);O=O(:,I);C=C(:,I) ;
ylim = [-50 60] ;
plot(X,'r')
% plot(dt,[nanmean(X,2) nanmean(Y,2) nanmean(Z,2) nanmedian(O,2) -JLWF(gradient(nanmedian(O,2),1/833),100)*150/57])

time = dt-dt(1) ;
time60 = (1/60):(1/60):time(end) ;
Acc = interp1(time,nanmean(C,2),time60)' ;

figure(2);clf
[R] = Motion_Cueing_VR(time60,Acc/100,Acc*0,Acc*0) ;

figure(3) ;clf

subplot(211)
hold on
plot(time,[nanmean(C,2) nanmean(X,2) nanmean(Z,2)])
% plot(time,[nanmean(C,2) nanmean(X,2) nanmean(Z,2) nanmean(Y,2)])
axis([0 5 ylim])
xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
title('Measured')
legend('Joystick','Visual','Vestibular')
box on; grid on

subplot(212)
% plot(time60,[Acc R.Rat_X_GIA*100 Acc*NaN sind(R.Moog_TiltX_Pos)*981])
plot(time60,[Acc R.Rat_X_GIA*100 Acc*NaN])
axis([0 5 ylim])
xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
legend('Joystick','Visual','Vestibular')
title('Simulated')
box on; grid on

% subplot(313)
% plot(time60,[JLWF([R.Moog_X_Acc*100 100*R.Moog_TiltX_Acc*pi/180*1.22 R.finalMoog_X_Acc*100-R.Moog_X_Acc*100],2)])
% axis([0 5 ylim])
% xlabel('time (s)');ylabel('Acc./G/GIA (cm/s^2)')
% title('GIA = A+G-(Tangential A)')
% box on; grid on
% legend('Desired Linear A','Tangential A','Linear A Error (Moog limits)')
% subplot(313)
% plot(time60,[(cumsum(Acc))/(60) (cumsum(R.Rat_X_GIA*100))/(60)])
% axis([0 5 0 50]);xlabel('time (s)');ylabel('Total displacement (cm)')

%% Load MC variables
f = fopen('MC_Variables3','r') ;
clear tline
j = 1 ;
finish = 0 ;
while ~finish
    for i = 1:56
        str = fgets(f) ;
        if str==-1, finish=1;end
        [~,k] = find(str=='=') ;
        str = [str(1:(k-2)) '(' num2str(j) ',1)' str(k-1:end-2) ';'] ;
        eval(str) ;
    end
    fgets(f) ;
    fgets(f) ;
    j=j+1 
end
fclose(f)
n=j-1 ;
Moog_Pitch = sind(Moog_TiltXY_Pos.x(1:n).*cos(Moog_Yaw_Pos(1:n))+Moog_TiltXY_Pos.y(1:n).*sin(Moog_Yaw_Pos(1:n)))*981 ; 


plot([gradient(gradient(VR_XY_Pos.x(1:n)))*3600 gradient(VR_XY_Vel.x(1:n))*60 ...
      VR_XY_Acc.x(1:n) VR_XY_Acc.x(1:n)-gradient(gradient(VR_XY_Pos.x(1:n)))*3600 VR_XY_GIAerror.x(1:n)])
  set(gca,'YLim',[-1 1])

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

%% Regression
Z = sind(JLWF([trials.GyroVel]',100)) ; Z = gradient(Z)*833 ;
Y = JLWF([trials.GyroAcc]',200) ;
X = JLWF([trials.DesiredAcc]',200) ;
JS = [trials.JS]' ;
phase = [trials.phase]' ;

% X=[trials.AlignedAcc] ; X=X(:,1:2:end) ;X=X(:) ;
% Y=[trials.AlignedGyroAcc] ; Y=Y(:,1:2:end) ;Y=Y(:);

lag = 150 ;
Y = Y(lag:end,:);Z = Z(lag:end,:);JS=JS(lag:end,:) ;phase=phase(lag:end,:) ;

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

title('From Janna data')
set(gcf,'Position',[849 910 711 588],'PaperPositionMode','auto','Color','w','InvertHardCopy','off');

%%
JS=[] ;X = [] ; Y = [] ;Z=[] ;A=[];B=[] ;C=[];

ii=ii+1
for ii = 1:5:length(trials)
    try
        VX = trials(ii).continuous.v ;
        Omega = trials(ii).continuous.w ;
        acc = [trials(ii).continuous.yac trials(ii).continuous.xac] ;
        vel = [trials(ii).continuous.vyaw trials(ii).continuous.vpit trials(ii).continuous.vrol] ;
        time = (1:length(VX))/833.333 ;
        
        % vx = cumsum(acc) ;
        % a = mean(vx(1:500,:));b=mean(vx(end-500:end,:)) ;
        % offset = (b-a)/(length(time)-500) ;
        acc(:,1)=acc(:,1)-Offset_Acc(1);acc(:,2)=acc(:,2)-Offset_Acc(2);
        
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
trials(ii).prs
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

%% test Moog Commands
%% ONLY moog commands
    moogpos = dlmread('moogTraj525');
    label = [{'forward'} {'lateral'} {'Z'} {'yaw'} {'pitch'} {'roll'}];
    ts = moogpos(:,1); ts = ts-ts(1); ts = ts+ts(2); % in miliseconds
    ts = ts*.001;

    diffts = [0; diff(ts)];
    ind = find(diffts > .02); % find indices of sample difference bigger than 20 ms
    lostsignal = diffts(ind);
    
    moogvel = diff(moogpos,1,1)/prs.dt; moogvel = [zeros(1,size(moogpos,2)); moogvel];
    moogacc = diff(moogpos,2,1)/prs.dt; moogacc = [zeros(2,size(moogpos,2)); moogacc];
    moogjerk = diff(moogpos,3,1)/prs.dt; moogjerk = [zeros(3,size(moogpos,2)); moogjerk];

    
    figure;
        subplot(4,1,1);plot(ts,moogpos(:,2));hold on;
        plot(ts(ind),lostsignal*2,'.');
        hline(0);xlabel('time (s)');ylabel('volts');
  title('forward position of Moog');
    
        subplot(4,1,2);plot(ts,moogvel(:,2));hold on;
        plot(ts(ind),lostsignal*2,'.');
        hline(0);xlabel('time (s)');ylabel('volts');
        
title('forward velocity of Moog');

        subplot(4,1,3);plot(ts,moogacc(:,2));hold on;
        plot(ts(ind),lostsignal,'.');        
        hline(0);xlabel('time (s)');ylabel('volts');
title('forward acceleration of Moog');
    
        subplot(4,1,4);
        plot(ts,moogjerk(:,2));hold on;
        plot(ts(ind),lostsignal,'.');        
        hline(0);xlabel('time (s)');ylabel('volts');
title('forward jerk of Moog');

    %% Moog commands and SMR file
data = ImportSMR('kl513.smr'); 
[trials,ch] = AddSMRData(data,prs); 
timeSMR = [1:length(ch.tsi)]*prs.dt;

if size(moogpos,2)==6 % without timestamp
    
    moogpos = dlmread('moogTraj513');
    label = [{'forward'} {'lateral'} {'Z'} {'yaw'} {'pitch'} {'roll'}];
    figure;
    for i = 1:size(moogpos,2)
        subplot(size(moogpos,2),1,i);plot(moogpos(:,i));hold on;
        title(label{i});
    end
    suptitle('position of Moog');
    
    moogvel = diff(moogpos,1,1)/prs.dt; moogvel = [zeros(1,size(moogpos,2)); moogvel];
    moogacc = diff(moogpos,2,1)/prs.dt; moogacc = [zeros(2,size(moogpos,2)); moogacc];
    moogjerk = diff(moogpos,3,1)/prs.dt; moogjerk = [zeros(3,size(moogpos,2)); moogjerk];

    figure;
    for i = 1:size(moogvel,2)
        subplot(size(moogvel,2),1,i);plot(moogvel(:,i));hold on;
        title(label{i});
    end
    suptitle('velocity of Moog');
    figure;
    for i = 1:size(moogacc,2)
        subplot(size(moogacc,2),1,i);plot(moogacc(:,i));hold on;
        title(label{i});
    end
    suptitle('acceleration of Moog');
    
    for j = 54:57
        D = dlmread('MC_Variables513','') ;
        figure;plot(D(:,j)) ; title(num2str(j))
    end
    
elseif size(moogpos,2)==7 % with timesamp
    
    
    moogpos = dlmread('moogTraj513');
    label = [{'forward'} {'lateral'} {'Z'} {'yaw'} {'pitch'} {'roll'}];
    ts = moogpos(:,1); ts = ts-ts(1); ts = ts+ts(2); % in miliseconds
    ts = ts*.001;

    diffts = [0; diff(ts)];
    ind = find(diffts > .02); % find indices of sample difference bigger than 20 ms
    lostsignal = diffts(ind);

    
    figure;
    for i = 2:size(moogpos,2) % skip the first column, the timestamp
        subplot(size(moogpos,2)-1,1,i-1);plot(ts,moogpos(:,i));title(label{i-1});hold on;
        plot(ts(ind),lostsignal*2,'.');
        plot(timeSMR,ch.tsi/3);hline(0);
    end
    suptitle('position of Moog');
    
    moogvel = diff(moogpos,1,1)/prs.dt; moogvel = [zeros(1,size(moogpos,2)); moogvel];
    moogacc = diff(moogpos,2,1)/prs.dt; moogacc = [zeros(2,size(moogpos,2)); moogacc];
    figure;
    for i = 2:size(moogvel,2)
        subplot(size(moogvel,2)-1,1,i-1);plot(ts,moogvel(:,i));title(label{i-1});hold on;
        plot(ts(ind),lostsignal*2,'.');
        plot(timeSMR,ch.tsi/3);hline(0);
        
    end
    suptitle('velocity of Moog');
    figure;
    for i = 2:size(moogacc,2)
        subplot(size(moogacc,2)-1,1,i-1);plot(ts,moogacc(:,i));title(label{i-1});hold on;
        plot(ts(ind),lostsignal,'.');        
        plot(timeSMR,ch.tsi/3);hline(0);        
    end
    suptitle('acceleration of Moog');
    
    for j = 54:57
        D = dlmread('MC_Variables513','') ;
        figure;plot(D(:,j)) ; title(num2str(j))
    end
end
%% Plot only for position and acceleration
data = ImportSMR('kl527.smr'); 
[trials,ch] = AddSMRData(data,prs);
timeSMR = [1:length(ch.tsi)]*prs.dt;

    diffts = [0; diff(ts)];
    ind = find(diffts > .02); % find indices of sample difference bigger than 20 ms
    lostsignal = diffts(ind);

moogpos = dlmread('moogTraj527');
label = [{'forward'} {'lateral'} {'Z'} {'yaw'} {'pitch'} {'roll'}];
ts = moogpos(:,1); ts = ts-ts(1); ts = ts+ts(2); % in miliseconds
ts = ts*.001;


moogvel = diff(moogpos,1,1)/prs.dt; moogvel = [zeros(1,size(moogpos,2)); moogvel];
moogacc = diff(moogpos,2,1)/prs.dt; moogacc = [zeros(2,size(moogpos,2)); moogacc];
moogjerk = diff(moogpos,3,1)/prs.dt; moogjerk = [zeros(3,size(moogpos,2)); moogjerk];
% position
figure;
plot(ts,moogpos(:,2));title(label{2-1});hold on;
plot(ts(ind),lostsignal,'.');
plot(timeSMR,ch.tsi/5);hline(0);
plot(timeSMR,(ch.ymp-ch.ymp(end))/1000);
legend('moog position command','time delays','button pushes','fwd position in SMR file');ylabel('volts');xlabel('time (s)');
suptitle('position of Moog');
% acceleration

acc = [0; diff(movmean(movmean(ch.v,20),20))/(prs.dt)];
figure;
plot(ts,moogacc(:,2));title(label{2-1});hold on;
plot(ts(ind),lostsignal,'.');
plot(timeSMR,ch.tsi/4);hline(0);
plot(timeSMR,movmean(acc,20)/200);
legend('moog acceleration command','time delays','button pushes','fwd acceleration in SMR file');ylabel('volts');xlabel('time (s)');
suptitle('acceleration of Moog');
% velocity
figure;
plot(ts,moogvel(:,2));title(label{2-1});hold on;
plot(ts(ind),lostsignal,'.');
plot(timeSMR,ch.tsi/2);hline(0);
plot(timeSMR,ch.v/20);
legend('moog velocity command','time delays','button pushes','fwd velocity in SMR file');ylabel('volts');xlabel('time (s)');
suptitle('velocity of Moog');
% jerk
figure;
plot(ts,moogjerk(:,2));hold on;
plot(ts(ind),lostsignal,'.');
plot(timeSMR,ch.tsi/20);hline(0);
hline(0);xlabel('time (s)');ylabel('volts');legend('moog jerk command','time delays','button pushes');
title('forward jerk of Moog');
suptitle('jerk of Moog');

