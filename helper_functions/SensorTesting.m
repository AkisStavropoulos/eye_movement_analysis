 

stimtype=gF2(trials,'prs.stimtype') ;


trials_vo = trials(stimtype==2) ;

trials = trials(stimtype~=2) ;


 

% X = [] ; Y = [] ;Z=[] ;A=[];B=[] ;

%

% for ii = 1:length(trials)

%    

% try

% 370 - switch mode
%%

x = [] ;

for i = 1:length(trials_vo)

   x = cat(1,x, [trials_vo(i).continuous.yac trials_vo(i).continuous.xac] ) ;

end

Offset_Acc = median(x) ;
figure;
%%

ii = ii+1 ;

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

   

GX = -cumsum(vel(:,2))*85*dt*981/57 ; % pitch

GY = -cumsum(vel(:,3))*85*dt*981/57 ; % roll



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

plot(time, [Omega -JLWF(vel(:,1)*85,50)]);title('Angular Vel')


subplot(2,4,5)

plot(time, [Acc_X JLWF(acc(:,1),300)*2000 GX]);legend('Desired','Actual','Tilt')

set(gca,'YLim',[-100 100])

title('Forward Acc')


subplot(2,4,6)

plot(time, [Acc_Y JLWF(acc(:,2),300)*2000 GY]);legend('Desired','Actual','Tilt')

set(gca,'YLim',[-100 100])

title('Lateral Acc')


% X = cat(1,X,[Acc_X Acc_Y]) ;

% Y = cat(1,Y,[GX GY]) ;

% Z = cat(1,Z,acc) ;

% A = cat(1,A,VX) ;

% B = cat(1,B,cumsum(acc(:,1))*dt) ;

%

% end

% end

trials(ii).prs


