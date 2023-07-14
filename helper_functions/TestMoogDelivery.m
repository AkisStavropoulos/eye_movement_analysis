%% Test Moog and sensors
default_prs;
dt = 1/60; % 60 Hz
HeadHeight = 23.5; % don't forget to write code to extract head center information from .LOG file!!
% EVERYTHING IS IN METERS!!!!
%% Load MC variables
filenum = input('Choose file number ending: ');
data = ImportSMR(['kl' num2str(filenum) '.smr']);
[~,ch] = AddSMRData(data,prs);
f = fopen(['MC_Variables' num2str(filenum)],'r') ;
data_mc = dlmread(['MC_Variables' num2str(filenum)]);
[~,mc] = AddMCData(data_mc,[]);

% extractMCvariables(mcvariables)

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
% Moog_Pitch = sind(Moog_TiltX_Pos(1:n).*cos(Moog_Yaw_Pos(1:n))+Moog_TiltY_Pos(1:n).*sin(Moog_Yaw_Pos(1:n)))*981 ;
%% Plots
% gravity components are in moog coordinates, should be transform for subject
% VR Y components are not really centrifugal because they're not in subject's coordinates
n = length(mc.Moog_TiltX_Pos);
indx = find(mc.flag);

distance_traveled = sqrt((ch.xmp - ch.xmp(indx(1))).^2 + (ch.ymp - ch.ymp(indx(1))).^2);

% forward
figure;suptitle(['Forward file ' num2str(filenum)]);
% plot(gradient(gradient(VR_X_Pos(1:n)))/(dt^2));
% plot(VR_X_Acc(1:n)-gradient(gradient(VR_X_Pos(1:n)))/(dt^2));
subplot(2,1,1);hold on; % VR
plot(mc.timestamp,mc.VR_X_Acc(1:n));  % corrected through the feedback loop
plot(mc.timestamp,gradient(mc.finalVR_X_Vel)/dt);
plot(mc.timestamp,mc.VR_X_GIAerror(1:n));
% plot(ch.ts,ch.ymp/1000); % distance in decameters
% plot(ch.ts,ch.xmp/1000); % distance in decameters
plot(ch.ts,distance_traveled/1000); % in decameters
vline(mc.flag(indx),'k');set(gca,'YLim',[-1 1]);title('VR');xlabel('time (s)');
legend('VR Desired Forward Acceleration','VR Achieved Forward Acceleration','VR forward GIA Error','distance covered (in dm)')

subplot(2,1,2);hold on; % Subject
plot(mc.timestamp,mc.Rat_X_Acc); % subject, corrected through the feedback loop
plot(mc.timestamp,mc.Rat_X_GIA); % subject
plot(mc.timestamp,mc.Rat_X_GIAerror); % subject
% plot(sind(Moog_TiltX_Pos)*9.81) % gravity component in subject's forward axis because of tilt
% plot(sind(Moog_TiltX_Acc*HeadHeight*(pi/180))); % acceleration of tilting back
vline(mc.flag(indx),'k');set(gca,'YLim',[-1 1]);title('Subject');xlabel('time (s)');
legend('Subj Desired Forward Acceleration','Subj Achieved Forward Acceleration','Subj forward GIA Error')
%    ...'Gravity Component in Forward Axis','Tilting Back Acceleration');

% centrifugal
figure;suptitle(['Centrifugal file ' num2str(filenum)]);
subplot(2,1,1);hold on; % VR
plot(mc.timestamp,mc.VR_Y_Acc(1:n)); % corrected through the feedback loop
plot(mc.timestamp,gradient(mc.finalVR_Y_Vel(1:n))/dt);
plot(mc.timestamp,mc.VR_Y_GIAerror(1:n));
vline(mc.flag(indx),'k');set(gca,'YLim',[-1 1]);title('VR');xlabel('time (s)');
legend('VR Desired Centrifugal Acceleration','VR Achieved Centrifugal Acceleration','VR Centrifugal GIA Error')

subplot(2,1,2);hold on; % Subject
plot(mc.timestamp,mc.Rat_Y_Acc); % subject, corrected through the feedback loop
plot(mc.timestamp,mc.Rat_Y_GIA); % subject
plot(mc.timestamp,mc.Rat_Y_GIAerror); % subject
% plot(sind(Moog_TiltY_Pos)*9.81) % gravity component in subject's forward axis because of tilt
% plot(sind(Moog_TiltY_Acc*HeadHeight*(pi/180))); % acceleration of tilting back
vline(mc.flag(indx),'k');set(gca,'YLim',[-1 1]);title('Subject');xlabel('time (s)');
legend('Subj Desired Centrifugal Acceleration','Subj Achieved Centrifugal GIA',' Subj Centrifugal GIA Error')
%     ...'Gravity Component in Centrifugal Axis','Centrifugal Tilt Acceleration');
%% sensors' filtering
% % acceleration sensors
% offsety = mean(ch.yac);
% offsetx = mean(ch.xac);
% [b,a]=butter(2,10/((1/ch.dt)/2),'low') ;
% yac = filter(b,a,[offsety*ones(50,1); ch.yac(end:-1:1,:); offsety*ones(50,1)]) ;
% xac = filter(b,a,[offsetx*ones(50,1); ch.xac(end:-1:1,:); offsetx*ones(50,1)]) ;
% yac = yac(51:end-50);
% xac = xac(51:end-50);
% yac = (yac-offsety)*1473; % to get cm/s^2
% xac = (xac-offsetx)*1473; % to get cm/s^2
% 
% % rotation sensors
% offsetroll = mean(ch.vrol);
% offsetyaw = mean(ch.vyaw);
% offsetpit = mean(ch.vpit);
% [b,a]=butter(2,10/(833.333/2),'low') ;
% roll = filter(b,a,[offsetroll*ones(20,1); ch.vrol(end:-1:1,:); offsetroll*ones(20,1)]) ;
% yaw = filter(b,a,[offsetyaw*ones(20,1); ch.vyaw(end:-1:1,:); offsetyaw*ones(20,1)]) ;
% pitch = filter(b,a,[offsetpit*ones(20,1); ch.vpit(end:-1:1,:); offsetpit*ones(20,1)]) ;
% roll = roll(21:end-20);
% yaw = yaw(21:end-20);
% pitch = pitch(21:end-20);
% roll = (roll-offsetroll)*80; % to get deg/s^2
% yaw = (yaw-offsetyaw)*80; % to get deg/s^2
% pitch = (pitch-offsetpit)*80; % to get deg/s^2
% 
%% fourrier analysis of sensor signals
% % Y = fft(X);
% % P2 = abs(Y/n);
% % P1 = P2(1:n/2+1);
% % P1(2:end-1) = 2*P1(2:end-1);
% % f = Fs*(0:(L/2))/L;
% % figure;plot(f,P1)
% 
% % fourrier of unfiltered signals
% offsets = [mean(ch.yac);mean(ch.xac);mean(ch.vrol);mean(ch.vyaw);mean(ch.vpit)];
% n = length(ch.yac);
% L = length(ch.yac);
% Fs = 1/ch.dt;
% f = Fs*(0:(n/2))/n;
% sensornames = [{'YAcc'} {'XAcc'} {'VRoll'} {'VYaw'} {'VPitch'}];
% sensignal = [ch.yac';ch.xac';ch.vrol';ch.vyaw';ch.vpit'] - offsets;
% Y = fft(sensignal,n,2);
% P2 = abs(Y/n);
% P1 = P2(:,1:n/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
% figure;
% for i=1:size(Y,1)
%     subplot(size(Y,1),1,i)
%     plot(f,P1(i,:),'r')
%     title(['Sensor ',sensornames{i}, ' in the Frequency Domain']);xlim([0 300]);ylim([0 .02])
% end
% 
% % fourrier of filtered signals
% offsets = [mean(yac);mean(xac);mean(roll);mean(yaw);mean(pitch)];
% n = length(ch.yac);
% L = length(ch.yac);
% Fs = 1/ch.dt;
% f = Fs*(0:(n/2))/n;
% sensornames = [{'YAcc'} {'XAcc'} {'VRoll'} {'VYaw'} {'VPitch'}];
% sensignal = [yac';xac';roll';yaw';pitch'];
% Y = fft(sensignal,n,2);
% P2 = abs(Y./n);
% P1 = P2(:,1:n/2+1);
% P1(:,2:end-1) = 2*P1(:,2:end-1);
% 
% figure;
% for i=1:size(Y,1)
%     subplot(size(Y,1),1,i)
%     plot(f,P1(i,:),'r')
%     title(['Filtered sensor ',sensornames{i}, ' in the Frequency Domain']);xlim([0 30]);ylim([0 10])
% end
% 
% % % fourrier for signal before start of experiment
% % fftlength = 15/ch.dt;
% % n = fftlength;
% % offset = mean(ch.yac(1:n));
% % Y = fft(ch.yac(1:n) - offset);
% % P2 = abs(Y/n);
% % P1 = P2(1:n/2+1);
% % P1(2:end-1) = 2*P1(2:end-1);
% % f = Fs*(0:(n/2))/n;
% % figure;plot(f,P1);title('spectrum of unfiltered "yac - offset", before start') % I get 2 frequencies: 23 and 270 Hz
%% fit sine wave
% signal = yac;
% indx = find(P1(1,:)>1,1,'last');
% fr = f(indx);
% t = 1:length(yac);
% sinclr = @(A)(signal - A*sin(2*pi*fr*t));
% A0 = 1; phase0 = .1;
% params0 = [A0,phase0];
% [params,fval] = fminsearch(sinclr,A0);
% 
% 
% 
