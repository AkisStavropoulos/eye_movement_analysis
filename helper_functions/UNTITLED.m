close all

for n = 1:10
h(n) = figure;
end

for i = 1:length(trials)


% figure(h(1));hold on;plot(trials(i).mc.timestamp);legend('timestamp')
% 
% figure(h(2));hold on;plot(trials(i).mc.flag);legend('flag')

figure(h(3));hold on;plot(trials(i).mc.JS_Yaw_Raw,-trials(i).mc.JS_X_Raw,'.');legend('joystick x raw','joystick yaw raw')

figure(h(4));hold on;plot(trials(i).mc.JS_Yaw_Raw,'.');legend('joystick yaw raw')

figure(h(5));hold on;plot(trials(i).mc.JS_X_Converted,'.');legend('joystick x converted')

figure(h(6));hold on;plot(trials(i).mc.JS_Yaw_Converted,'.');legend('joystick yaw converted')

% figure(h(7));hold on;plot(trials(i).mc.VR_Yaw_Vel,'k');legend('vr yaw velocity')
% 
% figure(h(8));hold on;plot(trials(i).mc.VR_Yaw_Pos,'k');legend('vr yaw position')
% 
% figure(h(9));hold on;plot(trials(i).mc.VR_X_Vel,'k');legend('vr x velocity')
% 
% figure(h(10));hold on;plot(trials(i).mc.VR_Y_Vel,'k');legend('vr y velocity')


end


for i = 14
figure;plot(data(:,i));title(['column ' num2str(i)]);
end

hold on;plot(100*data(:,58)),title([' column ' num2str(7) ' , column ' num2str(58)])
figure;plot(data(:,4),-data(:,3));

%% check random walk

for i = 1:length(trials)
    wmax(i) = trials(i).prs.wmax;
    vmax(i) = trials(i).prs.vmax;
    tau(i) = trials(i).prs.tau;
end
figure;plot(tau);
%% check GIA error and moog position / plot trajectories
indx_err = [];
for j = 1:length(trials)
    if find(trials(j).mc.Rat_X_GIAerror > .1,1)
        indx_err = [indx_err j];
    end
end
    
for i = indx_err
    dist = [];
    target_dist = [];
    moog_dist = [];
    moog_tilt = [];
    target_dist = .01*sqrt((trials(i).prs.fireflyposx-trials(i).continuous.xmp(1)).^2 + (trials(i).prs.fireflyposy-trials(i).continuous.ymp(1)).^2);
    dist = .01*sqrt((trials(i).continuous.ymp-trials(i).continuous.ymp(1)).^2 + (trials(i).continuous.xmp-trials(i).continuous.xmp(1)).^2);
    moog_disp = sqrt(trials(i).mc.Moog_X_Pos.^2 + trials(i).mc.Moog_Y_Pos.^2);
    moog_tilt = sqrt(trials(i).mc.Moog_TiltX_Pos.^2 + trials(i).mc.Moog_TiltY_Pos.^2);
    
    figure;subplot(1,2,1);plot(trials(i).continuous.ts,dist);hold on;
    plot(trials(i).mc.timestamp,trials(i).mc.Moog_X_Pos*10);
    plot(trials(i).mc.timestamp,trials(i).mc.Moog_Y_Pos*10);
    plot(trials(i).mc.timestamp,trials(i).mc.Rat_X_GIAerror*10);
    plot(trials(i).mc.timestamp,moog_disp*10);
    plot(trials(i).mc.timestamp,moog_tilt);
    hline(target_dist);
    legend('distance covered','moog pos X * 10','moog pos Y * 10','GIA error * 10',...
        'total moog displacement * 10','total moog tilt (deg?)','Location','northwest')
    title(['trial No. ' num2str(i)]);xlabel('time (s)');ylabel('distance (m)');grid on;
    
    subplot(1,2,2);plot(trials(i).continuous.xmp*.01,trials(i).continuous.ymp*.01);axis equal;
hold on;plot(trials(i).prs.fireflyposx*.01,trials(i).prs.fireflyposy*.01,'r*');title(['trajectory of trial No. ' num2str(i)]);
xlabel(' x axis (m)');ylabel('y axis-forward (m)')
end
%% plot GIA error as a function of distance
dist = [];
GIAerror = [];
figure;
for i = 1:length(trials)
        dist_temp = [];
        dist_temp = sqrt((trials(i).continuous.ymp-trials(i).continuous.ymp(1)).^2 + (trials(i).continuous.xmp-trials(i).continuous.xmp(1)).^2);
        dist(i) = dist_temp(end);
        GIAerror(i) = 100*max(trials(i).mc.Rat_X_GIAerror);
        hold on;bar(dist(i),GIAerror(i),8)
end
grid on;ylabel('max GIA error of trial (cm/s^2)');xlabel('distance per trial (cm)');
title(['GIA error as a function of distance for x = ' num2str(trials(i).prs.x) ' m, T = ' num2str(trials(i).prs.T) ' s']);

%% average max GIA error
medianGIAerr = median(GIAerror);
avgGIAerr = mean(GIAerror);
%% average target distance
for i = 1:length(trials)
    target_dist(i) = sqrt((trials(i).prs.fireflyposx-trials(i).continuous.xmp(1)).^2 + (trials(i).prs.fireflyposy-trials(i).continuous.ymp(1)).^2);
    trial_dur(i) = trials(i).continuous.ts(end);
end
figure;plot(target_dist,trial_dur,'ko');xlabel('target distance (cm)');ylabel('trial duration (s)');
ylim([0 22]);title('duration over distance');
avgtargetdist = mean(target_dist);
