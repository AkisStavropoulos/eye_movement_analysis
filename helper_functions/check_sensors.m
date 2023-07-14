%% check sensors
for i = 15:18
    indx = find(ts>=t.beg(i) & ts<=t.end(i));
    tt = ts(indx);
    figure;
    subplot(10,1,1);
    plot(tt,ch.xmp(indx));hold on;plot(tt,ch.ymp(indx));hold off;title('x and y position');
    subplot(10,1,2);
    plot(tt,ch.phi(indx));title('theta (deg)');
    subplot(10,1,3);
    plot(tt,ch.v(indx));title('velocity (cm/s)');
    subplot(10,1,4);
    plot(tt,ch.w(indx));title('angular velocity (deg/s)');
    subplot(10,1,5);
    plot(tt,ch.xac(indx));title('x-axis acceleration sensor reading');
    subplot(10,1,6);
    plot(tt,ch.yac(indx));title('y-axis acceleration sensor reading');
    subplot(10,1,7);
    plot(tt,ch.vrol(indx));title('roll velocity sensor reading');
    subplot(10,1,8);
    plot(tt,ch.vyaw(indx));title('yaw velocity sensor reading');
    subplot(10,1,9);
    plot(tt,ch.vpit(indx));title('pitch velocity sensor reading');xlabel('time(s)');
%     subplot(10,1,10);
%     plot(tt,cumsum(ch.yac(indx) - mean(ch.yac(indx)))*dt);title('integral of a_y');xlabel('time(s)');
    subplot(10,1,10); hold on;
    plot(tt,cumsum(ch.vpit(indx) - mean(ch.vpit(indx)))*dt);title('integral of pitch');xlabel('time(s)');
    plot(tt,cumsum(ch.vrol(indx) - mean(ch.vrol(indx)))*dt);title('integral of roll');xlabel('time(s)');
%     plot(tt,cumsum(ch.yac(indx) - mean(ch.yac(indx)))*dt);title('integral of a_y');xlabel('time(s)');
%     plot(tt,cumsum(ch.xac(indx) - mean(ch.xac(indx)))*dt);title('integral of a_x');xlabel('time(s)');
end


%%