%% check sensor integrals over all trials and correlations
roll_int = [];
pitch_int = [];
yaw_int = [];
xac_int = [];
yac_int = [];
trial_ind = 1;
for i = 1:length(t.beg)
        indx = find(ts>=t.beg(i) & ts<=t.end(i));
        tempr = cumsum(ch.vrol(indx) - mean(ch.vrol(indx)))*dt;
        roll_int = [roll_int; tempr];
        tempp = cumsum(ch.vpit(indx) - mean(ch.vpit(indx)))*dt;
        pitch_int = [pitch_int; tempp];
        tempy = cumsum(ch.vyaw(indx) - mean(ch.vyaw(indx)))*dt;
        yaw_int = [yaw_int; tempy];
        tempya = cumsum(ch.yac(indx) - mean(ch.yac(indx)))*dt;
        yac_int = [yac_int; tempya];
        tempxa = cumsum(ch.xac(indx) - mean(ch.xac(indx)))*dt;
        xac_int = [xac_int; tempxa];        
        ind = trial_ind(end) + length(indx);
        trial_ind = [trial_ind; ind];
end
tt = (1:length(roll_int))*dt;
figure;plot(tt,roll_int,'r');hold on;plot(tt,pitch_int,'b');...
    title('integral of pitch (b) and roll (r) velocity sensors');xlabel('time (s)'); vline(trial_ind*dt);
figure;plot(tt,roll_int,'r');hold on;plot(tt,yaw_int,'k');...
    title('integral of yaw (k) and roll (r) velocity sensors');xlabel('time (s)'); vline(trial_ind*dt);
figure;plot(tt,roll_int,'r');hold on;plot(tt,yac_int,'b');...
    title('integral of y acceleration (k) and roll (r) velocity sensors');xlabel('time (s)'); vline(trial_ind*dt);
figure;plot(tt,roll_int,'r');hold on;plot(tt,xac_int,'k');...
    title('integral of x acceleration (k) and roll (r) velocity sensors');xlabel('time (s)'); vline(trial_ind*dt);
figure;plot(tt,yaw_int,'r');hold on;plot(tt,xac_int,'k');...
    title('integral of y acceleration (k) and yaw (r) velocity sensors');xlabel('time (s)'); vline(trial_ind*dt);


ro_rp = corrcoef(roll_int,pitch_int);
ro_ry = corrcoef(roll_int,yaw_int);
    
ro_rya = corrcoef(roll_int,yac_int);
ro_rxa = corrcoef(roll_int,xac_int);

ro_pya = corrcoef(pitch_int,yac_int);        
ro_xaya = corrcoef(xac_int,yac_int);        
ro_yxa = corrcoef(xac_int,yaw_int);   
ro_yya = corrcoef(yac_int,yaw_int);   

yaw_pitch = pitch_int.*yaw_int;
ro_ypr = corrcoef(yaw_pitch,roll_int);

    
    
    
    
    
    
