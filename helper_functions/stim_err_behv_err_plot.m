function stim_err_behv_err_plot(trials,stimtype)
%% Plot final distance from target as a function of displacement due to GIAerror
sanity_check = 0;
dt = trials(1).mc.timestamp(2) - trials(1).mc.timestamp(1);
ind_rem = [];
for i = 1:length(trials)
    GIAerror{i} = [trials(i).mc.VR_X_GIAerror trials(i).mc.VR_Y_GIAerror];
    disp_temp = cumsum(cumsum(GIAerror{i},1)*dt,1)*dt;
    if ~isempty(disp_temp)
        stim_errX(i) = disp_temp(end,1);
        stim_errY(i) = disp_temp(end,2);
        disp_temp = sqrt(disp_temp(:,1).^2 + disp_temp(:,1).^2);
        stim_err(i) = disp_temp(end);
    else
        ind_rem = [ind_rem i];
    end
    
end
non_visual = sort([stimtype.s1 stimtype.s3]);
if ~isempty(ind_rem)
    ind_rem = intersect(non_visual,ind_rem);
    for i = 1:length(ind_rem)
        non_visual(non_visual == ind_rem(i)) = [];
    end
end
[r_tar,r_sub,~,~] = scatterDistAng(trials,0,0);
dist_err = r_sub - r_tar;
% [~,dist_err] = dist2target(trials,0,0); % choose one way of expressing distance from target

x = 100*stim_errX(non_visual)';
X = [ones(length(x),1) x];
y = dist_err(non_visual)';
b = regress(y,X);

figure;subplot(3,1,1);plot(100*stim_errX(non_visual),dist_err(non_visual),'.');grid on;
hold on;plot(x,x*b(2) + b(1),'r');title('stimulus X vs behavioral distance error');
xlabel('VR stimulus X error (cm)');ylabel('behavioral distance error (cm)');
suptitle(trials(1).prs.subject);

x = 100*stim_errY(non_visual)';
X = [ones(length(x),1) x];
y = dist_err(non_visual)';
b = regress(y,X);

subplot(3,1,2);plot(100*stim_errY(non_visual),dist_err(non_visual),'.');grid on;
hold on;plot(x,x*b(2) + b(1),'r');title('stimulus Y vs behavioral distance error');
xlabel('VR stimulus Y error (cm)');ylabel('behavioral distance error (cm)');
suptitle(trials(1).prs.subject);
x = 100*stim_err(non_visual)';
X = [ones(length(x),1) x];
y = dist_err(non_visual)';
b = regress(y,X);

subplot(3,1,3);plot(100*stim_err(non_visual),dist_err(non_visual),'.');grid on;
hold on;plot(x,x*b(2) + b(1),'r');title('stimulus  vs behavioral distance error');
xlabel('VR stimulus error (cm)');ylabel('behavioral distance error (cm)');
suptitle(trials(1).prs.subject);

%% sanity check
if sanity_check
figure;plot(trials(23).continuous.ts,trials(23).continuous.v);grid on;
hold on;plot(trials(23).mc.timestamp,100*trials(23).mc.Rat_X_GIAerror)
hold on;plot(trials(23).mc.timestamp,100*trials(23).mc.Rat_X_GIA)
hold on;plot(trials(23).mc.timestamp,100*cumsum(cumsum(trials(23).mc.Rat_X_GIAerror,1)*dt,1)*dt)
legend('Forward velocity (cm/s)','Forward GIA error (cm/s^2)','Forward GIA (cm/s^2)','Forward displacement error (cm)')
title(['Trial ' num2str(23) ', trajectory info'])
end