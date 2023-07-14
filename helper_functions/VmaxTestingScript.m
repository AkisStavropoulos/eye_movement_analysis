%% vmax testing script
% see what the maximum velocity of all trials is
for i = 1:length(trials)
    vmax(i) = max(trials(i).continuous.v);
    coef(i) = trials(i).prs.js_coef;
end

[coef,indx1] = sort(coef);
vmax = vmax(indx1);
% 
% figure;plot(coef,vmax);xlabel('joystick coefficient');ylabel('vmax')
% 
% figure; subplot(2,1,1);plot(vmax);ylabel('vmax');subplot(2,1,2);plot(coef);xlabel('number of trials');ylabel('JS coefficient')

%% plot all velocity traces with different colors (for each coefficient?)
cmap = jet(length(trials));
figure;
vtrial = [];tstrial = [];
for i = 1:length(trials)
    vtrial{i} = trials(i).continuous.v;
    tstrial{i} = trials(i).continuous.ts;
    plot(tstrial{i},vtrial{i},'Color',cmap(i,:));hold on;
    title('velocity of every trial');xlabel('time (s)');ylabel('velocity (cm/s)')
end

%% compute acceleration for each trial
figure;
for i = 1:length(trials)
    v = []; acc = [];
    v = trials(i).continuous.v;
    acc = [0 ; diff(v)/prs.dt];
    trials(i).continuous.acc = acc;
    amax(i) = max(acc);
    
    if amax(i) >= 10 % plot only trials that you moved
        plot(trials(i).continuous.ts,trials(i).continuous.acc,'Color',cmap(i,:));hold on;xlabel('time (s)');ylabel('linear acceleration (cm/s^2)');
        title('linear acceleration of every trial');
    end
    vmaxinp(i) = trials(i).prs.vmax;
end
% plot inputed vmax, reached vmax, amax for non-zero trials
indx = find(amax >= 10);
figure;subplot(3,2,1);plot(amax(indx),'*-');title('amax across trials');ylabel('amax (cm/s^2)');xlabel('trials')
subplot(3,2,3);plot(vmax(indx),'*-');title('vmax across trials');ylabel('vmax (cm/s)');xlabel('trials')
subplot(3,2,5);plot(vmaxinp(indx),'*-k');title('inputed (expected) vmax across trials');ylabel('inputed vmax (cm/s)');xlabel('trials')

%% compute expected acceleration for given vmax input
% for i = 1:length(trials)
% vmax(i) = [trials(i).prs.vmax];
% end
vmaxinp = vmaxinp(indx);
vtrial = vtrial(indx);
tstrial = tstrial(indx);

for i = 1:length(vmaxinp)
    [~,amax_exp(i)] = vmax2acc(.99,vmaxinp(i),10,1);
    hold on;plot(tstrial{i},vtrial{i},'r.');ylim([0 200]);xlim([0 35])
end
subplot(3,2,2);plot(amax_exp(indx),'r*-');title('expected amax given vmax input');xlabel('trials');ylabel('amax (cm/s^2)');

%% expected relationship between vmax and amax
figure;plot(amax_exp(indx),vmaxinp(indx),'r*-');title('amax over vmax');xlabel('vmax (cm/s)');ylabel('amax (cm/s^2)');