%% wmax testing script

% not so succesful

% see what the maximum velocity of all trials is
for i = 1:length(trials)
    wmax(i) = max(trials(i).continuous.w);
    wmin(i) = min(trials(i).continuous.w);
    coef(i) = trials(i).prs.js_coef;
end

[coef,indx1] = sort(coef);
wmax = wmax(indx1);
wmin = wmin(indx1);
% 
% figure;plot(coef,vmax);xlabel('joystick coefficient');ylabel('vmax')
% 
% figure; subplot(2,1,1);plot(vmax);ylabel('vmax');subplot(2,1,2);plot(coef);xlabel('number of trials');ylabel('JS coefficient')

%% plot all velocity traces with different colors (for each coefficient?)
cmap = jet(length(trials));
figure;
wtrial = [];tstrial = [];
for i = 1:length(trials)
    xtrial{i} = trials(i).continuous.xmp;
    ytrial{i} = trials(i).continuous.ymp;
    wtrial{i} = trials(i).continuous.w;
    tstrial{i} = trials(i).continuous.ts;
    plot(tstrial{i},wtrial{i},'Color',cmap(i,:));hold on;
    title('angular velocity of every trial');xlabel('time (s)');ylabel('angular velocity (deg/s)')
end

%% compute acceleration for each trial
figure;
for i = 1:length(trials)
    w = []; wacc = [];
    w = trials(i).continuous.w;
    wacc = [0 ; diff(w)/prs.dt];
    trials(i).continuous.wacc = wacc;
    wamax(i) = max(wacc);
    wamin(i) = min(wacc);
    
%     if wamax(i) >= 10 % plot only trials that you moved
        plot(trials(i).continuous.ts,trials(i).continuous.wacc,'Color',cmap(i,:));hold on;xlabel('time (s)');ylabel('angular acceleration (deg/s^2)');
        title('angular acceleration of every trial');
%     end
    wmaxinp(i) = trials(i).prs.wmax;
end
% plot inputed wmax, reached wmax, amax for non-zero trials
% indx = find(wamax >= 0);
wamax = wamax(indx1);
wamin = wamin(indx1);
wmax = wmax(indx1);
wmin = wmin(indx1);
wmaxinp = wmaxinp(indx1);
% 
% for i = 1:length(ytrial)
%     ydiff(i) = abs(ytrial{i}(1) - ytrial{i}(end));
% end
indx = 1:2:length(trials);

figure;subplot(3,2,1);plot(wamax,'*-');title('amax across trials');ylabel('amax (deg/s^2)');xlabel('trials');
hold on; plot(wamin,'r*-');legend('maximum angular acc','minimum angular acc');hold off;
subplot(3,2,3);plot(wmax,'*-');title('wmax across trials');ylabel('wmax (deg/s)');xlabel('trials');
hold on;plot(wmin,'r*-');legend('maximum angular v','minimum angular v');hold off;
subplot(3,2,5);plot(wmaxinp,'*-k');title('inputed (expected) wmax across trials');ylabel('inputed wmax (deg/s)');xlabel('trials')

%% compute expected acceleration for given vmax input
% for i = 1:length(trials)
% vmax(i) = [trials(i).prs.vmax];
% end
wmaxinp = wmaxinp(indx);
wtrial = wtrial(indx);
tstrial = tstrial(indx);

for i = 1:length(wmaxinp)
    [~,amax_exp(i)] = vmax2acc(.99,wmaxinp(i),10,1);
    hold on;plot(tstrial{i},wtrial{i},'r.');ylim([0 200]);xlim([0 35])
end
subplot(3,2,2);plot(amax_exp(indx),'r*-');title('expected amax given vmax input');xlabel('trials');ylabel('amax (cm/s^2)');

%% expected relationship between vmax and amax
figure;plot(amax_exp(indx),wmaxinp(indx),'r*-');title('amax over vmax');xlabel('vmax (cm/s)');ylabel('amax (cm/s^2)');