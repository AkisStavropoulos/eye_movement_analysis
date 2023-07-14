function errors_f_delta(trials,deltas)
%% Plot errors as a function of history
stats = [trials.stats]; 
d_err = [stats.d_err];
th_err = [stats.th_err];
%% plot distance error

x1 = deltas';
X1 = [x1 ones(length(x1),1)];
y = (d_err)';
[b]=regress(y,X1);
c1 = b(2);
b1 = b(1);
% figure;plot(deltas,d_err,'.');title(['distance error as a function of \tau difference, slope = ' num2str(round(b1))]);
% grid on;xlabel('\tau difference (s)');ylabel('distance error (cm)');
% hold on;plot(x1,x1*b1 + c1,'r');
% suptitle(trials(1).prs.subject);

% plot for absolute deltas

x2 = abs(deltas)';
X2 = [x2 ones(length(x2),1)];
y = (d_err)';
[b]=regress(y,X2);
c2 = b(2);
b2 = b(1);
% figure;plot(abs(deltas),d_err,'.');title(['distance error as a function of ABS \tau difference, slope = ' num2str(round(b2))]);
% grid on;xlabel('ABS \tau difference (s)');ylabel('distance error (cm)');
% hold on;plot(x2,x2*b2 + c2,'r');
% suptitle(trials(1).prs.subject);

% compare two plots
figure;
plot(abs(deltas),d_err,'.');hold on;plot(x2,x2*b2 + c2,'k');
plot(deltas,d_err,'.');plot(x1,x1*b1 + c1,'r');
title('distance error as a function of \tau difference');
grid on;xlabel('\tau difference (s)');ylabel('distance error (cm)');
legend('ABS \tau difference','ABS regression','SIGNED \tau difference','SIGNED regression');
suptitle(trials(1).prs.subject);
%% plot angular error
x1 = deltas';
X1 = [x1 ones(length(x1),1)];
y = (th_err)';
[b]=regress(y,X1);
c1 = b(2);
b1 = b(1);
% figure;plot(deltas,th_err,'.');title(['angular error as a function of \tau difference, slope = ' num2str(b1)]);
% grid on;xlabel('\tau difference (s)');ylabel('angular error (deg)');
% hold on;plot(x1,x1*b1 + c1,'r');
% ylim([-30 30]);
% suptitle(trials(1).prs.subject);

% plot for absolute deltas
x2 = abs(deltas)';
X2 = [x2 ones(length(x2),1)];
y = (th_err)';
[b]=regress(y,X2);
c2 = b(2);
b2 = b(1);
% figure;plot(abs(deltas),th_err,'.');title(['angular error as a function of ABS \tau difference, slope = ' num2str(b2)]);
% grid on;xlabel('\tau difference (s)');ylabel('angular error (deg)');
% hold on;plot(x2,x2*b2 + c2,'r');
% ylim([-30 30]);
% suptitle(trials(1).prs.subject);

% compare two plots
figure;
plot(abs(deltas),th_err,'.');hold on;plot(x2,x2*b2 + c2,'k');
plot(deltas,th_err,'.');plot(x1,x1*b1 + c1,'r');
ylim([-30 30]);
title('angular error as a function of \tau difference');
grid on;xlabel('\tau difference (s)');ylabel('angular error (deg)');
legend('ABS \tau difference','ABS regression','SIGNED \tau difference','SIGNED regression');
suptitle(trials(1).prs.subject);
