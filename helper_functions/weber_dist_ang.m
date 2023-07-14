function weber_dist_ang(subject,stype,tautype)
%% Test Weber law for target Distance and Angle
% input the field of the condition you want to test
trials = [subject.trials];
stimtype = [subject.stimtype];
tau = [subject.tau];
%% Distance
if ~isempty(stype) && ~isempty(tautype)
    indx = intersect(stimtype.(stype),tau.(tautype));
elseif isempty(stype) && ~isempty(tautype)
    indx = tau.(tautype);
elseif ~isempty(stype) && isempty(tautype)
    indx = stimtype.(stype);
else
    indx = 1:length(trials);
end

if ~isempty(stype)
    if strcmp(stype,'s1')
        stim = 'vestibular';
    elseif strcmp(stype,'s2')
        stim = 'visual';
    elseif strcmp(stype,'s3')
        stim = 'combined';
    end
else
    stim = 'all trials';
end

figure;
subplot(2,1,1);
d_err = [];
r_tar = [];
for i = 1:length(indx)
    d_err(i) = trials(indx(i)).stats.d_err;
    r_tar(i) = trials(indx(i)).prs.r_tar;
end
x = r_tar';
X = [x ones(length(x),1)];
y = d_err';
[b]=regress(y,X);
c = b(2);
b = b(1);
meany = mean(x*b + c);
sd = sqrt((d_err - meany).^2);
% plot(r_tar,d_err,'o');hold on;
% plot(x,x*b + c,'r','LineWidth',2);hold on;
% title('X');xlim([100 600]);xlabel('X [cm]');ylabel('distance error [cm]');hold off;
%
p = polyfit(r_tar,sd,2);
f = polyval(p,r_tar);
plot(r_tar,sd,'o');
hold on;plot(r_tar,f,'r.');grid on;
xlabel('target distance [cm]');ylabel('STD of distance error [cm]');title('SD of distance error over distance');
%

% sqrt
subplot(2,1,2);
d_err = [];
r_tar_sqrt = [];
for i = 1:length(indx)
    d_err(i) = trials(indx(i)).stats.d_err;
    r_tar_sqrt(i) = sqrt(trials(indx(i)).prs.r_tar);
end
p1 = polyfit(r_tar_sqrt,d_err,2);
f1 = polyval(p1,r_tar_sqrt);
sd = sqrt((d_err - f1).^2);
% plot(r_tar_sqrt,d_err,'o');hold on;
% plot(r_tar_sqrt,f1,'.r','LineWidth',2);
% title('sqrt(X)');xlim([0 25]);xlabel('sqrt(X)');ylabel('distance error');hold off;
% 
p = polyfit(r_tar_sqrt,sd,2);
f = polyval(p,r_tar_sqrt);
plot(r_tar_sqrt,sd,'o');
hold on;plot(r_tar_sqrt,f,'r.');grid on;
xlabel('SQRT of target distance [cm]');ylabel('STD of distance error [cm]');title('SD of distance error over SQRT of distance');
suptitle([trials(1).prs.subject ' - ' stim ' ' tautype]);

    fig(1) = gcf;
    fig(1).Units = 'centimeters';

%% Angle
% indx = intersect(stimtype.s2,tau.bin1);

figure;
subplot(2,1,1);
th_err = [];
th_tar = [];
for i = 1:length(indx)
    th_err(i) = trials(indx(i)).stats.th_err;
    th_tar(i) = trials(indx(i)).prs.th_tar;
end
x = th_tar';
X = [x ones(length(x),1)];
y = th_err';
[b]=regress(y,X);
c = b(2);
b = b(1);
meany = 0;%mean(x*b + c);
sd = abs(th_err);%sqrt((th_err - meany).^2);
% plot(th_tar,th_err,'o');hold on;
% plot(x,x*b + c,'r','LineWidth',2);hold on;
% title('X');xlim([-50 50]);xlabel('X [cm]');ylabel('distance error [cm]');hold off;
%
p = polyfit(abs(th_tar),sd,2);
f = polyval(p,abs(th_tar));
plot(abs(th_tar),sd,'o');
hold on;plot(abs(th_tar),f,'r.');grid on;
xlabel('target angle [deg]');ylabel('STD of angular error [deg]');title('SD of angular error over angle');
%

%sqrt
subplot(2,1,2);
th_err = [];
th_tar_sqrt = [];
for i = 1:length(indx)
    th_err(i) = trials(indx(i)).stats.th_err;
    th_tar_sqrt(i) = sqrt(abs(trials(indx(i)).prs.th_tar));
end
x = th_tar_sqrt';
X = [x ones(length(x),1)];
y = th_err';
[b]=regress(y,X);
c = b(2);
b = b(1);
meany = 0;%mean(x*b + c);
sd = abs(th_err);%sqrt((th_err - meany).^2);
% plot(r_tar_sqrt,d_err,'o');hold on;
% plot(x,x*b + c,'r','LineWidth',2);
% title('sqrt(X)');xlim([0 25]);xlabel('sqrt(X)');ylabel('distance error');hold off;
% 
p = polyfit(th_tar_sqrt,sd,2);
f = polyval(p,th_tar_sqrt);
plot(th_tar_sqrt,sd,'o');
hold on;plot(th_tar_sqrt,f,'r.');grid on;
xlabel('SQRT of target angle [deg]');ylabel('STD of angular error [deg]');title('SD of angular error over SQRT of angle');
suptitle([trials(1).prs.subject ' - ' stim ' ' tautype]);

fig(2) = gcf;
fig(2).Units = 'centimeters';
% place figures next to one another
count = 0;
for i = 1:length(fig)
    figdim = fig(i).Position(3:4);
   fig(i).Position = [0+(count*(figdim(1))) 15 figdim];  % 50: RL position , 18: max top position
    count = count + 1;    
end

    