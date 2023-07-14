function [avg,sig] = all_jerk_hist(subject,N,logsc,plt)
%% Plot distribution of cumulative jerk of every trial and every subject
% logsc: 1 to plot in logscale, 0 for normal scale
% plt: 1 to plot 0 for no

count = 0;
% thresh = 60;
% N = 50;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count +1;
        jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
        acc(count) = subject(i).trials(j).continuous.cumacc(end);
    end
end
% calculate mean and variance
avg = nanmean(log(jerk));
sig = nanstd(log(jerk));

if logsc
xbins = logspace(1,5,N);
loglabl = 'LOG - ';
logtitl = ' - logspace';
jerklims = [0 10*10^4];
else
xbins = N;
loglabl = [];
logtitl = [];
jerklims = [0 2*10^4];
end

if plt
figure;
whitebg([1 1 1]);

subplot(2,1,1);hold on;
hist(jerk,xbins);xlabel([loglabl 'cumulative jerk^2' ]);xlim(jerklims);
title(['Distribution of cum jerk^2 of all trials' logtitl]);
if logsc
set(gca,'xScale','log');
vline(exp(avg),'r');vline([exp(avg - sig) ; exp(avg + sig)],'k');
end

subplot(2,1,2);hold on;
hist(acc,N);xlabel('cumulative Acceleration');title('Distribution of cum acceleration of all trials');hold on;
xlim([0 1500]);
if length(subject) == 1
    suptitle(subject.name) 
end
end
%% divide in coasters and brakers (nominees)
% count = 0;
% brakers = [];
% coasters = [];
% thresh = input('Where would you like to set your threshold, sir? ');
% for i = 1:length(subject)
%     for j = 1:length(subject(i).trials)
%         count = count +1;
%         if jerk(count) >= thresh
%             brakers = [brakers count];
%         else
%             coasters = [coasters count];
%         end
%     end
% end
% vline(thresh);