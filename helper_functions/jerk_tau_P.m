function jerk_tau_P(subject,N)
%% heatmap of the joint distribution of jerk and tau: P(tau,J)
figure;
whitebg([1 1 1]);
jerk = [];
tau = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        tau(count) = subject(i).trials(j).prs.tau;
        jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
    end
end

[H,b] = hist3([jerk' tau'],'Nbins',[N N]);
imagesc(b{2}([1 end]),b{1}([1 end]),H);

c = colorbar('Ticks',round(linspace(0,max(H(:)),3),3),...
    'TickLabels',round(linspace(0,max(H(:)),3)./sum(H(:)),3));
c.Label.String = 'P(\tau,J^2)';
c.Label.FontSize = 12;
axis xy  tight;
xlim([min(tau) 5]);ylim([min(jerk) 10000]);
title('Joint Prob distribution of \tau and Jerk^2, all trials');xlabel('\tau [s]');ylabel('J [s^-^3]');
if length(subject) == 1
    suptitle(subject.name) 
end

%
%     x = [ceil(min(jerk)) ceil(max(jerk))];
%     y = [ceil(min(tau)) ceil(max(tau))];
%     plot(jerk,tau,'k.');grid on;
%     xlabel('Cumulative jerk');ylabel('\tau [s]');
