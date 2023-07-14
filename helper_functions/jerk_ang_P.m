function jerk_ang_P(subject,N)
%% heatmap of the joint distribution of jerk and target angle: P(theta_tar,J)
figure;

jerk = [];
theta = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        theta(count) = abs(subject(i).trials(j).prs.th_tar);
        jerk(count) = subject(i).trials(j).continuous.jerk(end);
        if theta(count) > 40
            theta(count) = nan;
        end
    end
end


[H,b] = hist3([jerk' theta'],'Nbins',[N N]);
imagesc(b{2}([1 end]),b{1}([1 end]),H);

c = colorbar('Ticks',round(linspace(0,max(H(:)),3),3),...
    'TickLabels',round(linspace(0,max(H(:)),3)./sum(H(:)),3));
c.Label.String = 'P(\theta_t_a_r,J)';
c.Label.FontSize = 12;
ylim([0 50]);xlim([0 150*10^2]);
axis xy tight;
title('Joint Prob distribution of \theta_t_a_r_g_e_t and Jerk, all trials');ylabel('J^2 ');xlabel('\theta_t_a_r [deg]');
if length(subject) == 1
    suptitle(subject.name) 
end

%
%     x = [ceil(min(jerk)) ceil(max(jerk))];
%     y = [ceil(min(tau)) ceil(max(tau))];
%     plot(jerk,tau,'k.');grid on;
%     xlabel('Cumulative jerk');ylabel('\tau [s]');
