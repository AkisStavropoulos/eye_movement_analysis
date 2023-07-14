function jerk_tau_ang(subject)
%% heatmap of the The relationship between jerk and tau and target angle
figure;

jerk = [];
tau = [];
ang = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        tau(count) = subject(i).trials(j).prs.tau;
        ang(count) = abs(subject(i).trials(j).prs.th_tar);
        jerk(count) = subject(i).trials(j).continuous.jerk(end);
    end
end

% [H,b] = hist3([dist' jerk'],'Nbins',[N N]);
% imagesc(b{2}([1 end]),b{1}([1 end]),H);
h = plot3(tau,ang,jerk,'.','MarkerSize',5);
grid on;

% c = colorbar('Ticks',round(linspace(0,max(H(:)),3),3),...
%     'TickLabels',round(linspace(0,max(H(:)),3)./sum(H(:)),3));
% c.Label.String = 'P(\tau,J)';
% c.Label.FontSize = 12;
ylim([0 40]);xlim([0 8]);
% axis xy  tight;
title('Jerk vs tau and target angle, all trials');
xlabel('\tau [s]');ylabel('target angle [deg]');zlabel('cumulative Jerk [s^-^3]');
if length(subject) == 1
    suptitle(subject.name) 
end
