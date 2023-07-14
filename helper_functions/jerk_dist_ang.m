function jerk_dist_ang(subject)
%% heatmap of the The relationship between jerk and target distance and angle
figure;

jerk = [];
dist = [];
ang = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        dist(count) = subject(i).trials(j).prs.r_tar;
        ang(count) = abs(subject(i).trials(j).prs.th_tar);
        jerk(count) = subject(i).trials(j).continuous.jerk(end);
    end
end

% [H,b] = hist3([dist' jerk'],'Nbins',[N N]);
% imagesc(b{2}([1 end]),b{1}([1 end]),H);
plot3(dist,ang,jerk,'.','MarkerSize',5);

% c = colorbar('Ticks',round(linspace(0,max(H(:)),3),3),...
%     'TickLabels',round(linspace(0,max(H(:)),3)./sum(H(:)),3));
% c.Label.String = 'P(\tau,J)';
% c.Label.FontSize = 12;
ylim([0 40]);xlim([200 600]);
% axis xy  tight;
title('Jerk vs target distance and angle, all trials');
xlabel('target distance [cm]');ylabel('target angle [deg]');zlabel('cumulative Jerk [s^-^3]');
if length(subject) == 1
    suptitle(subject.name) 
end
