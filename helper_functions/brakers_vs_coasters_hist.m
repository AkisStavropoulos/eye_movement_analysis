function brakers_vs_coasters_hist(subject,brakers,coasters,M,thresh)
%% Plot joystick profiles of all trials, and the distribution of joystick inputs
% M: number of joystick input bins
edges = linspace(-1,1,M+1);
ind = 0;
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        ind = ind + 1;
        u{ind} = (subject(i).trials(j).mc.JS_X_Raw./subject(i).trials(j).prs.vmax)';
        if find(u{ind} > 1.1)
            count = count + 1;
            disp([subject(i).trials(j).prs.subject ' - trial ' num2str(j) ', max = ' num2str(max(abs(u{ind})))]);
            subplot(1,2,1);hold on;
            plot(subject(i).trials(j).mc.JS_X_Raw);grid on;hline([subject(i).trials(j).prs.vmax -subject(i).trials(j).prs.vmax]);
            
            subplot(1,2,2);hold on;
            plot(u{ind});
        end
    end
end

% brakers
figure;subplot(1,2,1);hold on;
for n = 1:length(brakers)
plot([u{brakers(n)}]);
end
title('Braker nominees');ylabel('JS forw/back input');xlabel('time [frames]');ylim([-2 2]);xlim([0 3000]);

tempa = [];
tempa = [u{brakers}];
tempa(tempa == 0) = [];

hold off;
subplot(1,2,2);h = histogram(tempa,edges,'Normalization','probability');
h.Orientation = 'horizontal';
hold on;hline(0);
title(['distribution of joystick input, Jerk threshold = ' num2str(thresh)]);
ylabel('joystick input');xlabel('Probability');xlim([0 1]);ylim([-2 2]);grid on;

if length(subject) == 1
    suptitle(subject.name) 
end

% coasters 
figure;subplot(1,2,1);hold on;
for n = 1:length(coasters)
plot([u{coasters(n)}]);
end
title('Coaster nominees');ylabel('JS forw/back input');xlabel('time [frames]');ylim([-2 2]);xlim([0 3000]);

tempb = [];
tempb = [u{coasters}];
tempb(tempb == 0) = [];

hold off;
subplot(1,2,2);h = histogram(tempb,edges,'Normalization','probability');
h.Orientation = 'horizontal';
hold on;hline(0);
title(['distribution of joystick input, Jerk threshold = ' num2str(thresh)]);
ylabel('joystick input');xlabel('Probability');xlim([0 1]);ylim([-2 2]);grid on;

if length(subject) == 1
    suptitle(subject.name) 
end
