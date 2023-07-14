function  [rho,pval] = jerk_dist_all(subject,regr)
%% heatmap of the The relationship between jerk and target distance and angle
% regr: choose what type of regression you want to plot

figure;
whitebg([1 1 1]);

jerk = [];
dist = [];
count = 0;
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        count = count + 1;
        dist(count) = subject(i).trials(j).prs.r_tar;
        jerk(count) = subject(i).trials(j).continuous.cumjerk(end);
    end
end
[rho,pval] = nancorr(dist',jerk');

% [H,b] = hist3([dist' jerk'],'Nbins',[N N]);
% imagesc(b{2}([1 end]),b{1}([1 end]),H);
plot(dist,jerk,'.','MarkerSize',5);

% c = colorbar('Ticks',round(linspace(0,max(H(:)),3),3),...
%     'TickLabels',round(linspace(0,max(H(:)),3)./sum(H(:)),3));
% c.Label.String = 'P(\tau,J)';
% c.Label.FontSize = 12;
xlim([200 600]);ylim([0 150*10^2]);
% axis xy  tight;
title(['Jerk^2 vs d_t_a_r_g_e_t, all trials - \rho = ' num2str(rho) ', pval = ' num2str(pval)]);
xlabel('target distance [cm]');ylabel('cumulative Jerk^2');

if strcmp(regr,'poly')
    p = polyfit(dist,jerk,2);
    f = polyval(p,dist);
    [~,I] = sort(dist);
    hold on;plot(dist(I),f(I),'LineWidth',2.5);
elseif strcmp(regr,'lin')
    [x,y,b,c] = find_regress(dist,jerk);
    hold on;plot(x,b*x + c,'LineWidth',2.5);
elseif strcmp(regr,'movavg')
    N = 50;
    [dist,I] = sort(dist);
    hold on;
    M = movmean(jerk(I),N,'omitnan');
    nanind = find(isnan(jerk(I)));
    X = dist;
    X(nanind) = nan;
    plot(X,M,'LineWidth',2.5);
end


if length(subject) == 1
    suptitle(subject.name) 
end
