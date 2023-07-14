%% multivariate linear regression of mean and std of tau history
condition = {'vestibular','visual','combined'};

for n = 1:length(stnames)-1
    stats = [];
    th_err = [];
    hist_std = [];
    mean_tau = [];
    
    for j = 1:length(stimtype.(stnames{n}))
        mean_tau(j) = mean(trials(stimtype.(stnames{n})(j)).stats.tau_history);
        hist_std(j) = std(trials(stimtype.(stnames{n})(j)).stats.tau_history);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err = [stats.d_err];
    x1 = hist_std';
    x2 = mean_tau';
    y = th_err';
    X = [ones(size(x1)) x1 x2 x1.*x2];
    b = regress(y,X);
    figure;
    scatter3(x1,x2,y,'filled');hold on;
    x1fit = min(x1):.1:max(x1);
    x2fit = min(x2):.1:max(x2);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
    mesh(X1FIT,X2FIT,YFIT);
    xlabel('SD of \tau');xlim([0 2.5]);
    ylabel('Mean \tau');ylim([0 5]);
    zlabel('Distance error (cm)');zlim([-250 250]);
    title([condition{n} ' - Multi regression: \epsilon_d (\mu(\tau),\sigma(\tau))']);
    suptitle(trials(1).prs.subject);
    
end

%% see error as a function of std of recent tau history
% std_tau = sqrt(sum(tauhistory - mean(tau_history))^2)
figure;
colr = colormap(jet);
colr = colr(round(linspace(1,length(colr),length(stnames))),:);
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);
for n = 1:length(stnames)-1
    stats = [];
    th_err = [];
    hist_std = [];
    for j = 1:length(stimtype.(stnames{n}))
        
        hist_std(j) = std(trials(stimtype.(stnames{n})(j)).stats.tau_history);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err = [stats.d_err];
    [x,y,b,c] = find_regress(hist_std,th_err);
    plot(hist_std,th_err,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('\sigma(\tau)');ylabel('error (cm)');
title('Distance error as a function of SD of recent \tau history');legend(leg_input);

%% see error as a function of mean of recent tau history
figure;
colr = colormap(jet);
colr = colr(round(linspace(1,length(colr),length(stnames))),:);
condition = {'vestibular','visual','combined'};
leg_input = [condition; condition];
leg_input = leg_input(:);

for n = 1:length(stnames)-1
    stats = [];
    th_err = [];
    mean_tau = [];
    for j = 1:length(stimtype.(stnames{n}))
        mean_tau(j) = mean(trials(stimtype.(stnames{n})(j)).stats.tau_history);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    th_err = [stats.d_err];
    [x,y,b,c] = find_regress(mean_tau,th_err);
    plot(mean_tau,th_err,'.','Color',colr(n,:));hold on;plot(x,b*x + c,'Color',colr(n,:));
end
grid on;xlabel('mean(\tau)');ylabel('error (cm)');
title('Distance error as a function of MEAN of recent \tau history');legend(leg_input);

%% error as a function of std of recent delta history
for i = 1:length(trials)
    hist_mean(i) = mean(trials(i).stats.delta_history);
end
stats = [trials.stats];
th_err = [stats.d_err];

[x,y,b,c] = find_regress(hist_mean,th_err);


figure;plot(hist_mean,th_err,'.');hold on;plot(x,b*x + c,'r');grid on;

%%
for n = 1:bins
indx.(['bin' num2str(n)]) = [];
end

for i = 2:length(trials)
    hist_bins = trials(i).stats.hist_bins;
    for n = 1:bins
        if hist_bins(1) == n
            indx.(['bin' num2str(n)]) = [indx.(['bin' num2str(n)]) i];
        end
    end
end

stats = [trials.stats];
th_err = [stats.d_err];
figure;hold on;
for n = 1:bins
    bin_err = th_err(indx.(['bin' num2str(n)]));
    mean_bin_err(n) = mean(bin_err);
    std_bin_err(n) = std(bin_err);
    errorbar(n,mean_bin_err(n),std_bin_err(n));
end
plot(1:n,mean_bin_err,'o-');
xlim([0 6]);title('distance error as a function of increasing delta \tau bins, -3 to 3');
ylabel('distance error (cm)');xlabel('delta \tau bins');

% now plot as a function of the previous delta