function plot_mean_std_history_error(trials,stimtype)

%% multivariate linear regression of mean and std of tau history
stnames = fieldnames(stimtype);
condition = {'vestibular','visual','combined'};

for n = 1:length(stnames)-1
    stats = [];
    d_err = [];
    hist_std = [];
    mean_tau = [];
    
    for j = 1:length(stimtype.(stnames{n}))
        mean_tau(j) = mean(trials(stimtype.(stnames{n})(j)).stats.tau_history);
        hist_std(j) = std(trials(stimtype.(stnames{n})(j)).stats.tau_history);
    end
    stats = [trials(stimtype.(stnames{n})).stats];
    d_err = [stats.d_err];
    x1 = hist_std';
    x2 = mean_tau';
    y = d_err';
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

