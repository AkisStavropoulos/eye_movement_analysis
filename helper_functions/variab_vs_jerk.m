function [rho,pval] = variab_vs_jerk(subject,local,N)
%% Plot trial by trial variance vs jerk
% z-score both std and jerk ??
% choose local or global mean
% local = 1;
% N = 50; % trials window for moving mean
for i = 1:length(subject)
    % compute mean bias 
    dx = [];
    dy = [];
    dj = [];
    jerk = [];
    de_glob = [];
    de_loc = [];
    loc_avg_x = [];
    loc_avg_y = [];
    for j = 1:length(subject(i).trials)
        % compute euclideian error
        dx(j) = subject(i).trials(j).continuous.xmp(end) - subject(i).trials(j).prs.fireflyposx;
        dy(j) = subject(i).trials(j).continuous.ymp(end) - subject(i).trials(j).prs.fireflyposy;
        
        % extract total jerk for each trial
        jerk(j) = subject(i).trials(j).continuous.jerk(end);

    end
    % calculate local mean
    loc_avg_x{i} = movmean(dx,N); loc_avg_x{i}([1:N/2,end-N/2:end]) = nan;
    loc_avg_y{i} = movmean(dy,N); loc_avg_y{i}([1:N/2,end-N/2:end]) = nan;

    % subtract local mean
    loc_nb_dx{i} = dx - loc_avg_x{i};
    loc_nb_dy{i} = dy - loc_avg_y{i};
    
    % calculate global mean
    glob_avg_dx(i) = mean(dx);
    glob_avg_dy(i) = mean(dy);
    
    % subtract global mean
    glob_nb_dx{i} = dx - glob_avg_dx(i);
    glob_nb_dy{i} = dy - glob_avg_dy(i);
    
    % compute moving average of total jerk
    dj = movmean(jerk,N);
    
    % calculate variability for local ang global mean
    for j = 1:length(subject(i).trials)
        if j > N/2 && j < length(subject(i).trials) - N/2
            de_glob(j) = sqrt(mean(glob_nb_dx{i}(j-N/2:j+N/2).^2 + glob_nb_dy{i}(j-N/2:j+N/2).^2));
            de_loc(j) = sqrt(mean(loc_nb_dx{i}(j-N/2:j+N/2).^2 + loc_nb_dy{i}(j-N/2:j+N/2).^2));
        else
            de_glob(j) = nan;
            de_loc(j) = nan;
        end
    end
        
    if local
        % compute correlation and plot for local
    [rho(i),pval(i)] = nancorr(dj',de_loc');
    figure;hold on;
    colr = parula(length(dj));
    for n = 1:length(dj)
        plot(dj(n),de_loc(n),'Marker','o','MarkerFaceColor',colr(n,:));
    end
    xlim([0 8000]);ylim([0 300]);
    xlabel('jerk');ylabel('variability [cm]');title(['LOCAL variability vs jerk - ' subject(i).name]);
    c = colorbar('Ticks',linspace(0,1,5),'TickLabels',[1:ceil(length(dj)/5):length(dj)]);
    c.Label.String = 'Trial No.';
    suptitle(['\rho = ' num2str(rho(i)) ', pval = ' num2str(pval(i))]);
    %     figure;plot(dj);hold on;plot(de);title(['variability and jerk - ' subject(i).name]);
    %     jerk_vs_exp(subject(i).trials,subject(i).stimtype,'lin');
    else
    % compute correlation and plot for global
    [rho(i),pval(i)] = nancorr(dj',de_glob');
    figure;hold on;
    colr = parula(length(dj));
    for n = 1:length(dj)
        plot(dj(n),de_glob(n),'Marker','o','MarkerFaceColor',colr(n,:));
    end
    xlim([0 8000]);ylim([0 300]);
    xlabel('jerk');ylabel('variability [cm]');title(['GLOBAL variability vs jerk - ' subject(i).name]);
    c = colorbar('Ticks',linspace(0,1,5),'TickLabels',[1:ceil(length(dj)/5):length(dj)]);
    c.Label.String = 'Trial No.';
    suptitle(['\rho = ' num2str(rho(i)) ', pval = ' num2str(pval(i))]);
    %     figure;plot(dj);hold on;plot(de);title(['variability and jerk - ' subject(i).name]);
    %     jerk_vs_exp(subject(i).trials,subject(i).stimtype,'lin');
    end
end

