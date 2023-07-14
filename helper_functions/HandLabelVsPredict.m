%% Compare Hand-labeling vs Predictive method    
%% Fix predicted Coast vector    
    pred_zero_ind = find(pred_coast{i} == 0);
    flag_nan_ind = find(isnan(trl_flag{i}));
        
    pred_coast{i}(pred_zero_ind) = nan;
    
    %% Fix stopping Phase vector
    stop_zero_ind = find(stop_phase{i} == 0);
    flag_nan_ind = find(isnan(trl_flag{i}));
    
    stop_phase{i}(stop_zero_ind) = nan;
    
    %% Plot histogram of StoppingPhase/PredictedCoast along with flags
    colr_flag = brewermap(length(subject),'Dark2');
    figure;
    N = 3;
    [h{i},c{i}] = hist(stop_phase{i}./pred_coast{i},N);
    
    plot(c{i},h{i}/sum(h{i}),'Color',colr(i,:),'LineWidth',2);grid on;
    title('Distribution of ratio of StoppingPhase / PredictedCoasting');
    xlim([0 2]);xlabel('StoppingPhase / PredictedCoast');
    ylabel('Probability')
    hold on;
    N = 3;
    [f{i},b{i}] = hist(trl_flag{i}/2,N);grid on;
    plot(b{i},f{i}/sum(f{i}),'m','LineWidth',2);grid on;
    vline(1,'r-');
    legend({'Prediction','Hand-labeling'})
    suptitle(subject(i).name);
    %%
    figure;
    N = 3;
    histogram(stop_phase{i}./pred_coast{i},N,'Normalization','pdf');grid on;
    title('Distribution of ratio of StoppingPhase / PredictedCoasting');
    xlim([0 2]);xlabel('StoppingPhase / PredictedCoast');
    
    hold on;
    N = 3;
    histogram(trl_flag{i}/2,N,'Normalization','pdf');grid on;
    vline(1,'r-','Coasting');vline(0.2,'r-','Braking');
    legend({'Prediction','Hand-labeling'})
    suptitle(subject(i).name);
    %% One-to-One comparison
    notnanindx = find((~isnan(trl_flag{i})));
    x = trl_flag{i} + 1;
    y = stop_phase{i}./pred_coast{i};
    
    flagbin = unique(x(notnanindx));
    ybin = cell(length(flagbin),1);
    binmean = [];
    binstd = [];
    for n = 1:length(flagbin)
        indx = [];
        indx = find(x == flagbin(n));
        ybin{n} = y(indx);
        binmean(n) = nanmedian(ybin{n});
        binstd(n) = nanstd(ybin{n});
        
    end
    figure;
    subplot(1,2,1);
    plot(x,y,'o','MarkerSize',6);
    hold on;
    errorbar(flagbin,binmean,binstd,binstd,'r.-','MarkerSize',16);
    grid on;xlim([0 max(x)+min(x)]);ylim([0 1.5]);
    
    N = 3;
    for n = 1:length(flagbin)
        bh = [];
        bc = [];
        [bh,bc] = hist(ybin{n},N,'Normalization','pdf');grid on;
        plot(bh/sum(bh),bc,'LineWidth',4);
        
    end
    
    grid on;xlim([0 max(x)+min(x)]);ylim([0 1.5]);
    xlabel('Hand-labeling');
    xticks(sort(unique(x(notnanindx))));
    xticklabels({'Braking','Moderate Braking','Coasting'});
    ylabel('StoppingPhase / PredictedCoast');
    title('The Battle of Cri-Tyrion (vs Hand-Labeling)');
    
    subplot(1,2,2);
    N = 3;
    histogram(y,N,'Normalization','pdf','Orientation','horizontal');grid on;
    title('Distribution of ratio of StoppingPhase / PredictedCoasting');
    ylim([0 1.5]);ylabel('StoppingPhase / PredictedCoast');
    xlabel('P(prediction) and P(hand-label)');
    hold on;
    N = 3;
    histogram(trl_flag{i}/2,N,'Normalization','pdf','Orientation','horizontal');grid on;
    suptitle(subject(i).name);
    hline(1,'r-','Coasting');hline(0.2,'r-','Braking');
    legend('Prediction','Hand-labeling');
    
    suptitle(subject(i).name);




