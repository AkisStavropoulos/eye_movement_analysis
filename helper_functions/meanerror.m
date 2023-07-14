function [mean_err,std_err] = meanerror(r_tar,err)
%% plot the mean error of the subject's trials

% % create bins of 50cm intervals
tar_min_d = 100; tar_max_d = 600;
tar_bins = tar_min_d:50:tar_max_d;

if ~isstruct(r_tar)
    % sort errors based on the distance of the intended target
    [r_tar_sort,indx1] = sort(r_tar);
    err_sort = err(indx1);
    % find mean and std of bins
    mean_err = [];
    std_err = [];
    for i = 1:length(tar_bins)-1
        indx2 = find(r_tar_sort >= tar_bins(i) & r_tar_sort <= tar_bins(i+1));
        bin_mean = mean(err_sort(indx2));
        bin_std = std(err_sort(indx2));
        mean_err = [mean_err bin_mean];
        std_err = [std_err bin_std];
    end
    
    figure;errorbar(tar_bins(2:end),mean_err,sqrt(std_err));xlim([0 800]);title('mean error');
    xlabel('initial target distance (cm)');ylabel('mean error (cm)');
else
    figure;
    condnames = fieldnames(err);
    for j = 1:length(condnames)
        indx1 = [];
        % sort errors based on the distance of the intended target
        [r_tar_sort(j),indx1] = sort(r_tar.(condnames{j}));
        err_sort(j) = err.(condnames{j})(indx1);
        % find mean and std of bins
        mean_err1 = [];
        std_err1 = [];
        for i = 1:length(tar_bins)-1
            indx2 = find(r_tar_sort >= tar_bins(i) & r_tar_sort <= tar_bins(i+1));
            bin_mean = mean(err_sort(indx2));
            bin_std = std(err_sort(indx2));
            mean_err1 = [mean_err1 bin_mean];
            std_err1 = [std_err1 bin_std];
        end
        mean_err.(condnames{j}) = mean_err1;
        std_err.(condnames{j}) = std_err1;
        errorbar(tar_bins(2:end),mean_err1,sqrt(std_err1));xlim([0 800]);title('mean error');
        xlabel('initial target distance (cm)');ylabel('mean error (cm)');hold on;
    end
        hold off;
end