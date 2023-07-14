function [mean_err,std_err] = meanbin_signed_error(trials,binsize,condition,plt)
%% plot the mean error of the subject's trials
% output is sorted!!
% % create bins of 50cm intervals
% plt: plot yes 1 or not 0

%% BIN TAUS AND SEE WHAT HAPPENS!!
prs = [trials.prs];
r_tar = [prs.r_tar];
stats = [trials.stats];
d_err = [stats.d_err];

%%
if nargin <4
    plt = 0;
end

tar_min_d = 200; tar_max_d = 600;
% binsize = 100;
tar_bins = tar_min_d:binsize:tar_max_d;

if isempty(condition) % for all trials
    % sort errors based on the distance of the intended target
    [r_tar_sort,indx1] = sort(r_tar);
    err_sort = d_err(indx1);
    % find mean and std of bins
    mean_err = [];
    std_err = [];
    for i = 1:length(tar_bins)-1
        indx2 = find(r_tar_sort >= tar_bins(i) & r_tar_sort <= tar_bins(i+1));
        bin_mean = mean(err_sort(indx2));
        bin_std = std(err_sort(indx2))/sqrt(numel(indx2));
        mean_err = [mean_err bin_mean];
        std_err = [std_err bin_std];
    end
    if plt
        figure;errorbar(tar_bins(2:end),mean_err,std_err);xlim([0 800]);title('mean signed Distance error');
        xlabel('initial target distance (cm)');ylabel('mean error (cm)');
    end
else
    %% name the conditions for the legend
    condnames = fieldnames(condition);
    if any(strcmp(condnames,'s1bin1')) % check for intersection first
        for i = 1:length(condnames)-1
            if ~strcmp(condnames{i}(1:3),'lim')
                condS{i} = condnames{i}(1:2);
                condJS{i} = condnames{i}(3:end);
            end
        end
        condS =  condS(~cellfun('isempty',condS));
        condS = unique(condS);
        condJS =  condJS(~cellfun('isempty',condJS));
        condJS = unique(condJS);
        if any(strcmp(condS,'s1')) && any(strcmp(condS,'s2'))
            condSnames = [{'vestibular'} {'visual'} {'combined'} ];
        end
        if any(strcmp(condJS,'bin1')) && any(strcmp(condJS,'bin2'))
            for i = 1:length(condJS)
                if strcmp(condJS{i}(1),'b')
                    condJSnames{i} = ['\tau ' condJS{i}];
                    condJSbins{i} = condJS{i};
                else
%                     condJSlims{i} = ['\tau ' condJS{i}];
                end
            end
        end
        condJS = unique(condJSbins);
        
    else
        if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2'))
            for i = 1:length(condnames)-1
                condS{i} = condnames{i};
            end
            condSnames = [{'vestibular'} {'visual'} {'combined'} ];
            condJS = {[]};
            condJSnames = {[]};
        elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
            for i = 1:length(condnames)-1
                if strcmp(condnames{i}(1:3),'bin')
                    condJS{i} = condnames{i};
                    condJSnames{i} = ['\tau ' condJS{i}];
                end
            end
            condS = {[]};
            condSnames = {[]};
        end
    end
    if plt
        h(1) = figure;
        if length(condnames) > 3
            h(2) =  figure;
        end
    end

    %%
    
    % create legend
    leg_inp_temp = {[condSnames{i} ' ' condJSnames{j}]};
    leg_inp = [leg_inp ; repmat(leg_inp_temp,2,1)];
    leg_input{i} = leg_inp(:);

    %     cmap = jet(length(condnames));
    for j = 1:length(condnames)
        indx1 = [];
        r_tar_sort = [];
        err_sort = [];
        % sort errors based on the distance of the intended target
        [r_tar_sort,indx1] = sort(r_tar(condnames{j}));
        err_sort = d_err(condnames{j}(indx1));
        % find mean and std of bins
        mean_err1 = [];
        std_err1 = [];
        for i = 1:length(tar_bins)-1
            indx2 = find(r_tar_sort >= tar_bins(i) & r_tar_sort <= tar_bins(i+1));
            bin_mean = mean(err_sort(indx2));
            bin_std = std(err_sort(indx2))/sqrt(numel(indx2));
            mean_err1 = [mean_err1 bin_mean];
            std_err1 = [std_err1 bin_std];
        end
        % save data
        mean_err.(condnames{j}) = mean_err1;
        std_err.(condnames{j}) = std_err1;
    end
    
    
    % plot data
    if plt
        sbpltcols = (length(condnames));
        colr = lines(sbpltcols);
        
        if any(strcmp(condnames{1},'s1')) || any(strcmp(condnames{1},'bin1')) % stimtype OR tau bins
            for j = 1:length(condnames)
                errorbar(tar_bins(2:end),mean_err.(condnames{j}),std_err.(condnames{j}),'Color',colr(j,:),'Marker','.');xlim([0 800]);
                suptitle(trials(1).prs.subject);title(['mean signed Distance error, target bins: ' num2str(binsize) ' cm']);
                xlabel('initial target distance (cm)');ylabel('mean error (cm)');
                legend(condition,'location','northwest');hold on;
            end
        else
            % stimtype AND jscoef
            Scond = unique(conditionS);
            JScond = unique(conditionJS);
            % plot taus across modalities
            figure(h(1));
            for i = 1:length(Scond)
                indx = [];
                for j = 1:length(condnames)
                    if strcmp(condnames{j}(1:2),Scond{i})
                        indx = [indx j];
                    end
                end
                subplot(1,length(Scond),i);hold on;
                for n = 1:length(indx)
                    errorbar(tar_bins(2:end),mean_err.(condnames{indx(n)}),std_err.(condnames{indx(n)}),'Color',colr(n,:),'Marker','.');
                    xlim([0 800]);ylim([0 350]);
                    suptitle([trials(1).prs.subject ' - target bins: ' num2str(binsize) ' cm']);title(['mean signed Distance error for ' conditionStim{i}]);
                    xlabel('initial target distance (cm)');ylabel('mean error (cm)');hold on;
                    legend(conditionJS{:},'location','northwest');
                    
                end
            end
            % plot modalities across taus
            figure(h(2));
            for i = 1:length(JScond)
                indx = [];
                for j = 1:length(condnames)
                    if strcmp(condnames{j}(3:end),JScond{i}(end-3:end))
                        indx = [indx j];
                    end
                end
                subplot(1,length(JScond),i);hold on;
                for n = 1:length(indx)
                    errorbar(tar_bins(2:end),mean_err.(condnames{indx(n)}),std_err.(condnames{indx(n)}),'Color',colr(n,:),'Marker','.');
                    xlim([0 800]);ylim([0 350]);
                    suptitle([trials(1).prs.subject ' - target bins: ' num2str(binsize) ' cm']);title(['mean signed error for ' conditionJS{i}]);
                    xlabel('initial target distance (cm)');ylabel('mean error (cm)');hold on;
                    legend(conditionStim{:},'location','northwest');
                end
            end
        end
    end
end