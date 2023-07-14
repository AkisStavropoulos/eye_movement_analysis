function [mean_err,std_err] = meanbinerror_th(trials,condition,binsize,plt)
%% plot the mean error of the subject's trials
% output is sorted!!
% % create bins of 50cm intervals
% plt: plot yes 1 or not 0

%% BIN TAUS AND SEE WHAT HAPPENS!!
prs = [trials.prs];   
th_tar = [prs.th_tar];
stats = [trials.stats];
th_err = [stats.th_err];


%%
if nargin <4
    plt = 0;
end

tar_min_ang = 0; tar_max_ang = 40;
% binsize = 10;
tar_bins = tar_min_ang:binsize:tar_max_ang;


if isempty(condition) % for all trials
    % sort errors based on the distance of the intended target
    [r_tar_sort,indx1] = sort(th_tar);
    err_sort = th_err(indx1);
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
        figure;errorbar(tar_bins(2:end),mean_err,std_err);xlim([0 800]);ylim([-300 300]);
        title('mean signed Angular error');grid on;
        xlabel('initial target angle [deg]');ylabel('mean error [deg]');
        suptitle([trials(1).prs.subject ' - all trials']);
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
            condS =  condS(~cellfun('isempty',condS));
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
            condJS =  condJS(~cellfun('isempty',condJS));
            condJSnames =  condJSnames(~cellfun('isempty',condJSnames));
            condS = {[]};
            condSnames = {[]};
        end
    end
    
    count = 1;
    for i = 1:length(condS)
        for j = 1:length(condJS)
            cond{count} = [condS{i} condJS{j}];
            leg_input{count} = [condSnames{i} ' and ' condJSnames{j}];
            count = count+1;
        end
    end
    
    if plt
        h(1) = figure;
        if length(cond) > 3
            h(2) =  figure;
        end
    end
%%
    for j = 1:length(cond)
        indx1 = [];
        r_tar_sort = [];
        err_sort = [];
        % sort errors based on the distance of the intended target
        [r_tar_sort,indx1] = sort(th_tar(condition.(cond{j})));
        err_sort = th_err(condition.(cond{j})(indx1));
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
        mean_err.(cond{j}) = mean_err1;
        std_err.(cond{j}) = std_err1;
    end
    
    
    % plot data
    if plt
        sbpltcols = (length(cond));
        colr = lines(sbpltcols);
        
        if any(strcmp(cond{1},'s1')) || any(strcmp(cond{1},'bin1')) % stimtype OR tau bins
            for j = 1:length(cond)
                errorbar(tar_bins(2:end),mean_err.(cond{j}),std_err.(cond{j}),'Color',colr(j,:),'Marker','.');xlim([0 800]);ylim([-300 300]);
                suptitle(trials(1).prs.subject);title(['mean signed Angular error, target bins: ' num2str(binsize) ' deg']);
                xlabel('initial target angle [deg]');ylabel('mean error [deg]');grid on;
                hold on;
            end
            if any(strcmp(cond{1},'s1')) 
                legend(condSnames,'location','northwest');
            elseif any(strcmp(cond{1},'bin1'))
                legend(condJSnames,'location','northwest');
            end
        else
            % fix modalities, plot tau bins
            figure(h(1));
            for k = 1:length(condS)
                subplot(1,length(condS),k);
                indx = [];
                count = 0;
                for n = 1:length(cond)
                    if strcmp(cond{n}(1:2),condS{k})
                        count = count + 1;
                        indx = [indx n];
                        errorbar(tar_bins(2:end),mean_err.(cond{n}),std_err.(cond{n}),'Color',colr(count,:),'Marker','.');
                        hold on;
                    end
                    xlim([0 800]);ylim([-300 300]);grid on;
                    title(['mean signed Angular error, target bins: ' num2str(binsize) ' deg']);
                    xlabel('initial target angle [deg]');ylabel('mean error [deg]');
                end
                legend(leg_input{indx},'location','northwest');
                suptitle(trials(1).prs.subject);
            end
            % fix tau bins, plot modalities
            figure(h(2));
            for k = 1:length(condJS)
                subplot(1,length(condJS),k);
                indx = [];
                count = 0;
                for n = 1:length(cond)
                    if strcmp(cond{n}(3:end),condJS{k})
                        count = count + 1;
                        indx = [indx n];
                        errorbar(tar_bins(2:end),mean_err.(cond{n}),std_err.(cond{n}),'Color',colr(count,:),'Marker','.');
                        hold on;
                        count = count + 1;
                    end
                    xlim([0 800]);ylim([-300 300]);grid on;
                    title(['mean signed Distance error, target bins: ' num2str(binsize) ' cm']);
                    xlabel('initial target distance (cm)');ylabel('mean error (cm)');
                end
                legend(leg_input{indx},'location','northwest');
                suptitle(trials(1).prs.subject);
            end
        end
    end
            
end         
