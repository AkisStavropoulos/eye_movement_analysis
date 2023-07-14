function [epsilon,sd_epsilon] = epsilonerror1(err,r_tar)
%% epsilon error
% error normalized by initial distance
% epsilon = 1/N * sum(error./r_tar)

if ~isstruct(r_tar)
    
    N = length(r_tar);
    epsilon = (1/N)*sum(err./r_tar);
    sd_epsilon = std(err./r_tar);
    
else
    condnames = fieldnames(err);
    
    % name the conditions
    % name the conditions for the legend
    if any(strcmp(condnames,'s1bin1')) % check for intersection first
        for i = 1:length(condnames)
            conditionS{i} = condnames{i}(1:2);
            conditionJS{i} = condnames{i}(3:end);
        end
        if any(strcmp(conditionS,'s1')) && any(strcmp(conditionS,'s2'))
            conditionStim = [{'vestibular'} {'visual'} {'combined'} ];
        end
        if any(strcmp(conditionJS,'bin1')) && any(strcmp(conditionJS,'bin2'))
            for i = 1:length(condnames)
                conditionJS{i} = ['\tau ' conditionJS{i}];
            end
        end
        count = 0;
        for j = 1:length(conditionStim)
            for n = 1:length(conditionJS)
                count = count + 1;
                condition{count} = [conditionStim{j} ' and ' conditionJS{n}];
            end
        end
    else
        if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2'))
            condition = [{'vestibular'} {'visual'} {'combined'} ];
        elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
            for i = 1:length(condnames)
                condition{i} = ['\tau ' condnames{i}];
            end
        end
    end
    figure;
    % compute epsilon for each condition
    for j = 1:length(condnames)
        
        N = length(err.(condnames{j}));
        epsi(j) = (1/N)*sum(err.(condnames{j})./r_tar.(condnames{j}));
        epsilon.(condnames{j}) = epsi(j);
        sd_epsi(j) = std(err.(condnames{j})./r_tar.(condnames{j}));
        sd_epsilon.(condnames{j}) = sd_epsi(j);
        
        % plot data
        if any(strcmp(condnames{1},'s1')) || any(strcmp(condnames{1},'bin1')) % stimtype OR tau bins
            
            if any(strcmp(condnames{1},'s1'))
                
                colr = ['k' 'b' 'm'];
                cond = 1:length(condnames);
                errorbar(cond(j),epsi(j),sd_epsi(j),[colr(j) 'x']);
%                 plot(cond(j),epsi(j),[colr(j) 'x']);
                ax = gca; ax.XTick = cond;
                ax.XTickLabel = {condition{:}};hold on;
                title('normalized error by initial distance');xlabel('condition');ylabel('epsilon value');
                xlim([0 4]);ylim([0 0.5]);legend(condition);
            else
                colr = lines(length(condnames));
                cond = 1:length(condnames);
                errorbar(cond(j),epsi(j),sd_epsi(j),'Color',colr(j,:),'Marker','x','LineStyle','None');
%                 plot(cond(j),epsi(j),'Color',colr(j,:),'Marker','x','LineStyle','None');
                ax = gca; ax.XTick = cond;
                ax.XTickLabel = {condition{:}};hold on;
                title('normalized error by initial distance');xlabel('condition');ylabel('epsilon value');
                xlim([0 4]);ylim([0 0.5]);legend(condition);
            end
            
        else
            % tau bins across modalities
            condJS = unique(conditionJS);
            colr1 = repmat(hsv(length(condJS)),length(condJS),1);
            cond = sort(repmat([1:length(conditionStim)],1,length(condJS)));
            subplot(2,1,1);
            errorbar(cond(j),epsi(j),sd_epsi(j),'Color',colr1(j,:),'Marker','x','LineStyle','None');
%             plot(cond(j),epsi(j),'Color',colr1(j,:),'Marker','x','LineStyle','None');
            ax = gca; ax.XTick = 1:length(conditionStim);  
            ax.XTickLabel = {conditionStim{:}};hold on;
            title('normalized error by initial distance');ylabel('epsilon value');xlim([0 4]);ylim([0 0.5]);
            legend(condJS{:});
            
            % modalities across tau bins
            colr2_temp = hsv(length(conditionStim));
            colr2 = [];
            for n = 1:size(colr2_temp,1)
            colr2  = [colr2 ; repmat(colr2_temp(n,:),length(conditionStim),1)];
            end
            
            cond = repmat([1:length(condJS)],1,length(conditionStim));
            subplot(2,1,2);
            h(j) = errorbar(cond(j),epsi(j),sd_epsi(j),'Color',colr2(j,:),'Marker','x','LineStyle','None');
%             h(j) = plot(cond(j),epsi(j),'Color',colr2(j,:),'Marker','x','LineStyle','None');
            ax = gca; ax.XTick = 1:length(condJS);
            ax.XTickLabel = {condJS{:}};xlim([0 4]);ylim([0 0.5]);
            ylabel('epsilon value');hold on;legend(h(1:3:end),conditionStim{:},'location','northwest');

            
        end
    end
end
hold off;