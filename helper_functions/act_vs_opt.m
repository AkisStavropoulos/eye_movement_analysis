function [d_s,d_brake,d_T] = act_vs_opt(T,s,T_trial,s_trial,condition,trials)
%% plot accross modalities
% comparison of optimal vs actual data
% d_s = difference between actual and optimal switchtime (s_trial - s)
% d_brake = difference between actual and optimal brake duration [(T_trial - s_trial) - (T - s)]
% d_T = difference between actual and optimal trial duration (T_trial - T)
condnames = fieldnames(condition);

if any(strcmp(condnames,'s1bin1')) % check for intersection first
    ind_keep = [];
    for i = 1:length(condnames) % throw out the limits fileds
        if ~strcmp(condnames{i}(1),'l')
            ind_keep = [ind_keep i];
        end
    end
    condnames = condnames(ind_keep);
    tau_bins = (length(condnames)-1)/3;
    
    for i = 1:length(condnames)-1
        conditionS{i} = condnames{i}(1:2);
        conditionJS{i} = condnames{i}(3:end);
    end
    if any(strcmp(conditionS,'s1')) && any(strcmp(conditionS,'s2'))
        conditionStim = [{'vestibular'} {'visual'} {'combined'} ];
    end
    if any(strcmp(conditionJS,'bin1')) && any(strcmp(conditionJS,'bin2'))
        for i = 1:length(condnames)-1
            conditionJS{i} = ['\tau ' conditionJS{i}];
        end
        conditionJS = unique(conditionJS);
    end
    count = 0;
    for j = 1:length(conditionStim)
        for n = 1:length(conditionJS)
            count = count + 1;
            cond{count} = [conditionStim{j} ' and ' conditionJS{n}];
        end
    end
else
    if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2')) % stimtype
        cond = [{'vestibular'} {'visual'} {'combined'} {'all trials'}];
        
    elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2')) % taus
        lim = [];
        bin = [];
        count = 1;
        for i = 1:length(condnames)-1
            
            if strcmp(condnames{i}(end-3:end-1),'lim')
                lim = [lim i];
            elseif strcmp(condnames{i}(end-3:end-1),'bin')
                bin = [bin i];
                cond{count} = ['\tau ' condnames{bin(count)}];
                count = count + 1;
            end
        end
        limnames = condnames(lim);
            condnames1 = condnames(bin);
        condnames = condnames1;
    end
end


d_s = [];
d_T = [];
d_brake = [];
sbpltcols = length(condnames);
for n = 1:length(condnames)*2
    h(n) = figure;
end
for j = 1:length(cond)
    figure(h(j));
    indx = condition.(condnames{j});
    
    % switch time comparison
    d_s.(condnames{j}) = s_trial(indx)-s(indx);
    x = [1:length(d_s.(condnames{j}))]';
    X = [x ones(length(x),1)];
    y = d_s.(condnames{j})';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    subplot(1,3,1);plot(d_s.(condnames{j}));title(['actual - optimal switchtime']);
    hold on;plot(x,b*x + c, 'r');xlabel('trials');ylabel('switchtime diff (s)');
    ylim([-10 10]);grid on;
    % brake time comparison
    d_brake.(condnames{j}) = (T_trial(indx) - s_trial(indx)) - (T(indx) - s(indx));
    x = [1:length(d_brake.(condnames{j}))]';
    X = [x ones(length(x),1)];
    y = d_brake.(condnames{j})';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    subplot(1,3,3);plot(d_brake.(condnames{j}));title(['actual - optimal Tbrake']);
    hold on;plot(x,b*x + c, 'r');xlabel('trials');ylabel('Brake time diff (s)');
    ylim([-10 10]);grid on;
    % trial duration comparison
    d_T.(condnames{j}) = T_trial(indx)-T(indx);
    x = [1:length(d_T.(condnames{j}))]';
    X = [x ones(length(x),1)];
    y = d_T.(condnames{j})';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    subplot(1,3,2);plot(d_T.(condnames{j}));title(['actual - optimal trial duration']);
    hold on;plot(x,b*x + c, 'r');xlabel('trials');ylabel('Trial duration diff (s)');
    ylim([-10 10]);grid on;

    suptitle([trials(1).prs.subject ' - ' cond{j}]);
    
    figure(h(j+length(condnames)));
    colr = parula(length(indx));
    for n = 1:length(indx)
        
        % 2D plot of switchtime and brake time
        vline(0);hline(0);plot(d_brake.(condnames{j})(n),d_s.(condnames{j})(n),'.','Color',colr(n,:));hold on;grid on;
        title('Switchtime difference over Tbrake difference');xlabel('Tbrake diff');ylabel('switchtime diff');
        ylim([-10 10]);xlim([-10 10]);suptitle([trials(1).prs.subject ' - ' cond{j}]);
    end
    g = colorbar('Ticks',linspace(0,1,5),...
        'TickLabels',linspace(0,length(indx),5));
    ylabel(g,'Experience in trials')
end
