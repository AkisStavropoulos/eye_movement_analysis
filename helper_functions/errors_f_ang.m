function errors_f_ang(trials,condition,r_tar,r_sub,theta_tar,theta_sub)

%% plot radial and angular error as a function of target distance

if isempty(condition) % for all trials
    
    %     for i = 1:length(trials)
    %         taus(i) = trials(i).prs.tau;
    %     end
    pos = find(theta_tar >= 0);
    neg = find(theta_tar < 0);
    
    % radial error
    x = abs(theta_tar)';
    X = [x ones(length(x),1)];
    y = (r_sub - r_tar)';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    figure;
    subplot(2,1,1);plot(theta_tar(pos),r_sub(pos)-r_tar(pos),'.');
    hold on;plot(-theta_tar(neg),r_sub(neg)-r_tar(neg),'m.');
    hline(0);ylim([-300 300]);xlim([min(x) 50]);
    hold on;plot(x,x*b + c,'r');
    xlabel('target angle (deg)');ylabel('radial error (cm)');title('radial error as a function of target angle')
    legend('right turns','left turns','regression');
    % angular error
    x = abs(theta_tar');
    X = [x ones(length(x),1)];
    y = (abs(theta_sub) - abs(theta_tar))';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    subplot(2,1,2);plot(theta_tar(pos),abs(theta_sub(pos))-abs(theta_tar(pos)),'.');
    hold on;plot(-theta_tar(neg),abs(theta_sub(neg))-abs(theta_tar(neg)),'m.');
    hline(0);ylim([-30 30]);
    hold on;plot(x,x*b + c,'r');legend('right turns','left turns','regression');
    xlabel('target angle (deg)');ylabel('rotational error (deg)');title('angular error as a function of target angle');
    suptitle(trials(1).prs.subject);
else
    condnames = fieldnames(r_tar);
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
        conditionJS = unique(conditionJS);
        count = 0;
        for j = 1:length(conditionStim)
            for n = 1:length(conditionJS)
                count = count + 1;
                cond{count} = [conditionStim{j} ' and ' conditionJS{n}];
            end
        end
    else
        if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2'))
            cond = [{'vestibular'} {'visual'} {'combined'} ];
        elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
            for i = 1:length(condnames)
                cond{i} = ['\tau ' condnames{i}];
            end
        end
    end
    %%
    for j = 1:length(condnames)
        %         taus = [];
        x = [];
        X = [];
        y = [];
        pos = [];
        neg = [];
        %         for i = 1:length(condition.(condnames{j}))
        %             taus(i) = trials(condition.(condnames{j})(i)).prs.tau;
        %         end
        pos = find(theta_tar.(condnames{j}) >= 0);
        neg = find(theta_tar.(condnames{j}) < 0);
        figure;
        % radial error
        x = abs(theta_tar.(condnames{j}))';
        X = [x ones(length(x),1)];
        y = (r_sub.(condnames{j}) - r_tar.(condnames{j}))';
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        subplot(2,1,1);plot(theta_tar.(condnames{j})(pos), r_sub.(condnames{j})(pos) - r_tar.(condnames{j})(pos),'.');
        hold on;plot(-theta_tar.(condnames{j})(neg), r_sub.(condnames{j})(neg) - r_tar.(condnames{j})(neg),'m.');
        hline(0);ylim([-300 300]);xlim([min(x) 50]);
        hold on;plot(x,x*b + c,'r');legend('right turns','left turns','regression');
        xlabel('target angle (deg)');ylabel('radial error (cm)');title('radial error as a function of target angle')
        
        % angular error
        x = abs(theta_tar.(condnames{j}))';
        X = [x ones(length(x),1)];
        y = (abs(theta_sub.(condnames{j})) - abs(theta_tar.(condnames{j})))';
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        subplot(2,1,2);plot(theta_tar.(condnames{j})(pos), abs(theta_sub.(condnames{j})(pos)) - abs(theta_tar.(condnames{j})(pos)),'.');
        hold on;plot(-theta_tar.(condnames{j})(neg), abs(theta_sub.(condnames{j})(neg)) - abs(theta_tar.(condnames{j})(neg)),'m.');
        hline(0);ylim([-30 30]);xlim([min(x) 50]);
        hold on;plot(x,x*b + c,'r');legend('right turns','left turns','regression');
        xlabel('target angle (deg)');ylabel('rotational error (deg)');title('angular error as a function of target angle');
        suptitle([trials(1).prs.subject ' - ' cond{j}]);
        % obtain figures to rearrange them on the screen later
        fig(j) = gcf;
        fig(j).Units = 'centimeters';
        
        
    end
    % place figures on Right screen in vertical order
    if any(strcmp(condnames,'s1'))
        count = 0;
        for i = 1:length(fig)
            figdim = fig(i).Position(3:4);
            fig(i).Position = [-.5+(count*(figdim(1)-.5)) 2 figdim]; % 0: max left position , 2: max down position
            count = count + 1;
        end
    else
        figposx =  [-1.5 -1.5 -1.5 ...
            -1.5 -1.5 -1.5...
            50 50 50 ...
            50 50 50];
        figposy =  [15 15 15 ...
            1 1 1 ...
            15 15 15 ...
            3 3 3];
        count1 = 0;
        for n = 1:length(fig)
            figdim = fig(n).Position(3:4);
            fig(n).Position = [figposx(n)+(count1*(figdim(1)-2)) figposy(n) figdim]; % figposx: initial LR position , figposy: initial top-down position
            %         fig(i).Position = [-1.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
            count1 = count1 + 1;
            if ~rem(n,3); count1 = 0; end
        end
    end
end