function [meanbias,sigma,b] = biasvardecomp(trials,errvec,tarloc,binsize,plt)
%% bias and variance for each condition
% plt: plot yes 1 or no 0
if nargin < 3
    plt = 0;
end

% create binned space
% bins of 50 x 50 cm
maxd = 600;
mind = 100;
maxtheta = 42.5;

maxx = ceil(maxd*sind(maxtheta));
minx = -maxx;

maxy = maxd;
miny = floor(mind*cosd(maxtheta));

totedgex = maxx-minx;
remx = binsize*ceil(totedgex/binsize) - totedgex;
spacex = minx:binsize:(maxx + remx);

totedgey = maxy-miny;
remy = binsize*ceil(totedgey/binsize) - totedgey;
spacey = miny:binsize:(maxy + remy);
%%
if ~isstruct(errvec)
    
    b = [];
    residerr = [];
    tarnum = [];
    count = 0;
    for i = 1:length(spacey)-1
        
        tempy = [tarloc(2,:) >= spacey(i)].*[tarloc(2,:) <= spacey(i+1)];
        indxy = find(tempy);
        count = count + length(indxy); % to check if all targets were included
        if indxy
            
            for j = 1:length(spacex)-1
                
                tempx = [tarloc(1,:) >= spacex(j)].*[tarloc(1,:) <= spacex(j+1)];
                indxx = find(tempx);
                temp = tempx.*tempy;
                indx = find(temp); % column index of targets that are in this bin
                if indx
                    tarnum1 = length(indx);
                    tarnum = [tarnum tarnum1];
                    b_bin = mean(errvec(:,indx),2); % sum the columns
                    b = [b b_bin];
                    
                    residerr_bin = errvec(:,indx) - b_bin;
                    residerr = [residerr residerr_bin];
                end
            end
        end
    end
    % weighted sums instead of means, because each bin contains different number of targets
    M = length(residerr);
    weights = tarnum/M;
    
    sigma = sum(sum(residerr.^2,1))/(M);
    sigma = sqrt(sigma);
    
    b_temp = b.*repmat(weights,2,1);
    meanbx = sum(b_temp(1,:));
    meanby = sum(b_temp(2,:));
    meanbias = [meanbx ; meanby];
    %%
else
    condnames = fieldnames(errvec);
    
    if ~any(strcmp(condnames,'s1'))
        lim = [];
        bin = [];
        for i = 1:length(condnames)
            if strcmp(condnames{i}(end-3:end-1),'lim')
                lim = [lim i];
            elseif strcmp(condnames{i}(end-3:end-1),'bin')
                bin = [bin i];
            end
        end
        limnames = condnames(lim);
        condnames = condnames(bin);
    end
    
    if plt
        h1 = figure;
        h2 = figure;
        % name the conditions for plots
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
    end
    
    for k = 1:length(condnames)
        
        b1 = []; b_temp = [];
        residerr = [];
        tarnum = [];
        count = 0;
        
        for i = 1:length(spacey)-1
            
            tempy = [tarloc.(condnames{k})(2,:) >= spacey(i)].*[tarloc.(condnames{k})(2,:) <= spacey(i+1)];
            indxy = find(tempy);
            count = count + length(indxy); % to check if all targets were included
            if indxy
                
                for j = 1:length(spacex)-1
                    
                    tempx = [tarloc.(condnames{k})(1,:) >= spacex(j)].*[tarloc.(condnames{k})(1,:) <= spacex(j+1)];
                    indxx = find(tempx);
                    temp = tempx.*tempy;
                    indx = find(temp); % column index of targets that are in this bin
                    if indx
                        tarnum1 = length(indx);
                        tarnum = [tarnum tarnum1];
                        b_bin = mean(errvec.(condnames{k})(:,indx),2); % sum the columns
                        b1 = [b1 b_bin];
                        
                        residerr_bin = errvec.(condnames{k})(:,indx) - b_bin;
                        residerr = [residerr residerr_bin];
                    end
                end
            end
        end
        % weighted sums instead of means, because each bin contains different number of targets
        M = length(residerr);
        weights = tarnum/M;
        
        sigma1 = sum(sum(residerr.^2,1))/(M);
        sigma1 = sqrt(sigma1);
        sigma.(condnames{k}) = sigma1;
        
        b_temp = b1.*repmat(weights,2,1);
        meanbx = sum(b_temp(1,:));
        meanby = sum(b_temp(2,:));
        meanbias1 = [meanbx ; meanby];
        meanbias.(condnames{k}) = meanbias1;
        
        b.(condnames{k}) = b_temp;
        % plot data
        if plt
            if any(strcmp(condnames,'bin1')) % taus
                % std
                figure(h1);
                colr = hsv(length(condition));
                cond = 1:length(condition);
                plot(cond(k),sigma1,'Color',colr(k,:),'Marker','x','LineStyle','None');
                ax = gca; ax.XTick = cond;
                ax.XTickLabel = {condition{:}};hold on;
                title('standard deviation of residual error');xlabel('condition');ylabel('standard deviation');xlim([0 4]);ylim([0 100]);
                legend(condition);grid on;
                % bias
                figure(h2);
                plot(meanbx,meanby,'Color',colr(k,:),'Marker','x','LineStyle','None');title(['bias for \tau bins']);hold on;
                xlabel('x axis (cm)');ylabel('y axis (cm)');xlim([-300 300]);ylim([-300 300]);
                hline(0);vline(0);legend(condition);grid on;
                
            elseif any(strcmp(condnames,'s1')) % stimtype
                % std
                figure(h1);
                colr = hsv(length(condition));
                cond = 1:length(condition);
                plot(cond(k),sigma1,'Color',colr(k,:),'Marker','x','LineStyle','None');
                ax = gca; ax.XTick = cond;
                ax.XTickLabel = {condition{:}};hold on;
                title('standard deviation of residual error');xlabel('condition');ylabel('standard deviation');xlim([0 4]);ylim([0 100]);
                legend(condition);grid on;
                % bias
                figure(h2);
                plot(meanbx,meanby,'Color',colr(k,:),'Marker','x','LineStyle','None');title(['bias for modalities']);hold on;
                xlabel('x axis (cm)');ylabel('y axis (cm)');xlim([-300 300]);ylim([-300 300]);
                hline(0);vline(0);legend(condition);grid on;
                
            elseif any(strcmp(condnames,'s1bin1')) % stimtau
                
                % taus across modalities
                % std
                figure(h1);
                condS = unique(conditionS);
                condJS = unique(conditionJS);
                colr1 = repmat(hsv(length(condJS)),length(condJS),1);
                cond = sort(repmat([1:length(conditionStim)],1,length(condJS)));
                subplot(2,1,1);plot(cond(k),sigma1,'Color',colr1(k,:),'Marker','x','LineStyle','None');
                ax = gca; ax.XTick = 1:length(conditionStim);
                ax.XTickLabel = {conditionStim{:}};hold on;
                title('standard deviation of residual error');ylabel('standard deviation (cm)');xlim([0 4]);ylim([0 100]);
                legend(condJS{:});grid on;
                % bias
                figure(h2);
                for i = 1:length(condS)
                    if strcmp(condnames{k}(1:2),condS{i})
                        subplot(2,length(condS),i);hold on;
                        plot(meanbx,meanby,'Color',colr1(k,:),'Marker','x','LineStyle','None');title(['bias for ' conditionStim{i}]);hold on;
                        xlabel('x axis (cm)');ylabel('y axis (cm)');xlim([-300 300]);ylim([-300 300]);
                        hline(0);vline(0);legend(condJS{:});grid on;
                    end
                end
                % modalities across taus
                
                %std
                figure(h1);
                colr2_temp = hsv(length(conditionStim));
                colr2 = [];
                for n = 1:size(colr2_temp,1)
                    colr2 = [colr2 ; repmat(colr2_temp(n,:),length(conditionStim),1)];
                end
                cond = repmat([1:length(condJS)],1,length(conditionStim));
                subplot(2,1,2); h(k) = plot(cond(k),sigma1,'Color',colr2(k,:),'Marker','x','LineStyle','None');
                ax = gca; ax.XTick = 1:length(condJS);
                ax.XTickLabel = {condJS{:}};xlim([0 4]);ylim([0 100]);grid on;
                ylabel('standard deviation (cm)');hold on;legend(h(1:3:end),'vestibular','visual','combined','location','northwest');
                % bias
                figure(h2);
                for i = 1:length(condJS)                   
                    if strcmp(condnames{k}(3:end),condJS{i}(end-3:end))
                        subplot(2,length(condJS),i + length(condS));hold on;
                        p(k) = plot(meanbx,meanby,'Color',colr2(k,:),'Marker','x','LineStyle','None');title(['bias for ' conditionJS{i}]);hold on;
                        xlabel('x axis (cm)');ylabel('y axis (cm)');xlim([-300 300]);ylim([-300 300]);
                        hline(0);vline(0);legend(conditionStim{:});grid on;
                    end
                end
            end
        end
    end
figure(h1);suptitle([trials(1).prs.subject ' - target bins: ' num2str(binsize) ' cm^2']);
figure(h2);suptitle([trials(1).prs.subject ' - target bins: ' num2str(binsize) ' cm^2']);
end  
