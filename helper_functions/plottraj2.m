function plottraj2(trials,tau)

%% Plot trajectories
% different colors for each JS coefficient
condnames = fieldnames(tau);
x = -1000:1000;
y = x;

if any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
    lim = [];
    bin = [];
    for i = 1:length(condnames)-1
        if strcmp(condnames{i}(end-3:end-1),'lim')
            lim = [lim i];
        elseif strcmp(condnames{i}(end-3:end-1),'bin')
            bin = [bin i];
        end
        
    end
    limnames = condnames(lim);
    condnames = condnames(bin);
end
for i = 1:length(condnames)
    condition{i} = ['\tau ' condnames{i}];
end
    figure;

for i = 1:length(condnames)
    colr = lines(length(condnames));
    subplot(length(condnames),1,i)
    for j = 1:length(tau.(condnames{i}))
        indx = tau.(condnames{i})(j);
        plot(trials(indx).continuous.xmp,trials(indx).continuous.ymp,'Color',colr(i,:));hold on;
    end
    title(['trajectories: ' condition{i} ' : [' num2str(tau.(limnames{i})) ']']);hold off; ylim([-100 700]);xlim([-400 400]);xlabel('x');ylabel('y');grid on;
end
suptitle(trials(1).prs.subject);