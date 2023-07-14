function plotvel2D(trials,tau)

%% Plot trajectories
% different colors for each JS coefficient
condnames = fieldnames(tau);

dt = trials(1).continuous.ts(2) - trials(1).continuous.ts(1);


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
downsample = 8;
step = 2;
for i = 1:length(condnames)
    subplot(length(condnames),1,i)
    for j = 1:step:length(tau.(condnames{i}))
        n = tau.(condnames{i})(j);
        samples = 1:downsample:length(trials(n).continuous.v);
        colr = [];
            colr = parula(length(samples));
            for t = 1:length(samples)
                plot(trials(n).continuous.w(samples(t)),trials(n).continuous.v(samples(t)),'.','Color',colr(t,:));hold on;
            end
    end
    title(['Velocities over trial: ' condition{i} ' : [' num2str(tau.(limnames{i})) ']']);hold off; ylim([-20 100]);xlim([-15 15]);xlabel('angular velocity (deg/s)');ylabel('forward velocity (cm/s)');grid on;
    g = colorbar('Ticks',linspace(0,1,5),...
        'TickLabels',dt*linspace(0,length(trials(n).continuous.v),5));
    ylabel(g,'Real time of trial (s)')
    
end
suptitle(trials(1).prs.subject);



