function plotcontrol2d(trials,tau)

%% Plot joystick control
% different colors for each JS coefficient
condnames = fieldnames(tau);
x = -1000:1000;
y = x;
dt = (1/60);
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
downsample = 3;
step = 1;
for i = 1:length(condnames)
    subplot(length(condnames),1,i)
    for j = 1:step:length(tau.(condnames{i}))
        n = tau.(condnames{i})(j);
        samples = 1:downsample:length(trials(n).mc.JS_X_Raw);
        colr = [];
        colr = parula(length(samples));
        vmax = trials(n).prs.vmax;
        wmax = trials(n).prs.wmax;
        for t = 1:length(samples)
            plot(trials(n).mc.JS_Yaw_Raw(samples(t))/wmax,trials(n).mc.JS_X_Raw(samples(t))/vmax,'.-','Color',colr(t,:));hold on;
        end
    end
    title(['JS position over trial: ' condition{i} ' : [' num2str(tau.(limnames{i})) ']']);hold off; xlim([-1.1 1.1]);ylim([-1.1 1.1]);xlabel('JS rotational input');ylabel('JS forward input');grid on;
    g = colorbar('Ticks',linspace(0,1,5),...
        'TickLabels',dt*linspace(0,length(trials(n).mc.JS_X_Raw),5));
    ylabel(g,'Real time of trial (s)')
end
suptitle(trials(1).prs.subject);