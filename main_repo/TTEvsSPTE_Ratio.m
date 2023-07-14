function [kernel,ts] = TTEvsSPTE_Ratio(tracking,stability_constant)

%% Ratio between SPTE and TTE
% stability_constant: ratio stability constant

Nsubs = size(tracking,1);
Nstim = size(tracking,2);

colr = brewermap(Nstim+1,'Set1'); colr(3,:) = []; % red, blue, purple
if strcmp(tracking{1}.misc.sampleflag,'time')
    xaxislim = [0 12];  xaxislabel = 'time [s]';
elseif strcmp(tracking{1}.misc.sampleflag,'distperc')
    xaxislim = [0 1];  xaxislabel = 'distance %';
elseif strcmp(tracking{1}.misc.sampleflag,'distance2end')
    xaxislim = [-600 0];  xaxislabel = 'distance to end [cm]';
end

% extract minimum length of data series
if strcmp(tracking{1}.misc.sampleflag,'time')
    minlength = 12*60;
else
    minlength = min(cellfun(@(x) length(x.misc.ts), tracking(:)));
end

% get ratio
for i = 1:Nsubs
    for s = 1:Nstim
        ts = tracking{i,s}.misc.ts;
        SPTE = tracking{i,s}.eyepos.error.stop.all;
        TTE = tracking{i,s}.eyepos.error.targ.all;
        b{s}(i,:) = arrayfun(@(t) nanmean((stability_constant+TTE(:,t))./(stability_constant+SPTE(:,t))), 1:minlength);
                
    end
end

wn = 5;
b = cellfun(@(x) medfilt1(x,wn,[],2), b, 'uniformoutput',false);
kernel = b;

figure('name','All Subjects: Ratio between SPTE and TTE','numbertitle','off','position',[-10 690 970 290]);
for s = 1:Nstim
    subplot(1,Nstim,s); hold on; plot(ts(1:minlength),kernel{s},'color',colr(s,:));
    shadedErrorBar(ts(1:minlength),mean(kernel{s}),std(kernel{s})./sqrt(Nsubs),'lineprops',{'linewidth',3,'color',colr(s,:)}); axis([xaxislim 0 3]);
    hline(1,'k--'); xlabel(xaxislabel); ylabel('mean ratio kernel (TTE/SPTE)');
end


%% Relationship between bias and average Ratio of errors
if 0
[bias,skipsub] = get_bias(subject,params,0);
avg_coef = cell2mat(cellfun(@(x) mean(x(:,1:300),2), b,'uniformoutput',false));

rho = []; pval = [];
figure('name','Bias vs Errors Ratio','numbertitle','off'); hold on;
for s = 1:Nstim
    [rho(s),pval(s)] = corr(bias.r(:,s),avg_coef(:,s));
    plot(bias.r(:,s),avg_coef(:,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')   
end
xlabel('linear bias'); ylabel('average ratio of errors'); hline(1,'k--'); vline(1,'k--'); xlim([0.2 1.4]);
if 0
figure; hold on; 
for s = 1:Nstim
    plot(avg_coef(:,s),'color',colr(s,:));  
    plot(bias.r(:,s),'color',colr(s,:),'linewidth',3);
end
hline(1,'k--');
end
end
