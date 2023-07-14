function TrackingErrorsDifference(tracking)

%% Check average tracking errors difference (target vs stop position)

Nsubs = size(tracking,1);
Nstim = size(tracking,2);

for i = 1:Nsubs
%     if 0; figure('Name',['Subject ' num2str(i)],'numbertitle','off','position',[0 60 900 250]); end
    for s = 1:Nstim
        
        ts = tracking{i,s}.misc.ts;
        mu_targ = nanmean(tracking{i,s}.eyepos.error.targ.all);        sd_targ = nanstd(tracking{i,s}.eyepos.error.targ.all);
        mu_stop = nanmean(tracking{i,s}.eyepos.error.stop.all);        sd_stop = nanstd(tracking{i,s}.eyepos.error.stop.all);
        err_diff{i,s} = nanmean(tracking{i,s}.eyepos.error.stop.all - tracking{i,s}.eyepos.error.targ.all); 
        err_mu_diff{i,s} = (mu_stop - mu_targ);         err_mu_diff{i,s} = err_mu_diff{i,s} - err_mu_diff{i,s}(1);
        
%         if 0
%         subplot(1,Nstim,s); hold on;
%         plot(ts,mu_targ,'color',colr(s,:))
%         plot(ts,mu_stop,'color',colr(s,:)*0.5)
%         plot(ts,err_mu_diff{i,s},'r--','linewidth',1); xlim(xaxislim);
%         plot(ts,err_diff{i,s},'r','linewidth',3); 
%         xlabel(xaxislabel); ylabel('tracking error [deg]'); legend('target','stop','difference','location','northwest');
%         end
    end
end

% Plot
if strcmp(tracking{1}.misc.sampleflag,'time')
    xaxislim = [0 10];  xaxislabel = 'time [s]';
elseif strcmp(tracking{1}.misc.sampleflag,'distperc')
    xaxislim = [0 1];  xaxislabel = 'distance %';
elseif strcmp(tracking{1}.misc.sampleflag,'distance2end')
    xaxislim = [-600 0];  xaxislabel = 'distance to end [cm]';
end

colr = brewermap(Nstim,'Dark2');

if strcmp(tracking{1}.misc.sampleflag,'distance2end')
minlength = min(cellfun(@(x) length(x), err_mu_diff(:))) - 1;
err_mu_diff = cellfun(@(x) x(end-minlength:end)',err_mu_diff,'uniformoutput',false);
ts = ts(end-minlength:end);    
else
minlength = min(cellfun(@(x) length(x),err_mu_diff(:)));
err_mu_diff = cellfun(@(x) x(1:minlength)',err_mu_diff,'uniformoutput',false);
ts = ts(1:minlength);
end
figure('name','All Subjects: Difference between SPTE and TTE','numbertitle','off','position',[-10 690 970 290]);
for s = 1:Nstim
    subplot(1,Nstim,s); hold on;
    err_diff_mu(s,:) = nanmean([err_mu_diff{:,s}],2);
    err_diff_se(s,:) = nanstd([err_mu_diff{:,s}],[],2)./sqrt(Nsubs);
    
    plot(ts,[err_mu_diff{:,s}],'color',colr(s,:)); 
    shadedErrorBar(ts,err_diff_mu(s,:),err_diff_se(s,:),'lineprops',{'color',colr(s,:)}); hline(0,'k--');
    title('SPTE and TTE difference'); xlabel(xaxislabel); ylabel('tracking error difference [deg]');
end

%% Relationship between bias and Errors difference
if 0
[bias,skipsub] = get_bias(subject,params,0);
avg_err_diff = cellfun(@(x) mean(x(1:300)), err_mu_diff);

figure('name','Bias vs Errors difference','numbertitle','off'); hold on;
for s = 1:Nstim
    [rho(s),pval(s)] = corr(bias.r(:,s),avg_err_diff(:,s));
    plot(bias.r(:,s),avg_err_diff(:,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','none')   
end
xlabel('linear bias'); ylabel('average error difference'); hline(0,'k--'); vline(1,'k--'); xlim([0.2 1.4]);
end
