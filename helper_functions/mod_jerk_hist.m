function [avg,sig] = mod_jerk_hist(subject,N,logsc,plt)
%% Plot distribution of cumulative jerk of every trial and every subject conditioned on modalities
% logsc: 1 to plot in logscale, 0 for normal scale
if plt
figure;hold on;
whitebg([.5 .5 .5]);
end

stimtype = fieldnames(subject(1).stimtype);

stimnames = {'vestibular' , 'visual' , 'combined'};
colr = parula(length(stimnames));
for n = 1:length(stimtype)-1
    jerk = [];
    count = 0;
    for i = 1:length(subject)
        for j = 1:length(subject(i).stimtype.(stimtype{n}))
            count = count + 1;
            indx = subject(i).stimtype.(stimtype{n})(j);
            jerk(count) = subject(i).trials(indx).continuous.cumjerk(end);
        end
    end
    
    if logsc
        xbins = logspace(1,5,N);
        loglabl = 'LOG - ';
        logtitl = ' - logspace';
        jerklims = [0 10*10^4];
        % calculate mean and variance
        avg(1,n) = nanmean(log(jerk));
        sig(1,n) = nanstd(log(jerk));
    else
        xbins = N;
        loglabl = [];
        logtitl = [];
        jerklims = [0 2*10^4];
    end
    
    [y(n,:),x(n,:)] = hist(jerk,xbins);
    if plt
        plot(x(n,:),y(n,:)/count,'Color',colr(n,:),'LineWidth',2);grid on;
        xlabel([loglabl 'cumulative jerk^2' ]);xlim(jerklims);
        title(['Histogram of cumulative jerk^2 for different modalities' logtitl]);
        if logsc
            set(gca,'xScale','log');
            %         vline(exp(avg),'r');vline([exp(avg - sig) ; exp(avg + sig)],'k');
        end
        
        leg_input{n} = [stimnames{n}];
    end
end

if plt
    legend(leg_input{:});xlabel('cumulative jerk^2');ylabel('P(J^2|modality)');
    if length(subject) == 1
        suptitle(subject.name)
    end
end
