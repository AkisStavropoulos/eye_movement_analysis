function [avg,sig] = mod_tau_jerk_hist(subject,N,logsc,plt)
%% Plot distribution of cumulative jerk of every trial and every subject conditioned on modalities
% logsc: 1 to plot in logscale, 0 for normal scale

stimtype = fieldnames(subject(1).stimtype);
binnames = fieldnames(subject(1).tau.bins);
limnames = fieldnames(subject(1).tau.lims);
stimnames = {'vestibular' , 'visual' , 'combined'};

colr = parula(length(limnames));
for s = 1:length(stimtype)-1
    if plt
        figure;hold on; % one figure for each modality
        whitebg([.5 .5 .5]);
    end

    for t = 1:length(binnames)
        
        jerk = [];
        count = 0;
        
        for i = 1:length(subject)
            
            indx = intersect(subject(i).stimtype.(stimtype{s}),subject(i).tau.bins.(binnames{t})); % intersection
            
            for j = 1:length(indx)
                count = count + 1;
                jerk(count) = subject(i).trials(indx(j)).continuous.cumjerk(end);
            end
            
        end
        
        if logsc
            xbins = logspace(1,5,N);
            loglabl = 'LOG - ';
            logtitl = ' - logspace';
            jerklims = [0 10*10^4];
            % calculate mean and variance
%             avg(1,s) = nanmean(log(jerk));
%             sig(1,s) = nanstd(log(jerk));
        else
            xbins = N;
            loglabl = [];
            logtitl = [];
            jerklims = [0 2*10^4];
        end
        
        [y{s,t},x{s,t}] = hist(jerk,xbins);
        if plt
            plot(x{s,t},y{s,t}/count,'Color',colr(t,:),'LineWidth',2);grid on;
            xlabel([loglabl 'cumulative jerk^2' ]);xlim(jerklims);
            title(['Histogram of cumulative jerk^2 for different modalities and taus' logtitl]);
            if logsc
                set(gca,'xScale','log');
                %         vline(exp(avg),'r');vline([exp(avg - sig) ; exp(avg + sig)],'k');
            end
            
            leg_input{s,t} = [stimnames{s} ' ' binnames{t} ' - [' num2str(limnames{t}) ']'];
        end
    end
    if plt
        legend(leg_input{s,:});xlabel('cumulative jerk^2');ylabel('P(J^2|modality,\tau)');
        if length(subject) == 1
            suptitle(subject.name)
        end
    end

end

