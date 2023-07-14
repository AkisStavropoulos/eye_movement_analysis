function scattermod_js(trials,stimtype,tau,plt)
%% Scatterplot of distance and angle of target and subject's stop (arena coordinates)
% for different Stimulus type conditions (vestibular, visual, combined)
% polar transformation
stnames = fieldnames(stimtype);
stnames = stnames(1:end-1);

taunames = fieldnames(tau);
ind_lim = [];
ind_bin = [];
for i = 1:length(taunames)
    if strcmp(taunames{i}(1),'l')
        ind_lim = [ind_lim i];
    elseif strcmp(taunames{i}(1),'b')
        ind_bin = [ind_bin i];
    end
end
taulimnames = taunames(ind_lim);
taubinnames = taunames(ind_bin);

xx = -1000:1000;
yy = xx;
conditionS = [{'vestibular'} {'visual'} {'combined'}];
for k = 1:length(taubinnames)
conditionJS{k} = ['tau ' taunames{ind_bin(k)}];
end
prs = [trials.prs];
r_tar = [prs.r_tar];
r_sub2 = [prs.r_sub2];
th_tar = [prs.th_tar];
th_sub2 = [prs.th_sub2];

colr = brewermap(length(taubinnames),'Dark2');
if plt
    
    for j = 1:length(stnames)
        figure;
        for k = 1:length(taubinnames)
            indx = intersect(stimtype.(stnames{j}),tau.(taubinnames{k}));
            
            subplot(1,2,1);hold on;
            x = r_tar(indx)';
            y = r_sub2(indx)';
            %         X = [x ones(length(x),1)];
            %         [b]=regress(y,X); 
            %         c = b(2);
            %         b = b(1);
            %         plot(x,x*b + c,'Color',colr(j,:),'LineWidth',2);
            
            plot(x,y,'Color',colr(k,:),'Marker','.','LineStyle','none');
            title('scatterplot of distance traveled');grid on;
            
            subplot(1,2,2); hold on;
            x = th_tar(indx)';
            y = th_sub2(indx)';
            %         X = [x ones(length(x),1)];
            %         [b]=regress(y,X);
            %         c = b(2);
            %         b = b(1);
            %         plot(x,x*b + c,'Color',colr(j,:),'LineWidth',2);
            
            plot(x,y,'Color',colr(k,:),'Marker','.','LineStyle','none');
            title('scatterplot of \theta');grid on;
        end
        subplot(1,2,1);hold on;plot(xx,yy,'r--');
        ylim([0 800]);xlim([0 800]);ylabel('distance of subject (cm)');xlabel('distance of target (cm)');hold off;
        legend(conditionJS);
        
        subplot(1,2,2);hold on;plot(xx,yy,'r--');plot(xx,-yy,'r--');
        ylim([-80 80]);xlim([-80 80]);ylabel('\theta subject (degrees)');xlabel('\theta target (degrees)');hold off;
        legend(conditionJS);
        
        suptitle([trials(1).prs.subject ' - ' conditionS{j} ' condition']);
        % obtain figures to rearrange them on the screen later
        fig(j) = gcf;
        fig(j).Units = 'centimeters';
        
    end
    % place figures on Right screen in vertical order
    count = 0;
    for i = 1:length(fig)
        figdim = fig(i).Position(3:4);
        fig(i).Position = [-1.5+(count*(figdim(1)-2)) 12 figdim];  % 50: RL position , 18: max top position
        count = count + 1;
    end
end
