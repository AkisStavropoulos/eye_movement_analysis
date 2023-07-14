function scatterall(trials,stimtype,plt)
%% Scatterplot of distance and angle of target and subject's stop (arena coordinates)
% for different Stimulus type conditions (vestibular, visual, combined)
% polar transformation
stnames = fieldnames(stimtype);
stnames = stnames(1:3);
xx = -1000:1000;
yy = xx;
condition = [{'vestibular'} {'visual'} {'combined'} {'joystick'}];

prs = [trials.prs];
r_tar = [prs.r_tar];
r_sub = [prs.r_sub];
th_tar = [prs.th_tar];
th_sub = [prs.th_sub];

colr = brewermap(length(stnames),'Dark2');
if plt
    figure;
    
    for j = 1:length(stnames)
        subplot(1,2,1); hold on; axis equal;
        x = r_tar(stimtype.(stnames{j}))';
        y = r_sub(stimtype.(stnames{j}))';
        X = [x zeros(length(x),1)];
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        plot([0:800],[0:800]*b + c,'Color',colr(j,:),'LineWidth',2,'HandleVisibility','off');
        
        plot(x,y,'Color',colr(j,:),'Marker','.','MarkerSize',2,'LineStyle','none');
        title('scatterplot of distance traveled');%grid on;
        
        subplot(1,2,2); hold on; axis equal;
        x = th_tar(stimtype.(stnames{j}))';
        y = th_sub(stimtype.(stnames{j}))';
        X = [x zeros(length(x),1)];
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        plot([-80:80],[-80:80]*b + c,'Color',colr(j,:),'LineWidth',2,'HandleVisibility','off');
                
       plot(x,y,'Color',colr(j,:),'Marker','.','MarkerSize',2,'LineStyle','none');
        title('scatterplot of \theta');%grid on;
    end
end

subplot(1,2,1);hold on;plot(xx,yy,'k--');
ylim([0 800]);xlim([0 800]);ylabel('distance of subject (cm)');xlabel('distance of target (cm)');hold off;
        legend(condition{1:length(stnames)},'Location','northwest');

subplot(1,2,2);hold on;plot(xx,yy,'k--');%plot(xx,-yy,'r--');
ylim([-80 80]);xlim([-80 80]);ylabel('\theta subject (degrees)');xlabel('\theta target (degrees)');hold off;
       legend(condition{1:length(stnames)},'Location','northwest');

suptitle([trials(1).prs.subject ' - ' 'all modalities']);

