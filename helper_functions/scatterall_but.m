function scatterall_but(trials,stimtype,plt)
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
r_sub_but = [prs.r_sub_but];
th_tar = [prs.th_tar];
th_sub_but = [prs.th_sub_but];

colr = brewermap(length(stnames),'Dark2');
if plt
    figure;
    
    for j = 1:length(stnames)
        subplot(1,2,1);hold on;
        x = r_tar(stimtype.(stnames{j}))';
        y = r_sub_but(stimtype.(stnames{j}))';
%         X = [x ones(length(x),1)];
%         [b]=regress(y,X);
%         c = b(2);
%         b = b(1);
%         plot(x,x*b + c,'Color',colr(j,:),'LineWidth',2);
        
        plot(x,y,'Color',colr(j,:),'Marker','.','MarkerSize',8,'LineStyle','none');
        title('scatterplot of distance traveled');grid on;
        
        subplot(1,2,2); hold on;
        x = th_tar(stimtype.(stnames{j}))';
        y = th_sub_but(stimtype.(stnames{j}))';
%         X = [x ones(length(x),1)];
%         [b]=regress(y,X);
%         c = b(2);
%         b = b(1);
%         plot(x,x*b + c,'Color',colr(j,:),'LineWidth',2);
                
       plot(x,y,'Color',colr(j,:),'Marker','.','MarkerSize',8,'LineStyle','none');
        title('scatterplot of \theta');grid on;
    end
end

subplot(1,2,1);hold on;plot(xx,yy,'r--');
ylim([0 800]);xlim([0 800]);ylabel('distance of subject (cm)');xlabel('distance of target (cm)');hold off;
        legend(condition{1:length(stnames)});

subplot(1,2,2);hold on;plot(xx,yy,'r--');plot(xx,-yy,'r--');
ylim([-80 80]);xlim([-80 80]);ylabel('\theta subject (degrees)');xlabel('\theta target (degrees)');hold off;
       legend(condition{1:length(stnames)});

suptitle([trials(1).prs.subject ' - ' 'all modalities' ' - 1st button push']);

        fig = gcf;
        fig.Units = 'centimeters';
        

    % place figures on bottom center of screen
    count = 0;
    for i = 1:length(fig)
        figdim = fig(i).Position(3:4);
        fig(i).Position(2) = 1;
        count = count + 1;
    end


