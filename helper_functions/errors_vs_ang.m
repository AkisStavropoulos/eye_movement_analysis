function errors_vs_ang(trials,condition,r_tar,r_sub,theta_tar,theta_sub)
%% plot radial and angular error as a function of target angle

if isempty(condition) % for all trials
    
    pos = find(theta_tar >= 0);
    neg = find(theta_tar < 0);

    % radial error
    x = [theta_tar(pos) -theta_tar(neg)]';
    y = [r_sub(pos)-r_tar(pos) r_sub(neg)-r_tar(neg)]';
    X = [x ones(length(x),1)];
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    figure;
    subplot(2,1,1);plot(x,y,'b.');
    hline(0);ylim([-300 300]);xlim([min(x) 40]);
    hold on;plot(x,x*b + c,'r','LineWidth',2);
    xlabel('target angle (deg)');ylabel('radial error (cm)');title('radial error as a function of angle')
    grid on;
    % angular error
    x = [theta_tar(pos) -theta_tar(neg)]';
    y = [theta_sub(pos)-theta_tar(pos) -(theta_sub(neg)-theta_tar(neg))]';
    X = [x ones(length(x),1)];
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    
    subplot(2,1,2);plot(x,y,'b.');
    hline(0);ylim([-30 30]);
    hold on;plot(x,x*b + c,'r','LineWidth',2);
    xlabel('target angle (deg)');ylabel('rotational error (deg)');title('angular error as a function of angle');
    grid on;
    suptitle([trials(1).prs.subject ' - all trials']);
else
    condnames = fieldnames(theta_tar);
    % name the conditions for the legend
    if any(strcmp(condnames,'s1bin1')) % check for intersection first
        for i = 1:length(condnames)
            condS{i} = condnames{i}(1:2);
            condJS{i} = condnames{i}(3:end);
        end
        condS = unique(condS);
        condJS = unique(condJS);
        if any(strcmp(condS,'s1')) && any(strcmp(condS,'s2'))
            condSnames = [{'vestibular'} {'visual'} {'combined'} ];
        end
        if any(strcmp(condJS,'bin1')) && any(strcmp(condJS,'bin2'))
            for i = 1:length(condJS)
                condJSnames{i} = ['\tau ' condJS{i}];
            end
        end
        condJS = unique(condJS);
    else
        if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2'))
            for i = 1:length(condnames)
                condS{i} = condnames{i};
            end
            condSnames = [{'vestibular'} {'visual'} {'combined'} ];
            condJS = {[]};
            condJSnames = {[]};
        elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
            for i = 1:length(condnames)
                condJS{i} = condnames{i};
                condJSnames{i} = ['\tau ' condJS{i}];
            end
            condS = {[]};
            condSnames = {[]};
        end
    end
  %% Plot
  % create colormap
    indx = 0;
  if length(condJS) > 1 % for stimtau
      colr = colormap(lines); 
%       colr = colr(round(linspace(1,length(colr),length(condJS))),:);
  elseif length(condJS) == 1 % for stimtype
      figure;
      indx = indx + 1;
      colr = colormap(lines);
%       colr = colr(round(linspace(1,length(colr),length(condS))),:);
  end
  % start plotting
  leg_inp = [];
  for i = 1:length(condS)
      if length(condJS)>1
          figure;
          indx = indx + 1;
      end
      if length(condJS) > 1 && length(condS) > 1% for stimtau
          leg_inp = [];
      end
      count = 0;
      for j = 1:length(condJS)
          x = [];
          X = [];
          y = [];
          
          pos = find(theta_tar.([condS{i} condJS{j}]) >= 0);
          neg = find(theta_tar.([condS{i} condJS{j}]) < 0);
          % create legend
          leg_inp_temp = {[condSnames{i} ' ' condJSnames{j}]};
          leg_inp = [leg_inp ; repmat(leg_inp_temp,2,1)];
          leg_input{i} = leg_inp(:);
          % color counting
          if length(condS) == 1 % taus
              count = j;
          elseif length(condS) > 1 && length(condJS) == 1 % stimtype
              count = i;
          elseif length(condS) > 1 && length(condJS) > 1  % stimtau
              count = count + 1;
          end
          % radial error
          subplot(2,1,1);hold on;
          x = [theta_tar.([condS{i} condJS{j}])(pos) -theta_tar.([condS{i} condJS{j}])(neg)]';
          y = [r_sub.([condS{i} condJS{j}])(pos) - r_tar.([condS{i} condJS{j}])(pos)  r_sub.([condS{i} condJS{j}])(neg) - r_tar.([condS{i} condJS{j}])(neg)]';
          X = [x ones(length(x),1)];
          [b]=regress(y,X);
          c = b(2);
          b = b(1);
          plot(x, y,'.','Color',colr(count,:));
          hline(0);ylim([-300 300]);xlim([min(x) 50]);
          hold on;plot(x,x*b + c,'Color',colr(count,:),'LineWidth',2);
          xlabel('target angle (deg)');ylabel('radial error (cm)');title('radial error as a function of angle')
          legend(leg_input{i});grid on;
          % angular error
          subplot(2,1,2);hold on;
          x = [theta_tar.([condS{i} condJS{j}])(pos) -theta_tar.([condS{i} condJS{j}])(neg)]';
          y = [theta_sub.([condS{i} condJS{j}])(pos) - theta_tar.([condS{i} condJS{j}])(pos)...
              -(theta_sub.([condS{i} condJS{j}])(neg) - theta_tar.([condS{i} condJS{j}])(neg))]';
          X = [x ones(length(x),1)];
          [b]=regress(y,X);
          c = b(2);
          b = b(1);
          plot(x, y,'.','Color',colr(count,:));
          hline(0);ylim([-30 30]);xlim([min(x) 50]);
          hold on;plot(x,x*b + c,'Color',colr(count,:),'LineWidth',2);
          xlabel('target angle (deg)');ylabel('rotational error (deg)');title('angular error as a function of angle');
          legend(leg_input{i});grid on;
          % obtain figures to rearrange them on the screen later
          fig(indx) = gcf;
          fig(indx).Units = 'centimeters';
      end
  end
%   suptitle([trials(1).prs.subject]);
  % place figures on Right screen in vertical order
  fig = unique(fig);
  if any(strcmp(condnames,'s1'))
      count = 0;
      for i = 1:length(fig)
          figdim = fig(i).Position(3:4);
          fig(i).Position = [-.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
          count = count + 1;
          figure(fig(i));
          suptitle([trials(1).prs.subject]);
      end
  elseif any(strcmp(condnames,'bin1'))
      count = 0;
      for i = 1:length(fig)
          figdim = fig(i).Position(3:4);
          fig(i).Position = [-.5+(count*(figdim(1)-.5)) 2 figdim]; % 0: max left position , 2: max down position
          count = count + 1;
          figure(fig(i));
          suptitle([trials(1).prs.subject]);
      end
  else
      count = 0;
      for i = 1:length(fig)
          figdim = fig(i).Position(3:4);
          fig(i).Position = [-.5+(count*(figdim(1)-.5)) 2 figdim]; % 0: max left position , 2: max down position
          count = count + 1;
          figure(fig(i));
          suptitle([trials(1).prs.subject]);
      end
  end
  
end
