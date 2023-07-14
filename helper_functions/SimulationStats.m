%% Statistics of simulations

%% trend of increasing stopping distance from target as a function of x (ff_dist)
 for i = 1:size(ind,1) % or 1:length(x)
stop_dist_offset(i) = mean(dist([ind{i,:}])) - x(i);
 end
 
 figure;plot(x,stop_dist_offset,'-o');xlabel('x (cm)');ylabel('stopping distance offset (cm)');
 title('switch based, stop distance from target');grid on;%axis equal;
 
 %%  trial duration change as a function of T
 % average across x's
  for j = 1:size(ind,2) % or 1:length(T)
trl_dur(j) = mean(tt([ind{:,j}])) - T(j);
 end
 
 figure;plot(T,trl_dur,'-o');xlabel('T (s)');ylabel('trial duration offset (s)');
 title('switch based trial duration added time');grid on;ylim([0 1]);%axis equal;
 
 % average for every x
 for j = 1:size(ind,2) % or 1:length(T)
     for i = 1:size(ind,1) % or length(x)
         trl_dur2(i,j) = mean(tt([ind{i,j}])) - T(j);
     end
 end
 
 figure;plot(T,trl_dur2,'-o');xlabel('T (s)');ylabel('trial duration offset (s)');
 title('switch based trial duration added time');grid on;ylim([0 1]);%axis equal;

 %% find x, T that give the biggest allowed range of taus
  % create a parameters struct
  clear params;
 count = 1;
 for i = 1:length(x)
     for j = 1:length(T)
         params(count).x = x(i);
         params(count).T = T(j);
         count = count + 1;
     end
 end
 % find the range for each parameters pair
lbtau = [];ubtau = [];tau_range = [];
count = 1;
 for i = 1:size(ind,1)
     for j = 1:size(ind,2)
         if ~isempty(ind{i,j})
             
         lbtau(i,j) = min([bounded_params([ind{i,j}]).tau_input]);
         ubtau(i,j) = max([bounded_params([ind{i,j}]).tau_input]);
         tau_range(i,j) = ubtau(i,j) - lbtau(i,j);
         
         params(count).lbtau = lbtau(i,j);
         params(count).ubtau = ubtau(i,j);
         params(count).tau_range = tau_range(i,j);
         end
         count = count + 1;
     end
 end
 % plot
 cmap = colormap(jet);
 colorsteps = ceil(length(cmap)/length(x));
 cmap = cmap(1:colorsteps:end,:);
 figure;
 for i = 1:size(tau_range,1)
     hold on;plot(tau_range(i,:),'-o','Color',cmap(i,:));
     leg{i} = ['x = ' num2str(x(i)) ' cm'];
 end
 ylabel('tau range (s)');xlabel('values of T (s)');title('tau range for increasing x values');legend(leg{:});
 %
 remov = [];
 for i = 1:length(params)
     if ~isempty(params(i).tau_range)
         fval(i) = (params(i).T - params(i).tau_range)/params(i).T; % (T-tau_range)/T
         params(i).fval = fval(i);
     else remov = [remov i];
     end
 end
 params([remov]) = [];
 fval([remov]) = [];
 [winner,indx] = min(fval);
 
 figure;plot(fval);
 
 disp(['And the winner is... x = ' num2str(params(indx).x) ' cm, T = ' num2str(params(indx).T) ...
     ' s, with tau range = ' num2str(params(indx).ubtau - params(indx).lbtau) ' !!!'])