%% check switch time
% check simulations for s = tau.*log((1 + exp(T./tau))./2);
% for CoefSim script
tic;
count = 1;
sw = switchtime([bounded_params.tau_input],[bounded_params.T]);
Tpulse = sw;
% Tbrake = T - sw;

for i = 1:length(bounded_params)
    [vmax(i), ~,~,~,~, dist(i), vel{i},tt(i)] = coef2vmax2(bounded_params(i).jscoef, Tpulse(i), bounded_params(i).T-Tpulse(i), [], bounded_params(i).x, bounded_params(i).T, 0);
end
toc;
h7 = figure(7);clf;h8 = figure(8);clf;h9 = figure(9);clf
if exist('T_ind') &&  exist('x_ind')
    count = 1;
    % divide subplots in columns of 5 plots
    r = rem(size(ind,1),5);
    sbpltcols = (size(ind,1)-r)/5;
    if r
        sbpltcols = sbpltcols + 1;
    end
    
    for i = 1:size(ind,1)
        clear leg_input; clear leg_input_temp;
        for j = 1:size(ind,2)
            alpha = [bounded_params(ind{i,j}).jscoef];
            taus = [bounded_params(ind{i,j}).tau];
            vmaxs = [bounded_params(ind{i,j}).vmax];
            
            leg_input_temp(1,j) = {['x = ' num2str(x(i)) ', T = ' num2str(T(j)) ' s']};
            leg_input = {leg_input_temp{:}};
            
            if ~isempty(ind{i,j})
                
                figure(7);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(taus,tt(ind{i,j}),'.');title(['x = ' num2str(x(i)) ' cm']);
                xlabel('inputed tau (s)');ylabel('trial duration (s)');hline(T(j));
                suptitle('trial duration for different taus, for s switch time');
                
                figure(8);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(taus,dist(ind{i,j}),'.');title(['x = ' num2str(x(i)) ' cm']);
                xlabel('inputed tau (s)');ylabel('travel distance (cm)');
                suptitle('travel distance for different taus, for s switch time');
                hline(x(i));%ylim([x(i)-20 x(i)+20]);
                
                figure(9);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(tt(ind{i,j}),dist(ind{i,j}),'.');title(['x = ' num2str(x(i)) ' cm']);
                ylabel('distance (cm)');xlabel('travel duration (s)');
                suptitle('distance over duration for s switch time');
                hline(x(i));%ylim([x(i)-20 x(i)+20]);
            else leg_input_temp(j) = [];
            end
        end
        N = [];
        for n = 1:length(leg_input)
            if isempty(leg_input{n})
                N = [N n];
            end
        end
        leg_input(N) = [];
        figure(7);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','southwest','Box','off')
        figure(8);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','southwest','Box','off')
        figure(9);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','southwest','Box','off')
        
        count = count + 1;
    end
else
    h7 = figure(7);clf;h8 = figure(8);clf;h9 = figure(9);clf
    for i = 1:length(ind)
        alpha = [bounded_params(ind{i}).jscoef];
        taus = [bounded_params(ind{i}).tau_input];
        vmaxs = [bounded_params(ind{i}).vmax];
        
        if exist('T_ind')
            leg_input(i) = {['T = ' num2str(T(i)) ' s']};
            suptitle_name = ['x = ' num2str(x) ' cm'];
            figure(7);suptitle(suptitle_name);
            figure(8);suptitle(suptitle_name);
            figure(9);suptitle(suptitle_name);
        elseif exist('x_ind')
            leg_input(i) = {['x = ' num2str(x(i)) ' cm']};
            suptitle_name = ['T = ' num2str(T) ' s'];
            figure(7);suptitle(suptitle_name);
            figure(8);suptitle(suptitle_name);
            figure(9);suptitle(suptitle_name);
        else
            leg_input(i) = {['x = ' num2str(x(i)) ' cm' ', T = ' num2str(T(i)) ' s']};
            suptitle_name = ['x = ' num2str(x) ' cm' ', T = ' num2str(T) ' s'];
            figure(7);suptitle(suptitle_name);
            figure(8);suptitle(suptitle_name);
            figure(9);suptitle(suptitle_name);

        end
        
        figure(7);hold on;
        plot(taus,tt(ind{i}),'*');xlabel('inputed tau (s)');ylabel('trial duration (s)');hline(T);ylim([min(T)-1 max(T)+1]);
        title('trial duration for different taus, for s switch');legend(leg_input{:});
        
        figure(8);hold on;
        plot(taus,dist(ind{i}),'.');xlabel('inputed tau (s)');ylabel('travel distance (cm)');
        title('travel distance for different taus, for s switch time');legend(leg_input{:});
        hline(bounded_params(ind{i}(1)).x);%ylim([bounded_params(i).x-20 bounded_params(i).x+20]);
        
        figure(9);hold on;
        plot(tt(ind{i}),dist(ind{i}),'.');ylabel('distance (cm)');xlabel('travel duration (s)');
        title('distance over duration for s switch time');vline(bounded_params(ind{i}(1)).T);
        legend(leg_input{:});hline(bounded_params(ind{i}(1)).x);%ylim([bounded_params(i).x-20 bounded_params(i).x+20]);
    end
end
