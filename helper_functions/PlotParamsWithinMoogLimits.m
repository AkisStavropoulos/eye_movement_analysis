%% plot parameters within moog limitations
% for CoefSim script

h1 = figure(1);clf;h2 = figure(2);clf;h3 = figure(3);clf
if exist('T_ind') &&  exist('x_ind')
    count = 1;
    % divide subplots in columns of 5 plots
    r = rem(size(ind,1),5);
    sbpltcols = (size(ind,1)-r)/5;
    if r
        sbpltcols = sbpltcols + 1;
    end
    %
    for i = 1:size(ind,1)
        clear leg_input; clear leg_input_temp;
        for j = 1:size(ind,2)
            alpha = [bounded_params(ind{i,j}).jscoef];
            taus = [bounded_params(ind{i,j}).tau];
            vmaxs = [bounded_params(ind{i,j}).vmax];
            
            leg_input_temp(1,j) = {['x = ' num2str(x(i)) ', T = ' num2str(T(j)) ' s']};
            leg_input = {leg_input_temp{:}};
            
            if ~isempty(ind{i,j})
                figure(1);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(taus,vmaxs,'.');title(['x = ' num2str(x(i)) ' cm']);
                suptitle(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('tau (s)');ylabel('Vmax');
                
                figure(2);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(-dt./log(alpha),vmaxs,'.');title(['x = ' num2str(x(i)) ' cm']);
                suptitle(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('inputed tau');ylabel('Vmax');
                
                figure(3);grid on;
                subplot(5,sbpltcols,count);hold on;
                plot(alpha,vmaxs,'.');title(['x = ' num2str(x(i)) ' cm']);
                suptitle(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('Joystick coefficient');ylabel('Vmax');
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
        figure(1);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','northwest','Box','off')
        figure(2);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','northwest','Box','off')
        figure(3);subplot(5,sbpltcols,count);legend(leg_input{:},'Location','northwest','Box','off')

        count = count + 1;
    end
else
    for i = 1:length(ind)
        
        alpha = [bounded_params(ind{i}).jscoef];
        taus = [bounded_params(ind{i}).tau];
        vmaxs = [bounded_params(ind{i}).vmax];
        
        if exist('T_ind')
            leg_input(i) = {['T = ' num2str(T(i)) ' s']};
        elseif exist('x_ind')
            leg_input(i) = {['x = ' num2str(x(i)) ' s']};
        end
        
        figure(1);hold on;
        plot(taus,vmaxs,'.');legend(leg_input{:});
        title(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('tau (s)');ylabel('Vmax');
        
        figure(2);hold on;
        plot(-dt./log(alpha),vmaxs,'.');legend(leg_input{:});
        title(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('inputed tau');ylabel('Vmax');
        
        figure(3);hold on;
        plot(alpha,vmaxs,'.');legend(leg_input{:});
        title(['Allowed parameters for GIA error < ' num2str(thresh) ' m/s^2']);xlabel('Joystick coefficient');ylabel('Vmax');
    end
end
