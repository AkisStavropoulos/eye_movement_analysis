function plot_ff_ranges(mystruct,spec)
% input ff_range struct
% spec: specifywhich simulation you are plotting

count = 1;
difx = 0; difT = 0;
ubff_temp = [];
lbff_temp = [];
tau_temp = [];

for k = 1:length(mystruct)
    xx(count) = mystruct(k).x;
    TT(count) = mystruct(k).T;
        
    ubff_temp = [ubff_temp mystruct(k).ubff];
    lbff_temp = [lbff_temp mystruct(k).lbff];
    
    tau_temp = [tau_temp mystruct(k).tau];
    
    
    if k < length(mystruct)
        difx = mystruct(k).x - mystruct(k+1).x;
        difT = mystruct(k).T - mystruct(k+1).T;
    elseif k == length(mystruct)
        difx = -1;
        difT = -1;
    end
    if difx || difT % to save separately the range for each tau, and then find intersection across all taus
        
        ff_range{count} = [ubff_temp - lbff_temp];
        tau{count} = tau_temp;
        
        ubff{count} = ubff_temp;
        lbff{count} = lbff_temp;
        
        leg_input_temp(1,count) = {['T = ' num2str(TT(count)) ' s']};
        leg_input = {leg_input_temp{:}};
        
        tau_temp = [];
        ubff_temp = [];
        lbff_temp = [];
        count = count + 1;
    end
end


% plot
unique_x = unique(xx);
for i = 1:length(unique_x)
    ind{i} = find(xx == unique_x(i));
end
for n = 1:2*length(ind)
    h(n) = figure;
end
m = 0;
for j = 1:length(ind)
    
    % divide subplots in columns of 3 plots
    r = rem(length(ind{j}),3);
    sbpltcols = (length(ind{j})-r)/3;
    if r
        sbpltcols = sbpltcols + 1;
    end
    for i = 1:length(ind{j})
        
        figure(h(1+m));hold on;plot(tau{ind{j}(i)},ff_range{ind{j}(i)},'-o');legend(leg_input{ind{j}});
        xlabel('\tau (s)');ylabel('firefly range (cm)');grid on;
        title([spec ': firefly range as a function of \tau for x = ' num2str(xx(ind{j}(i))) ' cm']);ylim([0 max([ubff{:}])+200]);xlim([0 max([tau{:}])+2]);hold off;
        
        
        figure(h(2+m));subplot(3,sbpltcols,i);hold on;errorbar(tau{ind{j}(i)},lbff{ind{j}(i)},[],ff_range{ind{j}(i)},'o');
        plot(tau{ind{j}(i)},ff_range{ind{j}(i)},'-o');
        title(['x = ' num2str(xx(ind{j}(i))) ' m, T = ' num2str(TT(ind{j}(i))) ' s']);legend('LB and UB distance','absolute range');
        xlabel('\tau (s)');ylabel('firefly range (cm)');grid on;
        suptitle([spec ': firefly range as a function of \tau']);ylim([0 max([ubff{:}])+200]);xlim([0 max([tau{:}])+2]);hold off;
    end
    m = m + 2;    
end