function [indx] = plot_ranges(mystruct,ind)
% input params struct


for i = 1:size(ind,1)
    indx{i} = [ind{i,:}];
end

for i = 1:length(indx)
    if ~isempty(indx{i})
        xx(i) = mystruct(indx{i}(1)).x;
        TT{i} = [mystruct(indx{i}).T];
        
        tau_range{i} = [mystruct(indx{i}).tau_range];
        lbtau{i} = [mystruct(indx{i}).lbtau];
        ubtau{i} = [mystruct(indx{i}).ubtau];
        
        ff_range{i} = [mystruct(indx{i}).ff_range];
        lbff{i} = [mystruct(indx{i}).lbff];
        ubff{i} = [mystruct(indx{i}).ubff];
        
        leg_input_temp(1,i) = {['x = ' num2str(xx(i)) ' m']};
        leg_input = {leg_input_temp{:}};
        
    end
end

% divide subplots in columns of 3 plots
r = rem(length(indx),3);
sbpltcols = (length(indx)-r)/3;
if r
    sbpltcols = sbpltcols + 1;
end

for n = 1:4
    h(n) = figure;
end
% plot
for i = 1:length(xx)
    
    figure(h(1));hold on;plot(TT{i},tau_range{i},'-o');legend(leg_input{:});
    xlabel('T (s)');ylabel('\tau range');grid on;title('\tau range as a function of T');ylim([0 6]);xlim([3 15]);hold off;
    
%     figure(h(2));hold on;plot(TT{i},ff_range{i},'-o');legend(leg_input{:});
%     xlabel('T (s)');ylabel('firefly range');grid on;title('firefly range as a function of T');ylim([0 1100]);xlim([3 15]);hold off;
    
    
    
    figure(h(3));subplot(3,sbpltcols,i);hold on;errorbar(TT{i},lbtau{i},[],ubtau{i}-lbtau{i},'o');plot(TT{i},tau_range{i},'-o');
    title(['x = ' num2str(xx(i)) ' m']);legend('LB and UB \tau','absolute range');
    xlabel('T (s)');ylabel('\tau range');grid on;suptitle('\tau range as a function of T');ylim([0 6]);xlim([3 15]);hold off;
    
%     figure(h(4));subplot(3,sbpltcols,i);hold on;errorbar(TT{i},lbff{i},[],ubff{i}-lbff{i},'o');plot(TT{i},ff_range{i},'-o');
%     title(['x = ' num2str(xx(i)) ' m']);legend('LB and UB distance','absolute range');
%     xlabel('T (s)');ylabel('firefly distance range');grid on;suptitle('firefly range as a function of T');ylim([0 1100]);xlim([3 15]);hold off;
end
