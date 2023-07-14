function ind = get_params_indices(x,T,mystruct)

%% divide struct into categories based on variable
% find indices for same variable (x or T)
% for the simulations
% input params struct or ff_range struct

% rows: x
% cols: T

if length(T)>1
    for j = 1:length(T)
        indx2 = find([mystruct(:).T] == T(j));
        T_ind(j) = {indx2};
    end
end
if length(x)>1
    for i = 1:length(x)
        indx3 = find([mystruct(:).x] == x(i));
        x_ind(i) = {indx3};
    end
end

if exist('T_ind') &&  exist('x_ind')
    for i = 1:length(x_ind)
        for j = 1:length(T_ind)
            ind{i,j} = intersect(x_ind{i},T_ind{j});
        end
    end
elseif exist('T_ind')
    ind = T_ind;
elseif exist('x_ind')
    ind = x_ind;
else ind = {1:length(mystruct)};
end
rmsimrow = [];
rmsimcol = [];
for i = 1:size(ind,1)
    for j = 1:size(ind,2)
    if isempty(ind{i,j})
        rmsimrow = [rmsimrow i];
        rmsimcol = [rmsimcol j];
    end
    end
end
% if rmsimrow
%     for i = 1:length(rmsimrow)
%     ind(rmsimrow(i),rmsimcol(i)) = {nan};
%     end
% end
