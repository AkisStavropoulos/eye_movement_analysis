function x = nanify(x,indx)

%% Set values to NaN
% x = nanify(x,indx)
% indx can be an integer or logical

if nargin < 2
    indx = 1:numel(x);
end

x(indx) = nan;