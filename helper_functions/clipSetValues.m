function x = clipSetValues(x,indx,val)

%% Set values to val
% x = clipSetValues(x,indx)
% indx can be an integer or logical

if nargin < 2
    indx = 1:numel(x);
end

x(indx) = val;