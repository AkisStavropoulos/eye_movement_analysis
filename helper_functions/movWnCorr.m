function [rho,pval,nx] = movWnCorr(varmin,varmax,stp,wn,bin_var,x,y)

% calculate correlation with sliding window
% also outputs the indices of data vector corresponding to each bin
% data is centered on the min edge

nx = varmin:stp:varmax;

rho = nan(numel(nx),1); pval = nan(numel(nx),1); 
for k = 1:numel(nx)

binindx =  bin_var >= nx(k) & bin_var < nx(k)+wn ;

x_tmp = x(binindx);
y_tmp = y(binindx);

if ~isempty(x_tmp)
    [rho(k),pval(k)] = nancorr(x_tmp(:),y_tmp(:));
end

% binindx{k} = find(binindx);

end