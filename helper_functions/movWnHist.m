function [y,x,bin_indx] = movWnHist(tmin,tmax,stp,wn,data)

% calculate histogram with sliding window
% also outputs the indices of data vector corresponding to each bin

nx = tmin:stp:tmax;

y = []; x = []; bin_indx = [];
for k = 1:numel(nx)-1

x(k) = nanmean( [ nx(k) nx(k)+wn ] );

bin_data =  data >= nx(k) & data < nx(k)+wn ;

y(k) = sum(bin_data);

bin_indx{k} = find(bin_data);

end
