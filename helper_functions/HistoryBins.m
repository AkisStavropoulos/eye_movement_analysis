%% by tau bins
taunames = fieldnames(tau);
taubins = tau;
count = 1;
for j = 1:length(taunames)
    if strcmp(taunames{j}(1:3),'lim') || strcmp(taunames{j}(1:3),'all')
        taubins = rmfield(taubins,taunames{j});
        if ~strcmp(taunames{j}(1:3),'all')
            lims(count,:) = tau.(taunames{j}) ;
            count = count + 1;
        end
    end
end
triplet = [];
doublet = [];
singlet = [];
history = [];
doub_bins = [];
trip_bins = [];
sing_bins = [];
for i = 3:length(trials)
    % current tau
    cur_tau = trials(i).prs.tau;
    cur_bin = intersect(find(lims(:,1) <= cur_tau),find( lims(:,2) >= cur_tau));
    % 1 back tau
    back1 = i-1;
    back1_tau = trials(back1).prs.tau;
    back1_bin = intersect(find(lims(:,1) <= back1_tau),find( lims(:,2) >= back1_tau));
    % 2 back tau
    back2 = i-2;
    back2_tau = trials(back2).prs.tau;
    back2_bin = intersect(find(lims(:,1) <= back2_tau),find( lims(:,2) >= back2_tau));
        
    if cur_bin == back1_bin && cur_bin == back2_bin
       triplet  = [triplet i];
       trip_bins(length(triplet),:) = [cur_bin back1_bin back2_bin];
    elseif (cur_bin == back1_bin && cur_bin ~= back2_bin)
        doublet = [doublet i];
        doub_bins(length(doublet),:) = [cur_bin back1_bin back2_bin];
    elseif cur_bin ~= back1_bin
        singlet = [singlet i];
        sing_bins(length(singlet),:) = [cur_bin back1_bin back2_bin];
    end

end  
% each vector contains the trial indices for each type of history combination
% each 3 element vector contains the history of tau bins for the corresponding index
history.triplet = triplet';
history.trip_bins = trip_bins;
history.doublet = doublet';
history.doub_bins = doub_bins;
history.singlet = singlet';
history.sing_bins = sing_bins;
% sanity check
fucked = [];
for i = 1:size(trip_bins,1)
if size(unique(trip_bins(i,:))) > 1
    fucked = [fucked i];
end

end
%% By tau difference (delta)
bins = 4; % number of tau difference bins
indx = [];
taus = [trials.prs];
taus = [taus.tau];
deltas = [0 diff(taus)];

tau_history = [];
delta_history = [];
for i = 4:length(trials)
    
    tau_history(i,:) = taus(i:-1:i-3); % most recent is the first element in each row, row = indx of that element
    delta_history(i,:) = deltas(i:-1:i-2);
end

delta_history = abs(delta_history);
deltas_sort = sort(abs(deltas));
% divide based on cdf of tau differences
delta_cdf = cdf('Exponential',deltas_sort,mean(deltas_sort));

bin_size = max(delta_cdf)/(bins);
% or set bins: [0 1/2] , [1/2 5/6] , [5/6 1]

bin_ind = [];
count = 1;
for i = 1:bins-1 % if without -1, you'll get bins+1 number of bins (like dividers)
    bin_ind(i) = find(deltas_sort <= count*bin_size, 1,'last');
    count = count+1;
end

delta_lims = [min(deltas_sort) deltas_sort(bin_ind) max(deltas_sort)];
% bin deltas
for i = 1:length(deltas)
   bin_id = find(delta_lims <= abs(deltas(i)),1,'last'); 
   deltas_id(i) = bin_id;
end
% create struct
delta = [];
for i = 1:bins
    delta = setfield(delta,['lim' num2str(i)],[]);
    delta = setfield(delta,['bin' num2str(i)],[]);
end

for i = 1:bins
    % limits
    bin_lim(i,:) = [delta_lims(i) delta_lims(i+1)];
    delta.(['lim' num2str(i)]) = bin_lim(i,:);
    % indices
    indx_delta = intersect(find(deltas >= bin_lim(i,1)),find(deltas < bin_lim(i,2)));
    if i == length(delta_lims)-1 && ~isempty(find(deltas == bin_lim(i,2),1))
        indx_delta = [indx_delta find(deltas == bin_lim(i,2))];
    end
    delta.(['bin' num2str(i)]) = indx_delta;
    indx_delta = [];
end
%% Calculate the mean error of each history bin and plot

% create a function that takes as inputs the vector containing the indices
% and the trias struct
