function [trials,delta,deltas] = comp_delta(trials,bins,history_depth)

%% By tau difference (delta)
% bins = 4; % number of tau difference bins
indx = [];
taus = [trials.prs];
taus = [taus.tau];
deltas = [0 diff(taus)];
k = history_depth;

tau_history = [];
delta_history = [];
for i = 4:length(trials)
    
    tau_history(i,:) = taus(i:-1:i-k); % most recent is the first element in each row, row = indx of that element
    delta_history(i,:) = deltas(i:-1:i-(k-1));
end

% delta_history = abs(delta_history);
deltas_sort = sort(deltas);
% divide based on cdf of tau differences
delta_cdf = cdf('Normal',deltas_sort,mean(deltas_sort),std(deltas_sort));
% figure;plot(deltas_sort,delta_cdf);grid on;

bin_size = max(delta_cdf)/(bins);

bin_ind = [];
count = 1;
for i = 1:bins-1 % if without -1, you'll get bins+1 number of bins (like dividers)
    bin_ind(i) = find(deltas_sort <= count*bin_size, 1,'last');
    count = count+1;
end

% % or set bins: [0 1/2] , [1/2 5/6] , [5/6 1]
% bin_cuts = [1/2 5/6];
% bin_ind = [];
% for j = 1:length(bin_cuts)
%     bin_ind(j) = find(delta_cdf <= bin_cuts(j),1,'last');
% end

delta_lims = [min(deltas_sort) deltas_sort(bin_ind) max(deltas_sort)];
% or set manually
delta_lims = [min(deltas_sort) -1.25 -0.25 0.25 1.25 max(deltas_sort)];
% bin deltas
for i = 1:length(deltas)
   bin_id = find(delta_lims <= deltas(i),1,'last'); 
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
%% save on trials struct the recent history of taus
% most recent is the first element!!!
k = history_depth;
for i = 1:length(trials)
    if i - k > 0
        tau_history = taus(i:-1:i-k);
        trials(i).stats.tau_history = tau_history;
    else
        if i > 1
            tau_history = taus(i:-1:1);
            trials(i).stats.tau_history = tau_history;
        else
            trials(i).stats.tau_history = [];
        end
    end
end
%% save on trials struct the recent history of deltas
% most recent is the first element!!!
k = history_depth;
for i = 1:length(trials)
    if i - k > 0
        delta_history = deltas(i:-1:i-(k-1));
        trials(i).stats.delta_history = delta_history;
    else
        if i > 1
            delta_history = deltas(i:-1:2);
            trials(i).stats.delta_history = delta_history;
        else
            trials(i).stats.delta_history = [];
        end
    end
end
%% save on trials struct the corresponding bins of deltas
hist_depth = length(trials(end).stats.delta_history);

for i = 2:length(trials)
    delta_history = trials(i).stats.delta_history;
    if i - hist_depth > 0
        for k = 1:hist_depth
            for n = 1:size(bin_lim,1)
                if delta_history(k) <= bin_lim(n,2) && delta_history(k) >= bin_lim(n,1)
                    bin(k) = n;
                end
            end
        end
    else
        for k = 1:length(trials(i).stats.delta_history)
            for n = 1:size(bin_lim,1)
                if delta_history(k) <= bin_lim(n,2) && delta_history(k) >= bin_lim(n,1)
                    bin(k) = n;
                end
            end
        end
    end
    trials(i).stats.hist_bins = bin;
end
trials(1).stats.hist_bins = [];
            

