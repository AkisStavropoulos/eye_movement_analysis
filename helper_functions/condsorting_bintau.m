function [stimtype,tau,stimtau] = condsorting_bintau(trials,bins,limtype,tau_cdf,taus_x)
%% Create struct for condition sorting
% group taus
% choose trials based on conditions
% limtype: string 'hard' or 'flex' ONLY
while ~exist('hard')
    if strcmp(limtype,'hard')
        hard = 1;
        flex = 0;
    elseif strcmp(limtype,'flex')
        flex = 1;
        hard = 0;
    else
        limtype = input('Wrong limit type. Enter string ''hard'' or ''flex'' ONLY: ')    ;
    end
end

params = [trials.prs];
taus = [params.tau];
stimulus = [params.stimtype];

s = unique(stimulus);
%% Define tau limits
taus_sort = taus_x;
if flex
    % flexible limits
    smallest_tau = min(taus_sort);
    biggest_tau = max(taus_sort);
    range_taus = biggest_tau - smallest_tau;
    bin_size = range_taus/bins;
    
    count = 1;
    for i = 1:bins
        bin_ind(i) = find(taus_sort <= smallest_tau + count*bin_size, 1,'last');
        count = count+1;
    end
    tau_lims = [min(taus_sort) taus_sort(bin_ind) max(taus_sort)];
elseif hard
    % hard limits, by binning the cdf of all taus obtained in the experiment
    smallest_tau = min(taus_sort);
    biggest_tau = max(taus_sort);
    
    pctiles = linspace(min(tau_cdf),max(tau_cdf),bins+1);
    
    bin_ind = [];
    for i = 1:bins-1 % if without -1, you'll get bins+1 number of bins (like dividers)
        bin_ind(i) = find(tau_cdf <= pctiles(i+1),1,'last');    
    end   
        
    tau_lims = [smallest_tau taus_sort(bin_ind)' biggest_tau];
end
%% create structs for conditions
stimtype = [];
for i = 1:length(s)
    stimtype = setfield(stimtype,['s' num2str(s(i))],[]);
end
tau = [];
tau = setfield(tau,'lims',[]);
tau = setfield(tau,'bins',[]);
for i = 1:bins
    tau.lims = setfield(tau.lims,['lim' num2str(i)],[]);
    tau.bins = setfield(tau.bins,['bin' num2str(i)],[]);
end
%% fill the structs
% stimtype
for j = 1:length(trials)
    indx_s = find(s == trials(j).prs.stimtype);
    if indx_s
        stimtype.(['s' num2str(s(indx_s))]) = [stimtype.(['s' num2str(s(indx_s))]) j];
    end
end
% tau
for i = 1:bins
    % limits
    bin_lim(i,:) = [tau_lims(i) tau_lims(i+1)];
    tau.lims.(['lim' num2str(i)]) = bin_lim(i,:);
    % indices
    indx_tau = intersect(find(taus >= bin_lim(i,1)),find(taus < bin_lim(i,2)));
    if i == length(bins) && ~isempty(find(taus == bin_lim(i,2)))
        indx_tau = [indx_tau find(taus == bin_lim(i,2))];
    end
    tau.bins.(['bin' num2str(i)]) = indx_tau;
    indx_tau = [];
end
stnames = fieldnames(stimtype);

taubinnames = fieldnames(tau.bins);
taulimnames = fieldnames(tau.lims);
%% find the intersection of the conditions
for i = 1:length(stnames)
    
    for j = 1:length(taubinnames)
        
        indx = intersect(stimtype.(stnames{i}),tau.bins.(taubinnames{j}));
        stimtau.bins.([stnames{i} taubinnames{j}]) = indx;
        stimtau.lims.(taulimnames{j}) = tau.lims.(taulimnames{j});
        clear indx;
    end
end
    
stimtype.all = 1:length(trials);
tau.all = 1:length(trials);     
stimtau.all = 1:length(trials);