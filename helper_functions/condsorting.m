function [stimtype,tau,stimtau] = condsorting(trials)
%% Create struct for condition sorting
% Here, find consecutive trials with same tau for adaptation analysis
% choose trials based on conditions
params = [trials.prs];
taus = [params.tau];
stimulus = [params.stimtype];

s = unique(stimulus);
tau_temp = unique(taus);

% create structs for conditions
stimtype = [];
for i = 1:length(s)
    stimtype = setfield(stimtype,['s' num2str(s(i))],[]);
end


tau = [];
for i = 1:length(tau_temp)
    co = num2str(tau_temp(i));
    if length(co)>1
        if strcmp(co(2),'.')
            co(2) = '_';
        end
        
    end
    tau = setfield(tau,['tau' co],[]);
end

% fill the structs
for j = 1:length(trials)
    % stimtype
    indx_s = find(s == trials(j).prs.stimtype);
    if indx_s
    stimtype.(['s' num2str(s(indx_s))]) = [stimtype.(['s' num2str(s(indx_s))]) j];
    end
    % tau
    indx_tau = find(tau_temp == trials(j).prs.tau);
    if indx_tau
        co = num2str(tau_temp(indx_tau));
        if length(co)>1
            if strcmp(co(2),'.')
                co(2) = '_';
            end
        end
        tau.(['tau' co]) = [tau.(['tau' co]) j];
    end
    
end

% % fix for lack of data
% if  exist('jscoef.a95')
%     if length(tau.a0) < 200
%         
%         a0_95 = [tau.a0 tau.a95];
%         tau = rmfield(tau,'a0');
%         tau = rmfield(tau,'a95');
%         tau = setfield(tau,'a0_95',a0_95);
%         tau = orderfields(tau,{'a0_95' 'a975' 'a99'});
%         
%     elseif isempty(tau.a95)
%         
%         tau = rmfield(tau,'a95');
%         
%     end
% end

%
stnames = fieldnames(stimtype);
taunames = fieldnames(tau);

% find the intersection of the conditions
for i = 1:length(stnames)
    
    for j = 1:length(taunames)
        
        indx = intersect(stimtype.(stnames{i}),tau.(taunames{j}));
        stimtau.([stnames{i} taunames{j}]) = indx;
        clear indx;
    end
end
    
stimtype.all = 1:length(trials);
tau.all = 1:length(trials);     
stimtau.all = 1:length(trials);