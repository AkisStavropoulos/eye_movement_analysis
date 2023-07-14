function [trials,behv_stats] = extractMonkeySession(monkeyname,sessiondate,sessionfold)



if strcmpi(sessionfold(1),'Z')
    recSpec = 'U-probe';
    if strcmpi(monkeyname,'Quigley')
        recSpec = 'Sim_recordings';
elseif strcmpi(monkeyname,'Ody') || strcmpi(monkeyname,'Bruno') 
        recSpec = 'Utah Array';
elseif strcmpi(monkeyname,'Schro') && strcmpi(sessiondate,'Feb 08 2018')
        recSpec = 'Utah Array';
    end
elseif strcmpi(sessionfold(1),'W')
    recSpec = 'Stimulation';
end

datafold = fullfile(sessionfold,monkeyname,recSpec,sessiondate,'neural data\Pre-processing X E');

% load data
fname = dir(fullfile(datafold,'m*s*.mat'));
cd(datafold) 
load(fname(1).name)

trials = trials_behv;

% rename v_max field
for i = 1:numel(trials)
    if isfield(trials(i).prs,'v_max')
        trials(i).prs.vmax = trials(i).prs.v_max;
        trials(i).prs = rmfield(trials(i).prs,'v_max');
    end
    if isfield(trials(i).prs,'w_max')
        trials(i).prs.wmax = trials(i).prs.w_max;
        trials(i).prs = rmfield(trials(i).prs,'w_max');

    end
end

fprintf('Monkey: %s\nSession: %s \nNtrials = %d\n',monkeyname,fname(1).name,numel(trials));
