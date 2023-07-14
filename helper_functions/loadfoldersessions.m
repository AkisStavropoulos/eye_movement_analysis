function trials = loadfoldersessions(foldername)

%% Load one session of single firefly with motion cuing
if ~isempty(foldername)
    addpath(genpath(foldername));
end  
default_prs;
trials = AddTrials2Behaviour(prs); % includes some fixes
