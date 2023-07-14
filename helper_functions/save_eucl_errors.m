function trials = save_eucl_errors(trials,eucl_err)
%% save euclideian error
for i = 1:length(trials)
    trials(i).stats.eucl_err = eucl_err(i);
end