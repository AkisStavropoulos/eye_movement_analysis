%% NanMCFix
%% save empty MC files as nan instead of empty
for i = 1:length(subject)
    varnames = fieldnames(subject(i).trials(1).mc);
    for j = 1:length(subject(i).trials)
        if isempty(subject(i).trials(j).mc.timestamp)
            for n = 1:length(varnames)
                subject(i).trials(j).mc.(varnames{n}) = nan;
            end
        end
    end
end