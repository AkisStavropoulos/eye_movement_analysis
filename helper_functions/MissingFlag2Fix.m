%% MissingFlag2Fix
%% Missing MC file fix for flag2 field
for i = 1:length(subject)
    for j = 1:length(subject(i).trials)
        if numel(fieldnames(subject(i).trials(j).mc)) == 38
            subject(i).trials(j).mc.flag2 = [];
        end
    end
end
