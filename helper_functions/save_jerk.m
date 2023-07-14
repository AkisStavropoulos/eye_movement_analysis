function trials = save_jerk(trials,xjerk,yjerk,jerk,cumjerk)
%% save cum squared jerk in trials struct
if ~isempty(cumjerk)
    for i = 1:length(trials)
        trials(i).continuous.cumjerk = cumjerk{i};
    end
end

if ~isempty(xjerk)
    for i = 1:length(trials)
        trials(i).continuous.xjerk = xjerk{i};
    end
end

if ~isempty(yjerk)
    for i = 1:length(trials)
        trials(i).continuous.yjerk = yjerk{i};
    end
end

if ~isempty(jerk)
    for i = 1:length(trials)
        trials(i).continuous.jerk = jerk{i};
    end
end
