function trials = save_acc(trials,xacc,yacc,acc,cumacc)
%% save cum squared jerk in trials struct
if ~isempty(cumacc)
    for i = 1:length(trials)
        trials(i).continuous.cumacc = cumacc{i};
    end
end

if ~isempty(xacc)
    for i = 1:length(trials)
        trials(i).continuous.xacc = xacc{i};
    end
end

if ~isempty(yacc)
    for i = 1:length(trials)
        trials(i).continuous.yacc = yacc{i};
    end
end

if ~isempty(acc)
    for i = 1:length(trials)
        trials(i).continuous.acc = acc{i};
    end
end