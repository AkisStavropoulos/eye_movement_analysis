function [trl,mc] = AddJSData(data,t)
% organize the Motion Cueing variables that are saved from moogdots
% sampled at 60 Hz !!!
% Note that VR position is NOT subject position. It's the position of the SCREEN
% For SUBJECT's position, use SMR data
% If only one input provided , generates only MC variables
%% Put everything in a struct and extract variable names
if ~isempty(data)
    mc.timestamp = data(:,1);
    mc.timestamp = (mc.timestamp - mc.timestamp(1));
    
    mc.JS_X_Raw = data(:,2); % forward JS input scaled by vmax, looks converted
    mc.JS_Yaw_Raw = data(:,3); % rotational JS input scaled by wmax, looks converted
    varnames = fieldnames(mc);
    
    % interpolate missing timepoints
    ts = 0:1/60:mc.timestamp(end);
    for i=2:length(varnames)
        mc.(varnames{i}) = interp1(mc.timestamp,mc.(varnames{i}),ts,'previous');
    end
    mc.timestamp = interp1(mc.timestamp,mc.timestamp,ts);
    mc = structfun(@(x) x', mc, 'un', 0); % hold on;plot(mc.timestamp,mc.JS_X_Raw,'.')

    if ~isempty(t)
        %% extract trials 
        % dt = dt*prs.factor_downsample;
        for j=1:length(t.end)
            timeindx = mc.timestamp > t.beg(j) & mc.timestamp < t.end(j);
            trl(j).mc = structfun(@(x) x(timeindx), mc,'un',0);
        end
        
        emptytrl = [];
        for n = 1:length(trl)
            if ~isempty(trl(n).mc.timestamp)
                % start time at t = 0 for every trial
                trl(n).mc.timestamp =  trl(n).mc.timestamp - trl(n).mc.timestamp(1);
            else
                emptytrl = [emptytrl n];
            end
        end
        if emptytrl
            disp(['empty MC trials from trial ' num2str(min(emptytrl)) ' to trial ' num2str(max(emptytrl)) ' !!'])
        end
    end
    
else
    % if MC_Variables doesn't exist or empty
    mc.timestamp = nan;
    
    mc.JS_X_Raw = nan; % forward JS input scaled by vmax, looks converted
    mc.JS_Yaw_Raw = nan; % rotational JS input scaled by wmax, looks converted
    
    varnames = fieldnames(mc);
    
    for j=1:length(t.end)
        for i=1:length(varnames)
            trl(j).mc.(varnames{i}) = mc.(varnames{i});
            %             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
        end
    end
    
    disp('MC file is missing for this block !!')
    
end









