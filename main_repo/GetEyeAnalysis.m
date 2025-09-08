function [eye_fixation,eye_movement] = GetEyeAnalysis(subject,prs,params,upperbnd,BELIEF)
%%
% BELIEF: 1 to use subject's believed position
%         0 to use actual position
% upperbnd: compute upper bound of TTI based on subject's positional uncertainty
[poolindx,~] = get_poolindx(subject,params);
eye_fixation = [];  eye_movement = [];
for i = 1:length(subject) % [4 5 11] are bad
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        if upperbnd
            spatialstd = ComputeSpatialError(subject(i).trials(indx));
        else
            spatialstd = [];
        end
        % compute
        eye_fixation{i,s} = AnalyseEyefixation(subject(i).trials(indx),prs);
        eye_movement{i,s} = AnalyseEyemovement(eye_fixation{i,s},subject(i).trials(indx),spatialstd,prs,BELIEF);
    end
    disp(['......Subject = ' num2str(i)])
end
